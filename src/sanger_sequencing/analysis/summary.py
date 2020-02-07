# Copyright (c) 2018-2020 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""Functions that drive the main Sanger sequence analysis workflow."""


import logging
from typing import List, Tuple

from Bio.Data.CodonTable import ambiguous_dna_by_name
from Bio.SeqRecord import SeqRecord
from numpy import isnan, nanmean
from pandas import DataFrame, concat

from ..config import Configuration
from ..model import (
    ConflictReportInternal,
    ConflictTypeEnum,
    QualityEnum,
    SampleReportInternal,
)


__all__ = ("summarize_plasmid_conflicts", "concatenate_sample_reports")


logger = logging.getLogger(__name__)

CODON_TABLE = ambiguous_dna_by_name["Standard"].forward_table
START_CODONS = frozenset(ambiguous_dna_by_name["Standard"].start_codons)
STOP_CODONS = frozenset(ambiguous_dna_by_name["Standard"].stop_codons)


def concatenate_sample_reports(reports: List[SampleReportInternal]) -> DataFrame:
    """Concatenate many data frames into one."""
    data = []
    for sample in reports:
        df = sample.details
        if df is None:
            continue
        df["sample"] = sample.id
        df["primer"] = sample.primer
        data.append(df)
    if len(data) == 0:
        return DataFrame(
            columns=["plasmid_pos", "sample_pos", "snp", "plasmid_chr", "sample_chr",]
        )
    else:
        return concat(data, ignore_index=True, copy=False)


def determine_type(row) -> ConflictTypeEnum:
    if isnan(row.plasmid_pos) and isnan(row.sample_pos):
        msg = (
            "Detected a gap in the alignment. Only expecting single "
            "position conflicts."
        )
        logger.error(msg)
        raise ValueError(msg)
    elif isnan(row.plasmid_pos):
        return ConflictTypeEnum.INSERTION
    elif isnan(row.sample_pos):
        return ConflictTypeEnum.DELETION
    else:
        return ConflictTypeEnum.CHANGE


def determine_quality(region, threshold: float) -> QualityEnum:
    if nanmean(region["quality"]) >= threshold:
        return QualityEnum.HIGH
    else:
        return QualityEnum.LOW


def confirm_conflict(conflict_type, row, cover, threshold) -> Tuple[int, int]:
    num_confirmed = 0
    num_invalidated = 0
    for sample_id, sub in cover.groupby("sample", as_index=False, sort=False):
        if determine_quality(sub, threshold) == "low":
            logger.debug(
                "Ignoring low quality sample region for conflict confirmation."
            )
            continue
        # Index 1 should correspond to the mid-point and thus to the
        # conflict site. However, due to sequence read trimming we may hit
        # the end of the sequence alignment.
        if len(sub) < 3:
            logger.debug(
                "Ignoring incomplete sample region (less than three positions)."
            )
            continue
        cmp = sub.iloc[1]
        if not cmp.snp:
            logger.debug("Not a conflict site.")
            num_invalidated += 1
            continue
        if determine_type(cmp) != conflict_type:
            logger.debug("Different type of conflict site.")
            num_invalidated += 1
            continue
        if (cmp.sample_chr == row.sample_chr) and (cmp.plasmid_chr == row.plasmid_chr):
            num_confirmed += 1
        else:
            num_invalidated += 1
    return num_confirmed, num_invalidated


def determine_effects(
    row, plasmid, previous, following
) -> Tuple[List[Tuple[str, List[str]]], List[str]]:
    """
    Post-process conflicts and categorize them.

    Parameters
    ----------
    row
    plasmid
    previous
    following

    Returns
    -------
    list, list

    """
    features = []
    effects = []
    for feat in plasmid.features:
        # Process only features if the conflict position overlaps.
        if not (
            previous <= feat.location.end.position
            and feat.location.start.position <= following
        ):
            continue
        features.append((feat.type, feat.qualifiers.get("label", [])))
        if feat.type == "CDS":
            # Potential frame shift (usually rather a sequencing error).
            if isnan(row.plasmid_pos) or isnan(row.sample_pos):
                effects.append("Frame shift")
                continue
            # Does the new codon cause an amino acid change?
            seq = feat.extract(plasmid)
            # Translate plasmid position into feature position.
            feat_pos = int(row.plasmid_pos) - feat.location.start.position
            feat_pos -= 1  # Transform to 0-indexing.
            codon_pos = feat_pos % 3
            codon = list(seq[(feat_pos - codon_pos) : (feat_pos - codon_pos + 3)])
            if len(codon) < 3:
                logger.error("SNP at the beginning or end of CDS. " "Unknown effect.")
                effects.append("Unknown")
                continue
            cdn = "".join(codon)
            try:
                plasmid_aa = CODON_TABLE[cdn]
            except (KeyError, ValueError):
                if cdn in START_CODONS:
                    logger.warning("Start codon hit on plasmid.")
                    plasmid_aa = "START"
                elif cdn in STOP_CODONS:
                    logger.warning("Stop codon hit on plasmid.")
                    plasmid_aa = "STOP"
                else:
                    logger.error("Unknown codon '%s' on plasmid.", cdn)
                    plasmid_aa = "UNKNOWN"
            codon[codon_pos] = row.sample_chr
            cdn = "".join(codon)
            try:
                sample_aa = CODON_TABLE[cdn]
            except (KeyError, ValueError):
                if cdn in START_CODONS:
                    logger.warning("Change to start codon on sample.")
                    sample_aa = "START"
                elif cdn in STOP_CODONS:
                    logger.warning("Change to stop codon on sample.")
                    sample_aa = "STOP"
                else:
                    logger.error("Unknown codon '%s' on sample.", cdn)
                    sample_aa = "UNKNOWN"
            effects.append(f"{plasmid_aa} -> {sample_aa}")
    return features, effects


def summarize_plasmid_conflicts(
    sample: DataFrame, total: DataFrame, plasmid: SeqRecord
) -> List[ConflictReportInternal]:
    """
    Add useful information on sequence conflicts and their surroundings.

    There are three principle types of conflicts compared to a reference
    sequence:
    1. insertions
    2. deletions
    3. changes

    Depending on the quality in the neighborhood and information from other
    samples, each conflict may be (highly) likely, unresolved (i.e.,
    an unclear situation), or resolved (invalidated by other samples).

    """
    config = Configuration()
    conflicts = []
    if sample is None:
        return conflicts
    # Show what happens around a mismatch location on all samples.
    logger.info("Assessing %d conflicts.", total["snp"].sum())
    for row in sample.loc[sample["snp"], :].itertuples():
        conflict = ConflictReportInternal(
            plasmid_position=None if isnan(row.plasmid_pos) else int(row.plasmid_pos),
            sample_position=None if isnan(row.sample_pos) else int(row.sample_pos),
            plasmid_character=row.plasmid_chr,
            sample_character=row.sample_chr,
        )
        # The assumption here is that there will only ever be single position
        # conflicts and we can check the sequence before and after for more
        # information. This probably holds for a good read. Bad reads should
        # be repeated and not analyzed automatically.
        region = total.iloc[[row.Index - 1, row.Index, row.Index + 1], :]
        # Determine the kind of conflict.
        try:
            conflict.type = determine_type(row)
        except ValueError:
            continue
        # Determine the quality of the region.
        conflict.surrounding_quality = determine_quality(region, config.threshold)
        # Check for more information on other samples.
        # Due to a potential gap we take the plasmid index position before and
        # check rows in the total in-between that index and index + 2 which
        # should cover the gap.
        index = total[
            (total["sample"] != row.sample)
            & (total["plasmid_pos"] == region["plasmid_pos"].min())
        ].index
        index = [j for i in index for j in range(i, i + 3)]
        cover = total.loc[index, :]
        conflict.num_confirmed, conflict.num_invalidated = confirm_conflict(
            conflict.type, row, cover, config.threshold
        )
        # Add feature data.
        conflict.features_hit, conflict.effect = determine_effects(
            row, plasmid, region["plasmid_pos"].min(), region["plasmid_pos"].max(),
        )
        conflicts.append(conflict)
    return conflicts
