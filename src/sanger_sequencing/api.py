# -*- coding: utf-8 -*-

# Copyright (c) 2018 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Provide a high-level interface for the Sanger sequence analysis workflow.

Functionality to retrieve identifiers and sequence records from repositories
can be imported from the ``sanger_sequencing.clients`` subpackage.
"""

import logging
from typing import Dict, List

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from . import analysis, validation
from .helpers import log_errors


__all__ = ("sanger_verification", "plasmid_report", "sample_report")

LOGGER = logging.getLogger(__name__)


def sanger_verification(template: DataFrame,
                        plasmids: Dict[str, SeqRecord],
                        samples: Dict[str, SeqRecord]
                        ) -> List[Dict]:
    """
    Perform a complete Sanger verification for many plasmids and sample reads.

    Parameters
    ----------
    template : pandas.DataFrame
        A template table with three columns: plasmid, primer, sample which
        are all identifiers.
    plasmids : dict
        A mapping from plasmid identifiers to sequence records.
    samples : dict
        A mapping from sample identifiers to sequence records.

    Returns
    -------
    list
        A list of dictionaries that are each a plasmid report.

    Raises
    ------
    AssertionError
        All function arguments are extensively validated and may raise errors.

    See Also
    --------
    plasmid_report

    """
    LOGGER.info("Validate template.")
    errors = validation.validate_template(template)
    if len(errors) > 0:
        log_errors(errors)
        raise AssertionError("Invalid analysis template.")
    LOGGER.info("Validate plasmids.")
    for plasmid in plasmids.values():
        validation.validate_plasmid(plasmid, [])
    LOGGER.info("Validate samples.")
    for sample in samples.values():
        validation.validate_sample(sample)
    template = validation.drop_missing_records(template, plasmids, samples)
    LOGGER.info("Generate reports.")
    return [plasmid_report(plasmid_id, plasmids[plasmid_id], sub, samples)
            for plasmid_id, sub in template.groupby(
            "plasmid", as_index=False, sort=False)]


def plasmid_report(plasmid_id: str,
                   sequence: SeqRecord,
                   template: DataFrame,
                   samples: Dict[str, SeqRecord]) -> Dict:
    """
    Create an analysis report for a single plasmid and one or more reads.

    The plasmid report contains detailed reports on each sample read
    performed. It then evaluates each sequence alignment conflict using
    information from other samples as appropriate.

    Parameters
    ----------
    plasmid_id : str
        The plasmid identifier.
    sequence : Bio.SeqRecord.SeqRecord
        The plasmid's sequence record.
    template : pandas.DataFrame
        A part of the template table concerning this plasmid only.
    samples : dict
        A mapping from sample identifiers to sequence records.

    Returns
    -------
    dict
        An individual plasmid report.

    """
    LOGGER.info("Analyze plasmid '%s'.", plasmid_id)
    report = {
        "id": plasmid_id,
        "name": sequence.name,
        "samples": [
            sample_report(row.sample, samples[row.sample], row.primer,
                          plasmid_id, sequence)
            for row in template.itertuples(index=False)]
    }
    # Post-process reports in order to classify conflicts.
    LOGGER.debug("Concatenate the detailed sample reports.")
    total = analysis.concatenate_sample_reports(report["samples"])
    for rep in report["samples"]:
        rep["conflicts"] = analysis.summarize_plasmid_conflicts(
            rep["details"], total, sequence)
    return report


def sample_report(sample_id: str,
                  sample_sequence: SeqRecord,
                  primer_id: str,
                  plasmid_id: str,
                  plasmid_sequence: SeqRecord) -> Dict:
    """
    Create an analysis report for a single sample read.

    Parameters
    ----------
    sample_id : str
        The sample identifier.
    sample_sequence :  Bio.SeqRecord.SeqRecord
        The sample's sequence record.
    primer_id : str
        The primer identifier.
    plasmid_id : str
        The plasmid identifier.
    plasmid_sequence : Bio.SeqRecord.SeqRecord
        The plasmid's sequence record.

    Returns
    -------
    dict
        An individual sample report.

    """
    LOGGER.info("Analyze sample '%s'.", sample_id)
    report = {
        "id": sample_id,
        "primer": primer_id,
        "readLength": len(sample_sequence),
        "errors": []
    }
    # Convert to base `float` for JSON compatibility.
    try:
        start, trimmed_seq, quality_scores, end, median = analysis.trim_sample(
            sample_sequence)
    except ValueError as err:
        report["errors"].append(str(err))
        return report
    report["medianQuality"] = median
    report["cutBeginning"] = int(start)
    report["cutEnd"] = int(end)
    align = analysis.emboss_alignment(
        sample_id, trimmed_seq, plasmid_id, plasmid_sequence)
    report["details"] = analysis.alignment_to_table(
        align, quality_scores, start)
    return report
