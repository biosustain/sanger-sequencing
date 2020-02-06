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
import typing
from pathlib import Path

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from . import analysis, validation
from .config import Configuration
from .helpers import log_errors
from .reports import PlasmidReport, SampleReport, SangerReport


__all__ = ("sanger_report", "plasmid_report", "sample_report")


logger = logging.getLogger(__name__)


def sanger_report(
    template: DataFrame,
    plasmids: typing.Dict[str, SeqRecord],
    samples: typing.Dict[str, SeqRecord],
    threshold: typing.Optional[float] = None,
    output: typing.Optional[typing.Union[str, Path]] = None,
) -> SangerReport:
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
    threshold : float, optional
        Threshold on the Phred quality score used to ignore low quality regions
        at the beginning and end of a sample read (default 50). The Phred score
        scales typically between 0 and 62. A good Sanger sequencing read has a
        score of 55.
    output : PathLike, optional
        Output directory for alignment files (default current working
        directory).

    Returns
    -------
    sanger_sequencing.reports.SangerReport
        A collection of plasmid reports.

    Raises
    ------
    AssertionError
        All function arguments are extensively validated and may raise errors.

    See Also
    --------
    plasmid_report

    """
    kwargs = {}
    if threshold is not None:
        kwargs["threshold"] = threshold
    if output is not None:
        kwargs["output"] = output
    report = SangerReport(**kwargs)
    # Initialize global singleton with parameters.
    Configuration(threshold=report.threshold, output=report.output)
    logger.info("Validate template.")
    errors = validation.validate_template(template)
    if errors:
        log_errors(errors)
        raise AssertionError(
            "Invalid analysis template. Please see errors above for details."
        )
    logger.info("Validate plasmids.")
    for plasmid in plasmids.values():
        validation.validate_plasmid(plasmid, [])
    logger.info("Validate samples.")
    for sample in samples.values():
        validation.validate_sample(sample)
    template = validation.drop_missing_records(template, plasmids, samples)
    logger.info("Generate reports.")
    report.plasmids = [
        plasmid_report(plasmid_id, plasmids[plasmid_id], sub, samples)
        for plasmid_id, sub in template.groupby("plasmid", as_index=False, sort=False)
    ]
    return report


def plasmid_report(
    plasmid_id: str,
    sequence: SeqRecord,
    template: DataFrame,
    samples: typing.Dict[str, SeqRecord],
) -> PlasmidReport:
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
    report = PlasmidReport()
    logger.info("Analyze plasmid '%s'.", plasmid_id)
    report = {
        "id": plasmid_id,
        "name": sequence.name,
        "samples": [
            sample_report(
                row.sample, samples[row.sample], row.primer, plasmid_id, sequence,
            )
            for row in template.itertuples(index=False)
        ],
    }
    # Post-process reports in order to classify conflicts.
    logger.debug("Concatenate the detailed sample reports.")
    total = analysis.concatenate_sample_reports(report["samples"])
    for rep in report["samples"]:
        rep["conflicts"] = analysis.summarize_plasmid_conflicts(
            rep["details"], total, sequence
        )
    return report


def sample_report(
    sample_id: str,
    sample_sequence: SeqRecord,
    primer_id: str,
    plasmid_id: str,
    plasmid_sequence: SeqRecord,
) -> SampleReport:
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
    report = SampleReport()
    logger.info("Analyze sample '%s'.", sample_id)
    report = {
        "id": sample_id,
        "primer": primer_id,
        "readLength": len(sample_sequence),
        "errors": [],
    }
    # Convert to base `float` for JSON compatibility.
    try:
        start, trimmed_seq, quality_scores, end, median = analysis.trim_sample(
            sample_sequence
        )
    except ValueError as err:
        report["errors"].append(str(err))
        return report
    report["medianQuality"] = median
    report["cutBeginning"] = int(start)
    report["cutEnd"] = int(end)
    align = analysis.emboss_alignment(
        sample_id, trimmed_seq, plasmid_id, plasmid_sequence
    )
    report["details"] = analysis.alignment_to_table(align, quality_scores, start)
    return report
