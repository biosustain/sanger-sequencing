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

"""Provide Sanger sequence alignment functions and scoring."""

import logging
import re
from os.path import join
from typing import List, get_type_hints

from Bio import AlignIO
from Bio.Emboss.Applications import WaterCommandline
from Bio.SeqRecord import SeqRecord
from numpy import array, nan
from pandas import DataFrame

from ..config import Configuration


__all__ = ("emboss_alignment", "alignment_to_table")

LOGGER = logging.getLogger(__name__)


def alignment_to_table(align: AlignIO.MultipleSeqAlignment,
                       scores: array,
                       start: int) -> DataFrame:
    """Generate a table of the alignment."""
    report = list()
    b_coord = align.positions['bseq_start']
    a_coord = align.positions['aseq_start']
    for a_char, b_char in zip(align[0], align[1]):
        if a_char == b_char:
            report.append((a_coord, b_coord, False, a_char, b_char))
            a_coord += 1
            b_coord += 1
        elif b_char == '-':
            report.append((a_coord, nan, True, a_char, b_char))
            a_coord += 1
        elif a_char == '-':
            report.append((nan, b_coord, True, a_char, b_char))
            b_coord += 1
        else:
            report.append((a_coord, b_coord, True, a_char, b_char))
            a_coord += 1
            b_coord += 1
    df = DataFrame(report, columns=[
        'plasmid_pos', 'sample_pos', 'snp', 'plasmid_chr', 'sample_chr'])
    # Subtract 1 for Python indexing into scores.
    df["quality"] = (
        df.loc[df["sample_pos"].notnull(), "sample_pos"].astype(int) - 1).map(
        scores.__getitem__)
    # Re-align to original sample sequence position.
    df.loc[df["sample_pos"].notnull(), "sample_pos"] += start
    return df


def extract_emboss_positions(lines: List[str]):
    """Extract the alignment positions from the `water` output file."""
    asis_lines = [l.strip() for l in lines if l.startswith('asis')]
    start_pattern = re.compile(r'^asis\s+(\d+)')
    end_pattern = re.compile(r'(\d+)$')
    return {
        'aseq_start': int(start_pattern.match(asis_lines[0])[1]),
        'bseq_start': int(start_pattern.match(asis_lines[1])[1]),
        'aseq_end': int(end_pattern.search(asis_lines[-2])[1]),
        'bseq_end': int(end_pattern.search(asis_lines[-1])[1])
    }


def emboss_alignment(sample_id: str,
                     sample_sequence: SeqRecord,
                     plasmid_id: str,
                     plasmid_sequence: SeqRecord,
                     gap_open_penalty: float=2.0,
                     gap_extension_penalty: float=10.0,
                     tool: get_type_hints(WaterCommandline)=WaterCommandline
                     ) -> AlignIO.MultipleSeqAlignment:
    """
    Create an alignment between the known plasmid sequence and the Sanger read.

    We generally expect the sequences to be identical excluding single
    nucleotide errors coming from real mutations, replication errors, or bad
    reads. Read quality is not yet taken into account but will in future.

    Parameters
    ----------
    plasmid_sequence
    sample_sequence
    sample_id
    plasmid_id
    gap_open_penalty
    gap_extension_penalty
    tool

    Returns
    -------
    Bio.AlignIO.MultipleSeqAlignment
        The pairwise alignment.

    """
    config = Configuration()
    # Possibly the output of `aformat="markx10"` is easier to parse.
    cmd = tool(gapopen=gap_open_penalty, gapextend=gap_extension_penalty)
    cmd.asequence = f"asis:{plasmid_sequence.seq}"
    cmd.bsequence = f"asis:{sample_sequence.seq}"
    outfile = join(config.output, f"{sample_id}_{plasmid_id}.txt")
    cmd.outfile = outfile
    stdout, stderr = cmd()
    LOGGER.debug(stdout)
    LOGGER.debug(stderr)  # Unfortunately, normal messages are sent to stderr...

    align = AlignIO.read(cmd.outfile, "emboss")
    # Get coordinates manually because biopython ...
    identity = align.annotations['identity'] / len(sample_sequence)
    LOGGER.debug("Sequence identity is %0.2g.", identity)
    if identity < 0.9:
        LOGGER.info("Trying reverse complement!")
        cmd.bsequence = f"asis:{sample_sequence.reverse_complement().seq}"
        rev_outfile = join(
            config.output, "{sample_id}_{plasmid_id}_rev.txt")
        cmd.outfile = rev_outfile
        stdout, stderr = cmd()
        LOGGER.debug(stdout)
        LOGGER.debug(stderr)
        rev_align = AlignIO.read(cmd.outfile, "emboss")
        rev_identity = rev_align.annotations['identity'] / len(sample_sequence)
        LOGGER.debug("Complement sequence identity is %0.2g.", rev_identity)
    else:
        rev_identity = -1.0
        rev_align = None
        rev_outfile = None
    if identity > rev_identity:
        with open(outfile) as file_h:
            align.positions = extract_emboss_positions(file_h.readlines())
        LOGGER.debug(str(align.positions))
        return align
    else:
        with open(rev_outfile) as file_h:
            rev_align.positions = extract_emboss_positions(file_h.readlines())
        LOGGER.debug(str(rev_align.positions))
        return rev_align
