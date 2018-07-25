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

"""Provide sample analysis functions."""

import logging

from Bio.SeqRecord import SeqRecord
from numpy import arange, array, asarray, nanmedian

from ..config import Configuration


__all__ = ("trim_sample",)

LOGGER = logging.getLogger(__name__)


def trim_sample(seq: SeqRecord) -> (int, SeqRecord, array, int, float):
    """Cut off low quality ends of a Sanger sequencing record."""
    LOGGER.debug("Trim sample.")
    config = Configuration()
    scores = asarray(seq.letter_annotations['phred_quality'])
    median = float(nanmedian(scores))
    if median < config.threshold:
        message = (
            f"The median Phred quality ({median}) is below the "
            f"required threshold ({config.threshold}).")
        LOGGER.error(message)
        raise ValueError(message)
    mask = (scores >= config.threshold)
    index = arange(len(mask), dtype=int)
    min_i = index[mask][0]
    max_i = index[mask][-1] + 1  # Since Python excludes upper range limit.
    LOGGER.debug("Cutting %d nucleotides at the beginning and %d at the end.",
                 min_i + 1, len(seq) - max_i)
    return (
        min_i, seq[min_i:max_i], scores[min_i:max_i], len(seq) - max_i, median)
