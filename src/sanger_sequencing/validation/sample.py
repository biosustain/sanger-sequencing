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

"""Validate a sample sequence record."""

import logging

from Bio.SeqRecord import SeqRecord


__all__ = ("validate_sample",)

LOGGER = logging.getLogger(__name__)


def validate_sample(sample: SeqRecord):
    """
    Validate that the sample sequence record has Phred quality scores.

    Parameters
    ----------
    sample : Bio.SeqRecord.SeqRecord
        A sample sequence record that must contain the following annotations:

    Returns
    -------
    list
        List of errors that are themselves dictionaries with the keys 'code'
        and 'message'.

    """
    assert hasattr(sample, "letter_annotations")
    assert "phred_quality" in sample.letter_annotations
    assert len(sample.letter_annotations["phred_quality"]) == len(sample)
