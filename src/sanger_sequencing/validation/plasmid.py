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

"""Validate a plasmid sequence record."""

import logging
from typing import Iterable

from Bio.SeqRecord import SeqRecord


__all__ = ("validate_plasmid",)

LOGGER = logging.getLogger(__name__)


def validate_plasmid(plasmid: SeqRecord, primer_ids: Iterable[str]):
    """
    Validate that the plasmid sequence record has the required annotations.

    Parameters
    ----------
    plasmid : Bio.SeqRecord.SeqRecord
        A plasmid sequence record that must contain the following annotations.
    primer_ids : iterable
        The given primer identifiers coming from the analysis template should
        be annotated in the plasmid sequence record.

    Returns
    -------
    list
        List of errors that are themselves dictionaries with the keys 'code'
        and 'message'.

    """
    assert len(plasmid) > 0
    assert hasattr(plasmid, "features")
    for feat in plasmid.features:
        assert hasattr(feat, "location")
        assert hasattr(feat.location, "start")
        assert hasattr(feat.location.start, "position")
        assert 0 <= feat.location.start.position <= len(plasmid)
        assert hasattr(feat.location, "end")
        assert hasattr(feat.location.end, "position")
        assert 0 <= feat.location.end.position <= len(plasmid)
        assert hasattr(feat, "type")
        assert hasattr(feat, "qualifiers")
        assert "label" in feat.qualifiers
    for primer_id in primer_ids:
        # TODO: Test that all given primer identifiers are on the plasmid.
        pass
