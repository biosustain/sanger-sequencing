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

"""Verify the validation functions."""

from io import StringIO

import pytest
from pandas import DataFrame, read_csv

import sanger_sequencing.validation as validation


@pytest.fixture(scope="module")
def template():
    return read_csv(StringIO(
        """
        plasmid,primer,sample
        gfp,47,abcd
        gfp,47,efgh
        gfp,42,1
        gfp,42,2
        rfp,k,xyz
        rfp,l,uvw
        rfp,k,3
        rfp,l,4
        """
    ), skipinitialspace=True)


def test_validate_template(template):
    errors = validation.validate_template(template)
    assert len(errors) == 0


def test_validate_empty_template():
    errors = validation.validate_template(DataFrame())
    assert len(errors) == 1
    assert errors[0]["code"] == "empty"
    assert errors[0]["message"] == "There must be at least one data row."


@pytest.mark.parametrize("plasmids, samples, indeces", [
    ({"gfp": None},
     {"abcd": None, "efgh": None, "1": None, "2": None},
     [0, 1, 2, 3]),
    ({"rfp": None},
     {"xyz": None, "uvw": None, "3": None, "4": None},
     [4, 5, 6, 7]),
    ({"gfp": None},
     {"1": None, "2": None},
     [2, 3]),
    ({"rfp": None},
     {"xyz": None, "uvw": None},
     [4, 5])
])
def test_drop_missing_records(template, plasmids, samples, indeces):
    reduced = validation.drop_missing_records(template, plasmids, samples)
    assert (reduced.index == indeces).all()
