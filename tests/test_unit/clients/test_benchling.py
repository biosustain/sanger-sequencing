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

import pytest

from sanger_sequencing.clients import BenchlingClient


@pytest.fixture(scope="module")
def client():
    return BenchlingClient()


def test_get_plasmid_ids(client, mocker):
    assert False


def test_get_primer_ids(client, mocker):
    assert False


def test_get_plasmid_record(client, mocker):
    assert False
