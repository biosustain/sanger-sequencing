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
from Bio.SeqRecord import SeqRecord

from sanger_sequencing.clients import LabCollectorClient


TOKEN = "ABCDEFG"


@pytest.fixture(scope="module")
def client():
    return LabCollectorClient(
        api="https://teapot.com/webservice/v2/",
        token=TOKEN,
        timeout=5,
        cache_size=1,
    )


@pytest.mark.parametrize(
    "key, value",
    [
        ("X-LC-APP-Auth", TOKEN),
        ("Accept", "application/json"),
        ("Cache-Control", "no-cache"),
    ],
)
def test_headers(client, key, value):
    """The LabCollector API requires these headers to be set."""
    assert key in client.headers
    assert client.headers[key] == value


def test_get_plasmid_ids(client, mocker):
    pids = frozenset(range(10))
    get = mocker.patch(
        "sanger_sequencing.clients.LabCollectorClient._get_resource_ids",
        return_value=pids,
    )
    assert client.get_plasmid_ids() == pids
    get.assert_called_once_with(client.plasmids_resource)
    # Test caching. There must not be another call to _get_resource_ids.
    assert client.get_plasmid_ids() == pids
    assert get.call_count == 1


def test_get_primer_ids(client, mocker):
    pids = frozenset(range(5))
    get = mocker.patch(
        "sanger_sequencing.clients.LabCollectorClient._get_resource_ids",
        return_value=pids,
    )
    assert client.get_primer_ids() == pids
    get.assert_called_once_with(client.primers_resource)
    # Test caching. There must not be another call to _get_resource_ids.
    assert client.get_primer_ids() == pids
    assert get.call_count == 1


def test_get_plasmid_record(client, mocker):
    plasmid_id = "1234"
    expected = ("slick", SeqRecord("ATGC"))
    get = mocker.patch(
        "sanger_sequencing.clients.LabCollectorClient._get_sequence_record",
        return_value=expected,
    )
    result = client.get_plasmid_record(plasmid_id)
    assert result[0] == expected[0]
    assert result[1].seq == expected[1].seq
    get.assert_called_once_with(client.plasmids_resource, plasmid_id)
    # Test caching. There must not be another call to _get_sequence_record.
    assert client.get_plasmid_record(plasmid_id) == expected
    assert get.call_count == 1
