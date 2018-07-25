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

"""A basic client for the LabCollector API v2."""

import base64
from functools import lru_cache
from io import StringIO
from operator import itemgetter
from typing import FrozenSet, Tuple
from urllib.parse import urljoin

import requests
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .repository_client import RepositoryClient


__all__ = ("LabCollectorClient",)


class LabCollectorClient(RepositoryClient):
    """Provide a pythonic interface to the LabCollector API v2."""

    COUNT_PARAMETERS = {
        "fields": "count"
    }
    FIELD_PARAMETERS = {
        "fields": "Sequence_file_GenBank"
    }

    def __init__(self,
                 api: str,
                 token: str,
                 timeout: float=10.0,
                 cache_size: int=10_000,
                 **kwargs):
        """
        Initialize a Labcollector API v2 client.

        Parameters
        ----------
        api : str
            The base Labcollector API URL; typically similar to
            https://labcollector.your.domain/webservice/v2/.
        token : str
            A token that has read access to plasmids and primers.
        timeout : float
            A timeout in seconds for the GET requests.
        cache_size : int, optional
            The maximal size of the LRU cache. A trade-off between
            performance for already requested resources and the memory used.
        kwargs : optional
            Additional keyword arguments are used to update the default
            header for making requests to the API.

        """
        super().__init__()
        self.api = api
        if not self.api.endswith("/"):
            self.api += "/"
        self.token = token
        self.timeout = timeout
        self.headers = {
            'X-LC-APP-Auth': token,
            'Accept': "application/json",
            'Cache-Control': "no-cache"
        }
        self.headers.update(kwargs)
        self.plasmids_resource = urljoin(self.api, "plasmids")
        self.primers_resource = urljoin(self.api, "primers")
        # We want the instance methods to be cached rather than the class
        # methods.
        self.get_plasmid_ids = lru_cache(maxsize=1)(self.get_plasmid_ids)
        self.get_plasmid_record = lru_cache(maxsize=int(cache_size))(
            self.get_plasmid_record)
        self.get_primer_ids = lru_cache(maxsize=1)(self.get_primer_ids)

    def _get_resource_ids(self, endpoint: str) -> FrozenSet:
        response = requests.get(
            endpoint, headers=self.headers, params=self.COUNT_PARAMETERS,
            timeout=self.timeout)
        response.raise_for_status()
        # The counter seems to be the sequential identifier for objects
        # within the Labcollector database.
        count_getter = itemgetter("count")
        return frozenset(count_getter(elem) for elem in response.json())

    def get_plasmid_ids(self) -> FrozenSet:
        """Return a frozenset of all accessible plasmid identifiers."""
        return self._get_resource_ids(self.plasmids_resource)

    def get_primer_ids(self) -> FrozenSet:
        """Return a frozenset of all accessible primer identifiers."""
        return self._get_resource_ids(self.primers_resource)

    def _get_sequence_record(self,
                             endpoint: str,
                             record_id: str) -> Tuple[str, SeqRecord]:
        endpoint += f"/{record_id}"
        response = requests.get(
            endpoint, headers=self.headers, params=self.FIELD_PARAMETERS,
            timeout=self.timeout)
        response.raise_for_status()
        data = response.json()[0]
        name = data["Sequence_file_GenBank"][0]["name"]
        content = data["Sequence_file_GenBank"][0]["content"]
        return name, SeqIO.read(
            StringIO(base64.b64decode(content).decode("ascii")), "gb")

    def get_plasmid_record(self, plasmid_id: str) -> Tuple[str, SeqRecord]:
        """
        Return the name and sequence record for a specific plasmid.

        Parameters
        ----------
        plasmid_id : str
            The plasmid identifier of interest.

        Returns
        -------
        (str, Bio.SeqRecord.SeqRecord)
            The plasmid name and sequence record if any.

        """
        return self._get_sequence_record(self.plasmids_resource, plasmid_id)
