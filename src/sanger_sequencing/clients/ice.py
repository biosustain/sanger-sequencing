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

"""A basic client for the jbei ICE API."""

from .repository_client import RepositoryClient


__all__ = ("ICEClient",)


class ICEClient(RepositoryClient):
    """Provide a pythonic interface to the ICE API."""

    def __init__(self):
        """Initialize the ICE client."""
        super().__init__()

    def get_plasmid_ids(self):
        """Return a frozenset of all accessible plasmid identifiers."""
        raise NotImplementedError("Coming soon.")

    def get_plasmid_record(self, plasmid_id):
        """Return the name and sequence record for a specific plasmid."""
        raise NotImplementedError("Coming soon.")

    def get_primer_ids(self):
        """Return a frozenset of all accessible primer identifiers."""
        raise NotImplementedError("Coming soon.")
