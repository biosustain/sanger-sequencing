# Copyright (c) 2019 Novo Nordisk Foundation Center for Biosustainability,
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


"""Summarize results for an entire plasmid."""


__all__ = ("PlasmidReport",)


import typing

from pydantic.dataclasses import dataclass

from .sample_report import SampleReport


class PlasmidReportConfig:
    """
    Configure the `PlasmidReport` behavior.

    Please refer to https://pydantic-docs.helpmanual.io/#model-config for more
    details.

    """

    description = "Summarize results for an entire plasmid."
    validate_all = True
    validate_assignment = True
    fields = {
        "id": {"description": "The given identifier for the plasmid."},
        "name": {"description": "A human readable name for the plasmid."},
        "samples": {"description": "A collection of individual sample (read) reports."},
    }


@dataclass(config=PlasmidReportConfig)
class PlasmidReport:  # noqa: D101

    id: str
    name: str
    samples: typing.List[SampleReport] = ()
