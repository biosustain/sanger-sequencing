# Copyright (c) 2019-2020 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""Summarize results for an entire plasmid."""


import typing

from pydantic import BaseModel, Field

from .sample_report import SampleReport, SampleReportInternal


__all__ = ("PlasmidReport", "PlasmidReportInternal")


class BasePlasmidReport(BaseModel):
    """Define attributes for a base plasmid report."""

    id: str = Field(..., description="The given identifier for the plasmid.")
    name: str = Field(..., description="A human readable name for the plasmid.")

    class Config:
        """Configure the base plasmid report behavior."""

        description = "Summarize results for an entire plasmid."


class PlasmidReport(BasePlasmidReport):
    """Define attributes for a plasmid report used in external communication."""

    samples: typing.List[SampleReport] = Field(
        (), description="A collection of individual sample (read) reports."
    )

    class Config(BasePlasmidReport.Config):
        """Configure the plasmid report behavior."""

        allow_population_by_field_name = True
        orm_mode = True


class PlasmidReportInternal(BasePlasmidReport):
    """Define attributes for an internal plasmid report."""

    samples: typing.List[SampleReportInternal] = Field(
        (), description="A collection of individual internal sample (read) reports."
    )

    class Config(BasePlasmidReport.Config):
        """Configure the internal plasmid report behavior."""

        validate_all = True
        validate_assignment = True
