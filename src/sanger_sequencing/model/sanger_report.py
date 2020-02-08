# Copyright (c) 2018-2020 Novo Nordisk Foundation Center for Biosustainability,
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


"""Summarize results for multiple plasmids and samples."""


import typing
from pathlib import Path

import pydantic
from pydantic import BaseModel, Field

from .plasmid_report import PlasmidReport, PlasmidReportInternal


__all__ = ("SangerReport", "SangerReportInternal")


class BaseSangerReport(BaseModel):
    """Define attributes for a base Sanger read report."""

    threshold: pydantic.confloat(ge=0.0, le=62.0) = Field(
        50.0,
        description="Threshold on the Phred quality score used to ignore low quality "
        "regions at the beginning and end of a sample read. The Phred score"
        " scales typically between 0 and 62. A good Sanger sequencing "
        "read has a score of 55.",
    )

    class Config:
        """Configure the base Sanger report behavior."""

        description = "Configure and collect multiple plasmid reports."


class SangerReport(BaseSangerReport):
    """Define attributes for a Sanger read report used in external communication."""

    plasmids: typing.List[PlasmidReport] = Field(
        (), description="A collection of individual plasmid reports."
    )

    class Config(BaseSangerReport.Config):
        """Configure the Sanger report behavior."""

        allow_population_by_field_name = True
        orm_mode = True


class SangerReportInternal(BaseSangerReport):
    """Define attributes for an internal Sanger read report."""

    output: pydantic.DirectoryPath = Field(
        Path.cwd(), description="Output directory for alignment files."
    )
    plasmids: typing.List[PlasmidReportInternal] = Field(
        (), description="A collection of individual plasmid reports for internal use."
    )

    class Config(BaseSangerReport.Config):
        """Configure the internal Sanger report behavior."""

        validate_all = True
        validate_assignment = True
