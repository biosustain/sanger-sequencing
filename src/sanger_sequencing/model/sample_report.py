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


"""Summarize results for a single sample (read)."""


__all__ = ("SampleReport", "SampleReportInternal")


import typing

import pydantic
from pandas import DataFrame
from pydantic import BaseModel, Field

from .conflict_report import ConflictReport, ConflictReportInternal


class BaseSampleReport(BaseModel):
    """Define attributes for a base sample report."""

    id: str = Field(..., description="The given identifier for the sample.")
    primer: str = Field(
        ...,
        description="The primer identifier or name that was used for this sample read.",
    )
    read_length: pydantic.PositiveInt = Field(
        ...,
        alias="readLength",
        title="Read Length",
        description="The number of nucleotides of this sample read.",
    )
    median_quality: typing.Optional[pydantic.confloat(ge=0.0, le=62.0)] = Field(
        None,
        alias="medianQuality",
        title="Median Quality",
        description="The median Phred quality of the sample read.",
    )
    trim_start: typing.Optional[pydantic.conint(ge=0)] = Field(
        None,
        alias="trimStart",
        title="Trim Start",
        description="The number of nucleotides trimmed at the beginning of the "
        "sequence due to low Phred quality and before alignment.",
    )
    trim_end: typing.Optional[pydantic.conint(ge=0)] = Field(
        None,
        alias="trimEnd",
        title="Trim End",
        description="The number of nucleotides trimmed at the end of the sequence due "
        "to low Phred quality and before alignment.",
    )
    errors: typing.List[str] = Field(
        (), description="Errors in aligning the sample read."
    )

    class Config:
        """Configure the base sample report behavior."""

        description = "Summarize results for a sample read."


class SampleReport(BaseSampleReport):
    """Define attributes for a sample report used in external communication."""

    conflicts: typing.Optional[typing.List[ConflictReport]] = Field(
        (), description="A summary of the conflicts detected in this sample read."
    )

    class Config(BaseSampleReport.Config):
        """Configure the sample report behavior."""

        allow_population_by_field_name = True
        orm_mode = True


class SampleReportInternal(BaseSampleReport):
    """Define attributes for an internal plasmid report."""

    conflicts: typing.Optional[typing.List[ConflictReportInternal]] = Field(
        (), description="A summary of the conflicts detected in this sample read."
    )
    details: typing.Optional[DataFrame] = Field(
        None, description="A table of conflicts."
    )

    class Config(BaseSampleReport.Config):
        """Configure the internal sample report behavior."""

        # Allow the `details` field to be a `pandas.DataFrame`.
        arbitrary_types_allowed = True
        validate_all = True
        validate_assignment = True
