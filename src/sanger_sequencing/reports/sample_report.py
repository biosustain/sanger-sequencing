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


__all__ = ("SampleReport",)


import typing

import pydantic
from pydantic import BaseModel, Field


class SampleReport(BaseModel):  # noqa: D101

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
    median_quality: pydantic.confloat(ge=0.0, le=62.0) = Field(
        ...,
        alias="medianQuality",
        title="Median Quality",
        description="The median Phred quality of the sample read.",
    )
    trim_start: pydantic.conint(ge=0) = Field(
        ...,
        alias="trimStart",
        title="Trim Start",
        description="The number of nucleotides trimmed at the beginning of the "
                    "sequence due to low Phred quality and before alignment.",
    )
    trim_end: pydantic.conint(ge=0) = Field(
        ...,
        alias="trimEnd",
        title="Trim End",
        description="The number of nucleotides trimmed at the end of the sequence due "
                    "to low Phred quality and before alignment.",
    )
    details: typing.Any = Field(..., description="A table of conflicts.")
    errors: typing.List[str] = Field(
        (), description="Errors in aligning the sample read."
    )

    class Config:
        """
        Configure the `SampleReport` behavior.

        Please refer to https://pydantic-docs.helpmanual.io/#model-config for more
        details.

        """

        description = "Summarize results for a sample read."
        validate_all = True
        validate_assignment = True
