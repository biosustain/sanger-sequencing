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


"""Summarize results for multiple plasmids and samples."""


__all__ = ("SangerReport",)


import typing
from pathlib import Path

import pydantic
from pydantic import BaseModel, Field

from .plasmid_report import PlasmidReport


class SangerReport(BaseModel):  # noqa: D101

    threshold: pydantic.confloat(ge=0.0, le=62.0) = Field(
        50.0,
        description="Threshold on the Phred quality score used to ignore low quality "
        "regions at the beginning and end of a sample read. The Phred score"
        " scales typically between 0 and 62. A good Sanger sequencing "
        "read has a score of 55.",
    )
    output: pydantic.DirectoryPath = Field(
        Path.cwd(), description="Output directory for alignment files."
    )
    plasmids: typing.List[PlasmidReport] = Field(
        (), description="A collection of individual plasmid reports."
    )

    class Config:
        """
        Configure the `SangerReport` behavior.

        Please refer to https://pydantic-docs.helpmanual.io/#model-config for more
        details.

        """

        description = "Configure and collect multiple plasmid reports."
        validate_all = True
        validate_assignment = True
