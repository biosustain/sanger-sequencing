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
from tempfile import mkdtemp

import pydantic
from pydantic import Schema
from pydantic.dataclasses import dataclass

from .plasmid_report import PlasmidReport


class SangerReportConfig:
    """
    Configure the `SangerReport` behavior.

    Please refer to https://pydantic-docs.helpmanual.io/#model-config for more
    details.

    """

    description = "Configure and collect multiple plasmid reports."
    validate_all = True


@dataclass(config=SangerReportConfig)
class SangerReport:  # noqa: D101

    threshold: pydantic.confloat(ge=0.0, le=62.0) = Schema(
        default=50.0,
        description="Threshold on the Phred quality score. The Phred score "
        "scales typically between 0 and 62. A typical good Sanger sequencing "
        "read has a score of 55.",
    )
    output: pydantic.DirectoryPath = Schema(
        default=mkdtemp(), description="Output directory for alignment files."
    )
    plasmids: typing.List[PlasmidReport] = Schema(
        default=[], description="A collection of individual plasmid reports."
    )
