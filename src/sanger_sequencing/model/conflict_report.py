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


"""Summarize results for a single conflict."""


from enum import Enum
from typing import List, Optional, Tuple

import pydantic
from pydantic import BaseModel, Field


__all__ = (
    "ConflictReport",
    "ConflictReportInternal",
    "ConflictTypeEnum",
    "QualityEnum",
)


class ConflictTypeEnum(str, Enum):

    INSERTION = "insertion"
    DELETION = "deletion"
    CHANGE = "change"


class QualityEnum(str, Enum):

    HIGH = "high"
    LOW = "low"


class BaseConflictReport(BaseModel):
    """Define attributes for a base conflict report."""

    plasmid_character: str = Field(
        ..., alias="plasmidCharacter", title="Plasmid Character"
    )
    sample_character: str = Field(
        ..., alias="sampleCharacter", title="Sample Character"
    )
    plasmid_position: Optional[pydantic.PositiveInt] = Field(
        None, alias="plasmidPosition", title="Plasmid Position"
    )
    sample_position: Optional[pydantic.PositiveInt] = Field(
        None, alias="samplePosition", title="Sample Position"
    )
    type: Optional[ConflictTypeEnum] = None
    surrounding_quality: Optional[QualityEnum] = Field(
        None, alias="surroundingQuality", title="Surrounding Quality"
    )
    num_confirmed: Optional[pydantic.conint(ge=0)] = Field(
        None, alias="numConfirmed", title="Number Confirmed"
    )
    num_invalidated: Optional[pydantic.conint(ge=0)] = Field(
        None, alias="numInvalidated", title="Number Invalidated"
    )
    features_hit: List[Tuple[str, List[str]]] = Field(
        (), alias="featuresHit", title="Hit Features"
    )
    effect: List[str] = ()

    class Config:
        """Configure the base conflict report behavior."""

        description = "Summarize a detected sequencing conflict."
        allow_population_by_field_name = True


class ConflictReport(BaseConflictReport):
    """Define attributes for a conflict report used in external communication."""

    class Config(BaseConflictReport.Config):
        """Configure the external conflict report behavior."""

        orm_mode = True


class ConflictReportInternal(BaseConflictReport):
    """Define attributes for an internal conflict report."""

    class Config(BaseConflictReport.Config):
        """Configure the internal conflict report behavior."""

        validate_all = True
        validate_assignment = True
