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

"""Provide helper functions."""

import logging
from typing import Dict, Iterable

from depinfo import print_dependencies


__all__ = ("log_errors", "show_versions")

LOGGER = logging.getLogger(__name__)


def log_errors(errors: Iterable[Dict]):
    """Log all error objects from an iterable."""
    for err in errors:
        LOGGER.error("(%s) %s", err["code"], err["message"])


def show_versions():
    """Print all dependencies with their versions."""
    print_dependencies("sanger-sequencing")
