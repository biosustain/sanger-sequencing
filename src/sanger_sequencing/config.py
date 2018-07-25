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

"""Provide a global configuration singleton."""

import logging
from tempfile import mkdtemp


__all__ = ("Configuration",)

LOGGER = logging.getLogger(__name__)


class Singleton(type):
    """Implementation of the singleton pattern."""

    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(
                *args, **kwargs)
        return cls._instances[cls]


class Configuration(metaclass=Singleton):
    """
    Configure values that affect how the analysis is performed.

    Attributes
    ----------
    threshold : float
        Threshold on the Phred quality. The Phred score scales between 0 and
        62. A typical good Sanger sequencing read has a score of 55.
    output : str or pathlib.Path
        Output directory for alignment files.

    """

    def __init__(self, threshold=50.0, output=mkdtemp(), **kwargs):
        """
        Initialize the singleton configuration object.

        Parameters
        ----------
        threshold : float, optional
            Threshold on the Phred quality. The Phred score scales between 0 and
            62. A typical good Sanger sequencing read has a score of 55.
        output : str or pathlib.Path, optional
            Output directory for alignment files.
        kwargs : dict
            Passed to the ``super()`` class.

        """
        super().__init__(**kwargs)
        self.threshold = threshold
        self.output = output
