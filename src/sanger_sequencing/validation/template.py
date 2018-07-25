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

"""Validate an analysis template."""

import json
import logging
from typing import Dict, List

from Bio.SeqRecord import SeqRecord
from goodtables import validate
from importlib_resources import open_text
from pandas import DataFrame

from . import schemata


__all__ = ("validate_template", "drop_missing_records")

LOGGER = logging.getLogger(__name__)
with open_text(schemata, "template.json", encoding="utf-8") as file_handle:
    TEMPLATE_SCHEMA = json.load(file_handle)


def validate_template(template: DataFrame) -> List[Dict]:
    """
    Validate that the template respects the desired schema.

    Make use of goodtables [1]_ in order to validate the template data frame
    against the chosen specification.

    Parameters
    ----------
    template : pandas.DataFrame
        A data frame with three columns: 'plasmid', 'primer', 'sample'.

    Returns
    -------
    list
        List of errors that are themselves dictionaries with the keys 'code'
        and 'message'.

    References
    ----------
    .. [1] https://pypi.org/project/goodtables/

    """
    if len(template) == 0:
        return [
            {
                "code": "empty",
                "message": "There must be at least one data row."
            }
        ]
    # Convert the data frame into a format suitable for goodtables, the keys
    # of the first record are taken as the column headers.
    records = template.to_dict("records")
    result = validate(
        records, preset="table", headers=list(records[0]),
        schema=TEMPLATE_SCHEMA, order_fields=True)["tables"][0]
    return result.get("errors", [])


def drop_missing_records(
    template: DataFrame,
    plasmids: Dict[str, SeqRecord],
    samples: Dict[str, SeqRecord]
) -> DataFrame:
    """
    Drop rows with missing sequence records from the template.

    Template rows that lack a corresponding plasmid sequence or sample
    sequence are removed from the template thus ensuring a consistent analysis.

    Parameters
    ----------
    template : pandas.DataFrame
        A data frame with three columns: 'plasmid', 'primer', 'sample'.
    plasmids : dict
        A dictionary mapping plasmid names to their sequence records.
    samples : dict
        A dictionary mapping sample names to their sequence records.

    Returns
    -------
    pandas.DataFrame
        Return the template only for those rows where sequences exist.

    """
    plasmid_mask = template["plasmid"].isin(plasmids)
    if (~plasmid_mask).any():
        LOGGER.error(
            "The following plasmid(s) have no corresponding sequence record: "
            "%s.", ", ".join(template.loc[plasmid_mask, "plasmid"]))
    sample_mask = template["sample"].isin(samples)
    if (~sample_mask).any():
        LOGGER.error(
            "The following sample(s) have no corresponding sequence record: "
            "%s.", ", ".join(template.loc[plasmid_mask, "sample"]))
    return template.loc[plasmid_mask & sample_mask, :]
