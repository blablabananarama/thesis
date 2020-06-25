# Copyright 2018 Novo Nordisk Foundation Center for Biosustainability, DTU.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import json
import logging
import os

import requests

from simulations.app import app
from simulations.exceptions import PartNotFound
from simulations.metrics import API_REQUESTS
from simulations.utils import Singleton


logger = logging.getLogger(__name__)


class ICE(metaclass=Singleton):
    """ICE API client. Docs: http://ice.jbei.org/api/index.html"""

    def __init__(self):
        """On instantiation, request and store a session id for later use."""
        if os.environ["ENVIRONMENT"] in ("production", "staging"):
            self._update_session_id()
        else:
            # To speed up development, don't set a valid session id on init but rely on
            # the re-authentication logic should ICE be needed.
            self.SESSION_ID = ""

    def get_reaction_equations(self, genotype):
        """
        Request genotype part info from ICE and return reaction map information from the
        references field.
        """
        logger.info(f"Requesting genotype '{genotype}' from ICE")
        with API_REQUESTS.labels(
            "model", os.environ["ENVIRONMENT"], "ice", app.config["ICE_API"]
        ).time():
            response = requests.get(
                f"{app.config['ICE_API']}/rest/parts/{genotype}",
                headers=self._headers(),
            )

        # In case of authentication failure, get a new session id and re-try the
        # request. The ICE documentation says nothing about access token expiry, so it's
        # not clear whether this can or will ever occur.
        if response.status_code in (401, 403):
            self._update_session_id()
            with API_REQUESTS.labels(
                "model", os.environ["ENVIRONMENT"], "ice", app.config["ICE_API"]
            ).time():
                response = requests.get(
                    f"{app.config['ICE_API']}/rest/parts/{genotype}",
                    headers=self._headers(),
                )

        # If the part is not found, raise the corresponding exception
        if response.status_code == 404:
            raise PartNotFound()

        response.raise_for_status()
        result = response.json()

        # Extract reaction id -> string map from the expected format in the references
        # field. Note: This could be removed in the future, in favor of storing only the
        # reaction IDs in ICE, and looking up the corresponding stoichiometry/reaction
        # equations from the BiGG database instead.
        reactions = result["references"].split(",")
        reaction_tuples = [reaction.split(":") for reaction in reactions]
        reactions_map = {id.strip(): string.strip() for id, string in reaction_tuples}
        return reactions_map

    def _update_session_id(self):
        """
        Query ICE for a new access token. Note that this usually takes ~10 seconds!
        """
        logger.info(f"Requesting session token from ICE")
        response = requests.post(
            f"{app.config['ICE_API']}/rest/accesstokens",
            headers=self._headers(add_session_id=False),
            data=json.dumps(
                {
                    "email": app.config["ICE_USERNAME"],
                    "password": app.config["ICE_PASSWORD"],
                }
            ),
        )
        response.raise_for_status()
        self.SESSION_ID = response.json()["sessionId"]

    def _headers(self, add_session_id=True):
        """Return headers dict for use with ICE API requests."""
        headers = {"Content-Type": "application/json"}
        if add_session_id:
            headers.update({"X-ICE-Authentication-SessionId": self.SESSION_ID})
        return headers
