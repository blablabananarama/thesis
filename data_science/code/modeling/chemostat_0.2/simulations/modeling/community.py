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

import logging
import tempfile
import warnings

import cobra
import reframed

from simulations.modeling.reframed_helpers import (
    create_metabolite_id2name_mapping,
    generate_transactions,
)


logger = logging.getLogger(__name__)

METHODS = ["steadycom", "steadiercom"]


def simulate(wrappers, medium, method):
    """
    Run a SteadyCom community simulation.

    Parameters
    ----------
    wrappers: list(storage.ModelWrapper)
        A list of model wrappers containing cobrapy model instances.
    medium: list(str)
        A list of compound names. Exchange reaction identifiers are assumed to
        be formatted according to: "EX_{compound}_e"
    method: str
        The community simulation method. Currently accepted strings:
        "steadycom" or "steadiercom".
    """
    if method not in METHODS:
        raise ValueError(f"Unsupported community simulation method '{method}'")

    with warnings.catch_warnings(record=True) as reframed_warnings:
        logger.debug("Converting cobrapy models to reframed models")
        rf_models = []
        for model in [wrapper.model for wrapper in wrappers]:
            # The most funcational approach (albeit slow) seems to be to write and
            # reload SBML. reframed's cobrapy integration is currently pretty
            # minimal.
            with tempfile.NamedTemporaryFile() as file_:
                cobra.io.write_sbml_model(model, file_.name)
                # TODO: Consider accepting the flavor argument as a parameter
                # instead of always assuming BiGG.
                rf_models.append(reframed.load_cbmodel(file_.name, flavor="bigg"))

        logger.debug("Merging individual models to a community")
        community = reframed.Community("community", rf_models)

        logger.debug("Applying medium to the community")
        environment = reframed.Environment.from_compounds(
            medium, fmt_func=lambda x: f"R_EX_M_{x}_e"
        )
        environment.apply(community.merged_model, inplace=True)

        if method == "steadycom":
            logger.info(f"Simulating community model with SteadyCom")
            solution = reframed.SteadyCom(community)
        elif method == "steadiercom":
            logger.info(f"Simulating community model with SteadyCom")
            solution = reframed.SteadierCom(community)

        logger.debug(f"Formatting solution response")

        def model_id(original_id):
            """Map the models original name back to our platform internal DB IDs."""
            return next(
                wrapper.id for wrapper in wrappers if wrapper.model.id == original_id
            )

        # Calculate transactions (cross-feeding, uptake and secretion)
        logger.debug(f"Calculating transactions (cross-feeding, uptake and secretion)")
        ex_met_ids = solution.community.merged_model.get_external_metabolites()
        met_id2name = create_metabolite_id2name_mapping(ex_met_ids, community)
        transactions = generate_transactions(met_id2name, solution.exchange_map)

        # Convert the iterables to dictionaries for easier handling on the frontend
        abundance = [
            {"id": model_id(original_id), "value": abundance}
            for original_id, abundance in solution.abundance.items()
        ]
        cross_feeding = []
        for transaction in transactions:
            if transaction[0] == "medium":
                cross_feeding.append(
                    {
                        "from": "medium",
                        "to": model_id(transaction[1]),
                        "metabolite_id": transaction[2],
                        "metabolite_name": transaction[3],
                        "value": transaction[4],
                    }
                )
            elif transaction[1] == "medium":
                cross_feeding.append(
                    {
                        "from": model_id(transaction[0]),
                        "to": "medium",
                        "metabolite_id": transaction[2],
                        "metabolite_name": transaction[3],
                        "value": transaction[4],
                    }
                )
            else:
                cross_feeding.append(
                    {
                        "from": model_id(transaction[0]),
                        "to": model_id(transaction[1]),
                        "metabolite_id": transaction[2],
                        "metabolite_name": transaction[3],
                        "value": transaction[4],
                    }
                )
    return {
        "growth_rate": solution.growth,
        "abundance": abundance,
        "cross_feeding": cross_feeding,
        "warnings": [" ".join(w.message.args) for w in reframed_warnings],
    }
