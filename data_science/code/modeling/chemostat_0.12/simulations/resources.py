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

from cobra.exceptions import OptimizationError
from flask import Response, abort, jsonify, request
from flask_apispec import use_kwargs
from flask_apispec.extension import FlaskApiSpec
from prometheus_client import CONTENT_TYPE_LATEST, CollectorRegistry, generate_latest
from prometheus_client.multiprocess import MultiProcessCollector

from simulations import storage
from simulations.exceptions import Forbidden, ModelNotFound, Unauthorized
from simulations.modeling import community
from simulations.modeling.adapter import (
    apply_genotype,
    apply_measurements,
    apply_medium,
)
from simulations.modeling.operations import apply_operations
from simulations.modeling.simulations import simulate
from simulations.schemas import (
    CommunitySimulationRequest,
    ModificationRequest,
    SimulationRequest,
)


logger = logging.getLogger(__name__)


def init_app(app):
    """Register API resources on the provided Flask application."""
    app.add_url_rule("/healthz", view_func=healthz)
    app.add_url_rule("/metrics", view_func=metrics)

    app.add_url_rule(
        "/models/<int:model_id>/modify", view_func=model_modify, methods=["POST"]
    )
    app.add_url_rule("/simulate", view_func=model_simulate, methods=["POST"])
    app.add_url_rule(
        "/community/simulate", view_func=model_community_simulate, methods=["POST"]
    )

    docs = FlaskApiSpec(app)
    docs.register(model_modify, endpoint=model_modify.__name__)
    docs.register(model_simulate, endpoint=model_simulate.__name__)
    docs.register(model_community_simulate, endpoint=model_community_simulate.__name__)


@use_kwargs(ModificationRequest)
def model_modify(
    model_id,
    medium,
    genotype,
    fluxomics,
    metabolomics,
    proteomics,
    uptake_secretion_rates,
    molar_yields,
    growth_rate,
):
    if not request.is_json:
        abort(415, "Non-JSON request content is not supported")

    try:
        model_wrapper = storage.get(model_id)
    except Unauthorized as error:
        abort(401, error.message)  # noqa: B306
    except Forbidden as error:
        abort(403, error.message)  # noqa: B306
    except ModelNotFound as error:
        abort(404, error.message)  # noqa: B306

    # Use the context manager to undo all modifications to the shared model instance on
    # completion.
    with model_wrapper.model as model:
        # Build list of operations to perform on the model
        operations = []
        warnings = []
        errors = []
        if medium:
            results = apply_medium(model, model_wrapper.is_ec_model, medium)
            operations.extend(results[0])
            warnings.extend(results[1])
            errors.extend(results[2])

        if genotype:
            results = apply_genotype(model, genotype)
            operations.extend(results[0])
            warnings.extend(results[1])
            errors.extend(results[2])

        if (
            fluxomics
            or metabolomics
            or proteomics
            or uptake_secretion_rates
            or molar_yields
            or growth_rate
        ):
            results = apply_measurements(
                model,
                model_wrapper.biomass_reaction,
                model_wrapper.is_ec_model,
                fluxomics,
                metabolomics,
                proteomics,
                uptake_secretion_rates,
                molar_yields,
                growth_rate,
            )
            operations.extend(results[0])
            warnings.extend(results[1])
            errors.extend(results[2])

        if errors:
            # If any errors occured during modifications, discard generated operations
            # and return the error messages to the client for follow-up
            return {"errors": errors}, 400
        else:
            return {"operations": operations, "warnings": warnings}


@use_kwargs(SimulationRequest)
def model_simulate(model_id, method, objective_id, objective_direction, operations):
    try:
        model_wrapper = storage.get(model_id)
    except Unauthorized as error:
        abort(401, error.message)  # noqa: B306
    except Forbidden as error:
        abort(403, error.message)  # noqa: B306
    except ModelNotFound as error:
        abort(404, error.message)  # noqa: B306

    model = model_wrapper.model

    # Use the context manager to undo all modifications to the shared model instance on
    # completion.
    with model:
        apply_operations(model, operations)
        try:
            flux_distribution, growth_rate = simulate(
                model,
                model_wrapper.biomass_reaction,
                method,
                objective_id,
                objective_direction,
            )
        except OptimizationError:
            return jsonify({"status": model.solver.status})
        else:
            return jsonify(
                {
                    "status": model.solver.status,
                    "flux_distribution": flux_distribution,
                    "growth_rate": growth_rate,
                }
            )


@use_kwargs(CommunitySimulationRequest)
def model_community_simulate(model_ids, medium, method):
    try:
        model_wrappers = [storage.get(model_id) for model_id in model_ids]
    except Unauthorized as error:
        abort(401, error.message)  # noqa: B306
    except Forbidden as error:
        abort(403, error.message)  # noqa: B306
    except ModelNotFound as error:
        abort(404, error.message)  # noqa: B306

    return community.simulate(model_wrappers, medium, method)


def metrics():
    return Response(
        generate_latest(MultiProcessCollector(CollectorRegistry())),
        mimetype=CONTENT_TYPE_LATEST,
    )


def healthz():
    """
    HTTP endpoint for readiness probes.

    Return an empty response. This response will not be ready until the application has
    finished initializing, e.g., preloading models, which takes a few minutes.
    """
    return ""
