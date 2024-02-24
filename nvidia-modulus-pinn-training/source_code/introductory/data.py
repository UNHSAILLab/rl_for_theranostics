from sympy import Symbol
import numpy as np

import modulus
from modulus.hydra import to_absolute_path, ModulusConfig
from modulus.solver import Solver
from modulus.domain import Domain
from modulus.geometry.primitives_1d import Point1D, Line1D
from modulus.domain.constraint import (
    PointwiseConstraint,
)
from modulus.domain.inferencer import PointwiseInferencer
from modulus.key import Key
from modulus.node import Node
from modulus.models.fully_connected import FullyConnectedArch

@modulus.main(config_path="conf", config_name="config")
def run(cfg: ModulusConfig) -> None:

    # make list of nodes to unroll graph on
    u_net = FullyConnectedArch(
        input_keys=[Key("x")], output_keys=[Key("u")], nr_layers=3, layer_size=32
    )

    nodes = [u_net.make_node(name="u_network")]

    # add constraints to solver
    # make geometry
    x = Symbol("x")
    geo = Line1D(0, 1)

    # make domain
    domain = Domain()

    # data
    x_np = np.linspace(0, 0.3, 4)
    u_np = 0.5*(x_np-1)*x_np
    data = PointwiseConstraint.from_numpy(
        nodes=nodes,
        invar={"x": x_np.reshape(-1,1)},
        outvar={"u": u_np.reshape(-1,1)},
        batch_size=4
    )
    domain.add_constraint(data, "data")

    # add inferencer
    inference = PointwiseInferencer(
        nodes=nodes,
        invar={"x": np.linspace(0, 1.0, 100).reshape(-1,1)},
        output_names=["u"],
    )
    domain.add_inferencer(inference, "inf_data")

    # make solver
    slv = Solver(cfg, domain)

    # start solver
    slv.solve()


if __name__ == "__main__":
    run()
