import numpy as np
from sympy import Symbol, Function, Number, pi, sin

import modulus
from modulus.hydra import to_absolute_path, ModulusConfig
from modulus.solver import Solver
from modulus.domain import Domain
from modulus.geometry.primitives_1d import Point1D, Line1D
from modulus.domain.constraint import (
    PointwiseBoundaryConstraint,
    PointwiseInteriorConstraint,
)
from modulus.domain.inferencer import PointwiseInferencer
from modulus.key import Key
from modulus.node import Node
from modulus.models.fully_connected import FullyConnectedArch
from modulus.eq.pde import PDE


class CustomPDE(PDE):
    def __init__(self, f=1.0):
        # coordinates
        x = Symbol("x")

        # make input variables
        input_variables = {"x": x}

        # make u function
        u = Function("u")(*input_variables)

        # source term
        if type(f) is str:
            f = Function(f)(*input_variables)
        elif type(f) in [float, int]:
            f = Number(f)

        # set equations
        self.equations = {}
        self.equations["custom_pde"] = (
            u.diff(x, 2) - f
        )  # "custom_pde" key name will be used in constraints


@modulus.main(config_path="conf", config_name="config")
def run(cfg: ModulusConfig) -> None:

    # make list of nodes to unroll graph on
    eq = CustomPDE(f=1.0)
    u_net = FullyConnectedArch(
        input_keys=[Key("x"), Key("l")], output_keys=[Key("u")], nr_layers=3, layer_size=32
    )

    nodes = eq.make_nodes() + [u_net.make_node(name="u_network")]

    # add constraints to solver
    # make geometry
    x, l = Symbol("x"), Symbol("l")
    geo = Line1D(0, l)

    # make domain
    domain = Domain()

    # bcs
    bc = PointwiseBoundaryConstraint(
        nodes=nodes,
        geometry=geo,
        outvar={"u": 0},
        batch_size=10,
        parameterization={l: (1, 2)},
    )
    domain.add_constraint(bc, "bc")

    # interior
    interior = PointwiseInteriorConstraint(
        nodes=nodes,
        geometry=geo,
        outvar={"custom_pde": 0},
        batch_size=100,
        bounds={x: (0, 1)},
        parameterization={l: (1, 2)},
    )
    domain.add_constraint(interior, "interior")
    
    # add inferencer
    l_inf = np.linspace(1.0, 2.0, 11)
    for i in range(len(l_inf)):
        x_inf = np.linspace(0, l_inf[i], 100)
        inference = PointwiseInferencer(
            nodes=nodes,
            invar={"x": x_inf.reshape(-1,1), "l": l_inf[i]*np.ones_like(x_inf.reshape(-1, 1))},
            output_names=["u"],
        )
        domain.add_inferencer(inference, "inf_data_"+str(i))

    # make solver
    slv = Solver(cfg, domain)

    # start solver
    slv.solve()


if __name__ == "__main__":
    run()
