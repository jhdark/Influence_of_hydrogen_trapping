from fenics import *
import FESTIM
from properties import S_0_lipb, E_S_lipb

from parameters import parameters, id_lipb


def run_H_transport(parameters, S_0_lipb=S_0_lipb, E_S_lipb=E_S_lipb, log_level=20):
    """Runs a hydrogen transport simulation projecting a predefined velcoity field onto
    a given mesh, accounting for advection

    Args:
        parameters ([type]): [description]
        S_0_lipb (float, optional): pre-exponential factor of solubility in
            LiPb. Defaults to S_0_lipb.
        E_S_lipb (float, optional): activation energy of solubility in LiPb.
            Defaults to E_S_lipb.
        log_level (int, optional): [description]. Defaults to 20.

    Returns:
        [type]: [description]
    """
    # create a simulation with these normal parameters
    my_sim = FESTIM.Simulation(parameters, log_level=log_level)
    my_sim.initialise()

    # read the u_full function written in solve_NS_submesh.py
    mesh = my_sim.mesh
    V_ele = VectorElement("CG", mesh.ufl_cell(), 2)
    V_u = FunctionSpace(mesh, V_ele)

    mesh_cfd = "Results/velocity_fields/u_full_fine.xdmf"
    u = Function(V_u)
    XDMFFile(mesh_cfd).read_checkpoint(u, "u", -1)

    # modify the form F
    id_flow = id_lipb
    test_function_solute = split(my_sim.v)[0]
    solute = split(my_sim.u)[0]
    dx = my_sim.dx
    S_lipb = S_0_lipb * exp(-E_S_lipb / FESTIM.k_B / my_sim.T)
    my_sim.F += inner(dot(grad(S_lipb * solute), u), test_function_solute) * dx(id_flow)

    # run the simulation with the modified formulation
    return my_sim.run()


if __name__ == "__main__":
    run_H_transport(parameters)
