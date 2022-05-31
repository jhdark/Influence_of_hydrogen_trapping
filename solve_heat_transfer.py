from fenics import *
import FESTIM
import parameters

my_sim = FESTIM.Simulation(parameters.parameters, log_level=20)
my_sim.initialise()

T = my_sim.T
XDMFFile("Results/T.xdmf").write(T)
