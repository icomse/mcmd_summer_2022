import hoomd
import sys

trajfile = sys.argv[1]+'.gsd'
T = float(sys.argv[2])
cpu = hoomd.device.CPU()
sim = hoomd.Simulation(device=cpu,seed=0)
sim.create_state_from_gsd(filename='random.gsd') #N and V are set here
integrator = hoomd.md.Integrator(dt=0.005)
cell = hoomd.md.nlist.Cell(buffer = 0.4)
lj_potential = hoomd.md.pair.LJ(nlist=cell)
lj_potential.params[('A','A')] = dict(epsilon=1,sigma=1)
lj_potential.r_cut[('A','A')]=2.5
nvt = hoomd.md.methods.NVT(kT=T,filter=hoomd.filter.All(),tau=1.)
integrator.forces.append(lj_potential)
integrator.methods.append(nvt)
sim.operations.integrator = integrator
sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=1.5)
selection = hoomd.filter.All() # "which atoms"
logger = hoomd.logging.Logger() # will be used for "what's logged"
writer = hoomd.write.GSD(filename=trajfile, # "where to store"
                             trigger=hoomd.trigger.Periodic(100000), #when to store
                             mode='wb',
                             filter=selection) #filter=hoomd.filter.Null() to only store log
thermo_props = hoomd.md.compute.ThermodynamicQuantities(filter=selection) # What to store
logger.add(thermo_props)
logger.add(sim,quantities=['timestep','walltime','tps'])
writer.log = logger #need to tell our write which logger to use when it's logging info
sim.operations.computes.append(thermo_props) #tell our simulation to *compute* the thermo properties
sim.operations.writers.append(writer) # tell our simulation which writer(s) to use
sim.run(1e6)
