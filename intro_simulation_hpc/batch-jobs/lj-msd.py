#Intended for getting TPS and density data in batch slurm jobs as a function of T, P, -n , N_particles

import hoomd
import gsd.hoomd
import itertools
import numpy as np
import math
import freud
import sys

#key variables
m = 7 #increase for more atoms
N_particles = 4 * m**3 #helper for initialization
Temperature = float(sys.argv[2])

tau = 0.2
trajfile = 'nvt'+sys.argv[1]+'.gsd'
write_period = 1e4 
maxtime = 1e6

############################
# Turn on HOOMD and initialize configuration
gpu = hoomd.device.GPU()
sim = hoomd.Simulation(device=gpu,seed=0)

spacing = 1.3
K = math.ceil(N_particles**(1 / 3))
L = K * spacing
x = np.linspace(-L / 2, L / 2, K, endpoint=False)
position = list(itertools.product(x, repeat=3))

snapshot = gsd.hoomd.Snapshot()
snapshot.particles.N = N_particles
snapshot.particles.position = position[0:N_particles]
snapshot.particles.typeid = [0] * N_particles
snapshot.configuration.box = [L, L, L, 0, 0, 0]
snapshot.particles.types = ['C']
sim.create_state_from_snapshot(snapshot) #may need debugging

######################
#Potential and integrator setup
integrator = hoomd.md.Integrator(dt=0.005)
cell = hoomd.md.nlist.Cell(buffer = 0.4)
lj_potential = hoomd.md.pair.LJ(nlist=cell)
lj_potential.params[('C','C')] = dict(epsilon=1,sigma=1)
lj_potential.r_cut[('C','C')]=2.5
ensemble = hoomd.md.methods.NVT(kT=Temperature,filter=hoomd.filter.All(),tau=tau)
integrator.forces.append(lj_potential)
integrator.methods.append(ensemble)
sim.operations.integrator = integrator
sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=Temperature)

######################
# Logger setup and run
selection = hoomd.filter.All() # "which atoms"
logger = hoomd.logging.Logger() # will be used for "what's logged"
writer = hoomd.write.GSD(filename=trajfile, # "where to store"
                             trigger=hoomd.trigger.Periodic(int(write_period)), #when to store
                             mode='wb',
                             filter=selection) #filter=hoomd.filter.Null() to only store log
thermo_props = hoomd.md.compute.ThermodynamicQuantities(filter=selection) # What to store
logger.add(thermo_props)
logger.add(sim,quantities=['timestep','walltime','tps'])
writer.log = logger #need to tell our write which logger to use when it's logging info
sim.operations.computes.append(thermo_props) #tell our simulation to *compute* the thermo properties
sim.operations.writers.append(writer) # tell our simulation which writer(s) to use
sim.run(maxtime)

######################
#analyze
step = []
volume= []
tps = []
N = 0

with gsd.hoomd.open(trajfile,'rb') as traj:
    for frame in traj:
        step.append(frame.configuration.step)
        tps.append(frame.log['Simulation/tps'][0])
        volume.append(frame.log['md/compute/ThermodynamicQuantities/volume'][0])
step, tps, volume= np.array(step), np.array(tps), np.array(volume)
print("N={} T={}, V={}: TPS={}".format(N_particles, Temperature, volume.mean(),  tps.mean()) )

##############################
# MSD analysis
start = 0
stop = -1
with gsd.hoomd.open(trajfile, "rb") as trajectory:
    init_box = trajectory[start].configuration.box
    final_box = trajectory[stop].configuration.box
    assert all(
            [i == j for i,j in zip(init_box, final_box)]
            ), f"The box is not consistent over the range {start}:{stop}"
    positions = []
    images = []
    for frame in trajectory[start:stop]:
        atom_pos = frame.particles.position[:]
        atom_img = frame.particles.image[:]
        positions.append(atom_pos)
        images.append(atom_img)
    if np.count_nonzero(np.array(images)) == 0:
        warn(
            f"All of the images over the range {start}-{stop} "
            "are [0,0,0]. You may want to ensure this gsd file "
            "had the particle images written to it."
            )
    msd = freud.msd.MSD(box=init_box, mode='window')
    msd.compute(np.array(positions), np.array(images), reset=False)
    np.save("msd.npy", msd.particle_msd)

