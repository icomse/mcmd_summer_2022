import gsd
import gsd.hoomd
import numpy as np
def arrays_from_gsd(filename, keys= ['timestep','potential_energy']):
    with gsd.hoomd.open(filename,'rb') as traj:
        arrays = []
        for k in keys:
            a = []
            for frame in traj:
                if k == 'timestep':
                    a.append(frame.configuration.step)
                else:
                    a.append(frame.log['md/compute/ThermodynamicQuantities/'+k][0])
            arrays.append(np.array(a))
        return arrays  