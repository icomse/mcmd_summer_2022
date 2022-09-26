The following tests are implemented in this notebooks:
 
1. Tests of kinetic energy distribution
2. Tests of whether the ensemble is correct (NVT and NPT)
3. A test of equipartition of kinetic energies.

Full documentation:
-------------------
Please see https://physical-validation.readthedocs.io for the full documentation of the validation package.

You can also read our publications and cite our work here:

> Merz PT, Shirts MR. Testing for physical validity in molecular simulations. PLoS ONE 13(9): e0202764, (2018) https://doi.org/10.1371/journal.pone.0202764
> PT Merz, WT Hsu, MW Thompson, S Boothroyd, CC Walker, MR Shirts. physical_validation: A Python package to assess the physical validity of molecular simulation results. Journal of Open Source Software 7 (69), 3981, (2022) https://doi.org/10.21105/joss.03981

Directory structure:
--------------------

Top level:

  1. README.md (this file)
  2. [Python notebook](physical_validation_of_simulations.ipynb)
  3. [Slides of lecture](PhysicalValidation_Day2_July2022.pdf)

__Water Simulation__

	systems/
		water/
			ana_gromacs.py
			ana_flatfile.py
			ana_pyarray.py
			top/
			md_NVT_be_1/
			md_NVT_vr_1/
			md_NVT_be_2/
			md_NVT_vr_2/

			md_NPT_be_1/
			md_NPT_vr_1/
			md_NPT_be_2/
			md_NPT_vr_2/
			md_NPT_be_3/
			md_NPT_vr_3/
			md_NPT_be_4/
			md_NPT_vr_4/

The `md_*` folders contain results of GROMACS simulations for 900 water molecules:

 1. Ran in different ensembles:
    * NVT
    * NPT
 2. Using different thermostats: 
    * _vr_ - velocity-rescale
    * _be_ - Berendsen thermostat
 3. At different state points:
    * _NVT_ - 1) T; 2) T+dT 
    * _NPT_ - 1) (T, P); 2) (T+dT, P); 3) (T, P+dP); 4) (T+dT, P+dP) 

The `top/` folder contains the topology of the system. 

`python ana.py` performs the ensemble checks. The Python notebooks implements the same tests.

`python ana_flatfile.py` performs the same NVT ensemble checks on the velocity rescale data set, with the data loaded from an ASCII files with a list of energies. The Python notebooks implements the same tests.

`python ana_numpy.py` performs the same NVT ensemble checks on the velocity rescale data set, with the data loaded from a numpy array of energies. The Python notebooks implements the same tests.
 

__Trp-cage Protein__

	systems/
		trp-cage/
			ana.py
			be_1/
			be_2/
			vr_1/
			vr_2/

`python ana.py` performs the equipartition check on a Trp-cage protein. 

The `be_*` and `vr_*` folders contain results of simulations ran in NPT using different thermostatting schemes: 'vr' stands for velocity-rescale, 'be' for Berendsen thermostat, while '_1' denotes a simulation using a single thermostat for the entire system, and '_2' denotes a simulation using two separate thermostats, one for the protein and one for the solvent. 

Note that the solvent was stripped from the trajectory files to considerably reduce file size and execution time for the workshop, so only the solute kinetic energy is analyzed.  

[] (`tot` is the total kinetic energy, `tra` is translational kinetic energy of the group, `rni` is the rotational plus internal, `rot` is the rotational kinetic energy, and `int` is the interna kinetic energy. The Python notebooks implement the same tests.)


Some output notation used by the validation tools:
--------------------------------------------------
Before analyzing the distributions, the trajectories are analyzed to obtain uncorrelated samples and remove unequilibrated data. This consists of three steps:
  1. Equilibration: Any data points at the beginning of the trajectory which are not yet equilibrated (simulation has reached steady state) are ignored.
  2. Decorrelation: After equilibration, the remaining data points are analyzed for correlations in time. If samples are found to be correlated, only a decorrelated subset is used.
  3. Tail pruning (optional): To avoid numerical noise to strongly influence the fitting of small data sets, the extreme values of the simulation might be discarded.
  The result of this procedure is reported as `After equilibration, decorrelation and tail pruning, XX% (YYY frames) of original trajectory remain.`

