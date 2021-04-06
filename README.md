# Calculation of steady-state temperature profile
This code calculates the steady-state temperature distribution of nanoparticles on a semi-infinite substrate, embedded in a semi-infinite background medium by solving the discretized steady-state heat diffusion equation.

### Step 1: Installation & Set-up
* Compile the script `Lattice_Diffusion.c` by typing `make all`
* To run t-DDA, you will need a look-up table for the values of the lattice Green function evaluated at all points in the temperature calculation. Inside the folder `lattice_greenfunction` is a premade 20x20x20 grid, `Green_grid20.txt`. You must make sure your entire discretized shape can fit within the extent of the grid. To make your own, see the Mathematica script `lattice_greenfunction/integrate_bessels.nb`. A 20x20x20 grid should run in a few minutes. (I am working on finding a way to perform this integration in python for those who do not have access to Mathematica).

### Step 2: Preparing input files
Two input files are needed to calculate temperature profiles, a parameter file and an input file. The parameter file, called `var.par` in the examples, specifies system parameters. The input file, called `tdda_input` in the examples, can take on two forms. It can either specify the absorption cross secction or the electric fields and polarizations. The user can choose which input they'd prefer to provide. Each example is provided in the `example_sphere` folder.

### Step 3: Running `t-dda`
Once the input files are prepared and the lattice Green function values are tabulated, you are ready to run the code. Please note the following:
* All particles must start at lattice site x = -1 and continue into the negative x plane. This is because the semi-infinite substrate has been hard-coded to start at x = 0 and continue into the positive plane.
* If you are a new user of DDSCAT and/or do not wish to edit Draine's DDSCAT to produce the correct forms of the input files, you can use a version of DDSCAT we have already modified. This version can be found in [this repo](http://github.com/MasielloGroup/g-dda).

### Step 4: Citing t-DDA
Please cite [this paper](https://pubs.acs.org/doi/10.1021/jz500421z) when using this code. For some more (informal) notes, click [here](https://www.overleaf.com/read/mrdzxbrwspqt).

## Examples
The folder `example_sphere` includes two examples for how to run the `t-dda` given two different input files. Both examples calculate the temperature resulting from a 5 nm radius sphere driven with plane wave light. (Soon you will be able to verify the numerical results against the analytic solutions for a sphere given in the not yet created `analytic_sphere` folder.) 

### Temperatures from field and polarization inputs
This method is most appropriate for calculating the resulting thermal profiles of a nanoparticle system composed of multiple different shaped particles. It uses the electric fields and polarizations of each discrete dipole in a particle to form a power absorbed per dipole. This example is found in the folder `input_fp`. This folder assumes you have installed `g-dda`, linked above. If you have your own version of DDSCAT, feel free to use that and edit the scripts. 
* Starting from the optical calculation, make your shape file and edit the parameter file, `ddscat.par` according to the instructions in the g-dda repo.
* Follow the Jupyter notebook for a walk-through. (As of now, the notebook will create your shape, but will not edit ddscat.par)
* As you're playing around with this script, it will generate a lot output files. Feel free to `bash quick_remove.sh` to quickly delete all the output files to start over. 

### Temperatures from absorption cross section input
This method works by turning the absorption cross-section into a power absorbed per discrete dipole. This method will only work for single-particle temperature calculations. This example is found in the folder `input_acs`.
* Follow the Jupyter notebook for a walk-through.
* Similar to above, make sure you edit `ddscat.par` before running the code.