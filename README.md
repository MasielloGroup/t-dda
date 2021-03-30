# Calculation of steady-state temperature profile
This code calculates the steady-state temperature distribution of nanoparticles on a semi-infinite substrate, embedded in a semi-infinite background medium by solving the discretized steady-state heat diffusion equation.

### Step 1: Installation & Set-up
* Compile the script `Lattice_Diffusion.c` by typing `make all`
* To run t-DDA, you will need a look-up table for the values of the lattice Green function evaluated at all points in the temperature calculation. Inside the folder `lattice_greenfunction` is a premade 20x20x20 grid, `Green.num_20x20x20`. You must make sure your entire discretized shape can fit within the extent of the grid. (There is no error message yet for this mistake. It will cause the code to run indefinetly.) To make your own, see the python script `lattice_greenfunction/make_Green.py`. It is not well documented/commented  yet, and it takes a very long time to run. The 20x20x20 grid took approximately 20 hours to run. A  300x300x300 grid can be found at [this link](https://drive.google.com/file/d/1vtJJoeluL_DfTk8wFYkrMUdBL-KFVxul/view?usp=sharing).

### Step 2: Preparing input files
Two input files are needed to calculate temperature profiles, a parameter file and an input file. The parameter file, called `var.par` in the examples, specifies system parameters. The input file, called `tdda_input` in the examples, can take on two forms. It can either specify the absorption cross secction or the electric fields and polarizations. The user can choose which input they'd prefer to provide. Each example is provided in the `example_sphere` folder.

### Step 3: Running t-DDA
Once the input files are prepared and the lattice Green function values are tabulated, you are ready to run the code.

### Step 4: Citing t-DDA
Please cite [this paper](https://pubs.acs.org/doi/10.1021/jz500421z) when using this code. For some more (informal) notes, click [here](https://www.overleaf.com/read/mrdzxbrwspqt).

## Examples
The folder `example_sphere` includes two examples for how to run the t-DDA given two different input files. Both examples calculate the temperature resulting from a 5 nm radius sphere driven with plane wave light. You can verify the numerical results against the analytic solutions for a sphere given in the `analytic_sphere` folder. If you are a new user of DDSCAT and/or do not wish to edit Draine's DDSCAT to produce the correct forms of the input files, you can use a version of DDSCAT we have already modified. This version can be found in [this repo](https://github.com/MasielloGroup/g-dda).

### Temperatures from field and polarization inputs
This method is most appropriate for calculating the resulting thermal profiles of a nanoparticle system composed of multiple different shaped particles. It uses the electric fields and polarizations of each discrete dipole in a particle to form a power absorbed per dipole. This example is found in the folder `input_fp`. This folder assumes you have installed `g-dda`, linked above. If you have your own version of DDSCAT, feel free to use that and edit the scripts. 
* Starting from the optical calculation, make your shape file and edit the parameter file, `ddscat.par` according to the instructions in the g-dda repo. Note: the particle must start at lattice cite x = -1 and continue into the negative x plane. The substrate has been hard-coded to start at x = 0 and continue into the positive plane. 
* Follow the Jupyter notebook for a walk-through.
* As you're playing around with this script, it will generate a lot output files. Feel free to `bash quick_remove.sh` to quickly delete all the output files to start over. 

### Temperatures from absorption cross section input
This method works by turning the absorption cross-section into a power absorbed per discrete dipole. This method will only work for single-particle temperature calculations. This example is found in the folder `input_acs`.
* Follow the Jupyter notebook for a walk-through.