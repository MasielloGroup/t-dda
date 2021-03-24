# Calculation of steady-state temperature profile
This code calculates the steady-state temperature distribution of nanoparticles on a semi-infinite substrate, embedded in a semi-infinite background medium by solving the discretized steady-state heat diffusion equation.

## Step 1: Installation & Set-up
* Compile the script `Lattice_Diffusion.c` by typing `make all`
* To run t-DDA, you will need a look-up table for the values of the lattice Green function evaluated at all points in the temperature calculation. Inside the folder `lattice_greenfunction` is a premade 20x20x20 grid. A much larger 300x300x300 grid can be found at: https://drive.google.com/file/d/1s4z__Sj6Ze3STPznusq1IBdATIJjLmko/view?usp=sharing. To make your own, see the python script "make_Green.py". You must make sure your entire discretized shape can fit within the extent of the grid.

## Step 2: Preparing input files
Two input files are needed to calculate temperature profiles, a parameter file and an input file. The parameter file, called `var.par` in the examples, specifies system parameters. The input file, called `tdda_input` in the examples, can take on two forms. It can either specify the absorption cross secction or the electric fields and polarizations. The user can choose which input they'd prefer to provide. Each example is provided in the `example_sphere` folder. 
* var.par: Specifies the spatial extent and other system parameters.
* tdda_input: Specifies either the electric fields and polarizations for every (absorbing) point, or the heat powers for every (absorbing point. 

## Step 3: Running t-DDA
Once the input files are prepared and the lattice Green function values are tabulated, type `bash run_tdda.sh` to run the code. 

## Step 4: Citing t-DDA
Please cite [this paper](https://pubs.acs.org/doi/10.1021/jz500421z) when using this code. For some more (informal) notes, see: https://www.overleaf.com/read/mrdzxbrwspqt.

# Examples
The folder `example_sphere` includes two examples for how to run the t-DDA given two different input files. Both examples calculate the temperature resulting from a 5 nm radius sphere driven with plane wave light. You can verify the numerical results against the analytic solutions for a sphere given in the `analytic_sphere` folder. If you are a new user of DDSCAT and/or do not wish to edit Draine's DDSCAT to produce the correct forms of the input files, you can use a version of DDSCAT we have already modified. This version can be found in (this repo)(https://github.com/MasielloGroup/g-dda).

## Temperatures from field and polarization inputs
This method is most appropriate for calculating the resulting thermal profiles of a nanoparticle system composed of multiple different shaped particles. It uses the electric fields and polarizations of each discrete dipole in a particle to form a power absorbed per dipole. DDSCAT can be used to generate these values. If you do not wish to take the time to write scripts to format DDSCAT to do this, feel free to use a version of DDSCAT we have edited to do this. This version can be found at: 


## Temperatures from absorption cross section input