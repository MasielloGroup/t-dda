# Calculation of steady-state temperature profile
This code calculates the steady-state temperature distribution of nanoparticles on a semi-infinite substrate, embedded in a semi-infinite background medium by solving the discretized steady-state heat diffusion equation.

## Step 1: Installation & Set-up
* Complie the script "Lattice_Diffusion.c" by typing "make all"
* To run t-DDA, you will need a look-up table for the values of the lattice green function evaluated at all points in the temperature calculation. Inside the folder "lattice_greenfunction" is a premade 20x20x20 grid. A much larger 300x300x300 grid can be found at: https://drive.google.com/file/d/1s4z__Sj6Ze3STPznusq1IBdATIJjLmko/view?usp=sharing. To make your own, see the python script "make_Green.py"

## Step 2: Preparing input files
Two input files are needed to calculate temperature profiles
* var.par: Specifies the spatial extent and other system parameters.
* tdda_input: Specifies the electric fields and polarizations for every (absorbing) point

For specific requirements of these files, see https://www.overleaf.com/read/mrdzxbrwspqt.

## Step 3: Citing t-DDA
Please cite https://pubs.acs.org/doi/10.1021/jz500421z when using this code. For some more (informal) notes, see: https://www.overleaf.com/read/mrdzxbrwspqt.

