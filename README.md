### Calculation of steady-state temperature profile

# Step 1: Installation & Set-up
* Complie the script "Lattice_Diffusion.c" by typing "make all"
* To run t-dda, you will need a look-up table for the values of the lattice green function evaluated at all points in the temperature calculation. Inside the folder "lattice_greenfunction" is a premade 20x20x20 grid. A much larger 300x300x300 grid can be found at: https://drive.google.com/file/d/1s4z__Sj6Ze3STPznusq1IBdATIJjLmko/view?usp=sharing. To make your own, see the python script "make_Green.py"
