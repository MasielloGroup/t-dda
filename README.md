# Calculation of Steady-State Temperature Distributions
This code calculates the steady-state temperature distribution of nanoparticles on a semi-infinite substrate, embedded in a semi-infinite background medium by solving the discretized steady-state heat diffusion equation.

### Step 1: Installation & set-up
* Compile the script `Lattice_Diffusion.c` by typing `make all`
* To run `t-dda`, you will need a look-up table for the values of the lattice Green's function evaluated at all points in the temperature calculation. Inside the folder `lattice_greenfunction` is a premade 20x20x20 grid, `Green_grid20.txt` and a larger 100x100x100 grid. You **MUST** make sure your entire discretized shape can fit within the extent of the grid. To make your own, see the Mathematica script `lattice_greenfunction/integrate_bessels.nb`. A 20x20x20 grid should run in a few minutes. (I am working on finding a way to perform this integration in Python for those who do not have access to Mathematica).

### Step 2: Preparing input files
Two input files are needed to calculate temperature profiles, a parameter file and an input file. The parameter file, called `var.par` in the examples, specifies system parameters. The input file, called `tdda_input` in the examples, can take on two forms. It can either specify the absorption cross secctions or the electric fields and polarizations. The user can choose which input they'd prefer to provide. Each example is provided in the `example_sphere` folder.

### Step 3: Running `t-dda`
Once the input files are prepared and the lattice Green's function values are tabulated, you are ready to run the code. Please note the following:
* All particles must start at lattice site *x* = -1 and continue into the negative *x* plane. This is because the semi-infinite substrate has been hard-coded to start at *x* = 0 and continue into the positive plane.
* If you are a new user of `DDSCAT` and/or do not wish to edit Draine's `DDSCAT` to produce the correct forms of the input files, you can use a version of `DDSCAT` we have already modified. This version can be found in [this repo](http://github.com/MasielloGroup/g-dda).

### Step 4: Citing `t-dda`
Please cite [this paper](https://pubs.acs.org/doi/10.1021/jz500421z) when using this code. For some more (informal) notes, click [here](https://www.overleaf.com/read/mrdzxbrwspqt).

## Examples
The folder `example_sphere` includes two examples for how to run the `t-dda` given two different kinds of input files. One approach uses electric fields and polarizations as inputs (`input_fp`), and the other approach uses the absorption cross-section as an input (`input_acs`). In both folders, there contains a Jupyter notebook which walks through an example calculation for the temperature rise resulting from a 10 nm radius sphere driven with plane wave light. Additionally, `example_sphere` contains a folder (`analytic_sphere`) which has the files necessary to calculate the temperature rise from a sphere analytically. Lastly, `cluster_scripts` contains files which may be used to submit these calculations on a supercomputer.

### Parameter File
The following examples will use different input files, but identical parameter files. Named `var.par` in the examples, this parameter file must be set as follows.
>num_k: 1 				number of materials included in input file
>k_out: 0.3 			thermal conductivity of semi-infinite background [W/mK]
>k_in: 314 				thermal conductivity of material in input file [W/mK]
>k_sub: 0.3 			thermal conductivity of semi-infinite substrate [W/mK]
>lambda: 0.505e-06 		wavelength corresponding to input file
>n_m: 1. 				refractive index of semi-infinte background
>I_0: 1e+9 				intensity of laser light
>unit: 1.0 				lattice spacing

>d: 1 					step size
>x_min: -25 			minumum x value
>x_max: 1 				maximum x value
>y_min: -13 			minumum y value
>y_max: 13 				maximum y value
>z_min: -13 			minumum z value
>z_max: 13 				maximum z value

>x_plane: -12 			x plane to calculate temperatures in
>input_mode: 1 			fields / polarizations or abs. cross sections

### Temperatures from field and polarization inputs (`input_fp`)
This method uses the electric fields and polarizations of each discrete dipole to form a power absorbed per particle. A Jupyter notebook which goes through an example of calculating the temperatures using this approach is found in the folder `example_sphere/input_fp`.  To run the calculation, you must have installed `g-dda`, linked above. If you have your own version of `DDSCAT`, feel free to use that and edit the scripts. 
* Starting from the optical calculation, make your shape file and edit the parameter file, `ddscat.par` according to the instructions in the g-dda repo.
* Follow the Jupyter notebook for a walk-through. (As of now, the notebook will create your shape, but will not edit `ddscat.par`)
* As you're playing around with this script, it will generate a lot output files. Feel free to `bash quick_remove.sh` to quickly delete all the output files to start over. 

### Temperatures from absorption cross section input (`input_acs`)
This method works by turning the absorption cross-section into a power absorbed per discrte dipole. This method will only work for temperature calculations on particles of identical size. You may use it to study systems with many particles, they must just be of the same size. This example is found in the folder `example_sphere/input_acs`.
* Follow the Jupyter notebook for a walk-through.
* Similar to above, make sure you edit `ddscat.par` before running the code.

### Analytic Sphere Comparison
 In `example_sphere/analytic_sphere`, there is a Jupyter notebook which compares `t-dda` to the analytic solution. For this notebook to run, make sure you run and therefore create the files as described in `input_fp`.
