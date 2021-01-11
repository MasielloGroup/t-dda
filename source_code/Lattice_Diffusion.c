#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "lattice_prototypes.h"

//	This code implements the T-DDA. Author: Chris Baldwin
//	This program computes steady-state temperature distributions for an externally-illuminated target (represented by a finite list of points on a cubic lattice) immersed in a homogeneous background.
//	The syntax for running this program is as follows: 'Lattice_Diffusion <Green's function value list> <parameter list> <input file name> <output file name>'
//	What should be in each input is as follows:
//		-<Green's funion value list>: A list of the lattice Green's function values for a cubic lattice, as described in the paper "Application of the Lattice Green's Function for Calculating the Resistance of an Infinite Network of Resistors" by Jozsef Cserti, published in AJP, volume 68, pg. 896, in 2000.
//										The list should be formatted as follows: "%d %d %d %f\n", where the first three integers (which should be non-negative) give the position of each point relative to a source at the origin, and the fourth number gives the value of the lattice Green's function at that point.
//		-<parameter list>: 	A list of the parameter values to be used in the calculation. A description of the parameters is given below.
//		-<input file name>: The location of the file containing the target geometry and heat source.
//							This file should be formatted in one of two ways, depending on the value of 'input_mode' in the parameter list.
//							If input_mode == 1, then each line of this file should be of the form: "%f %f %f %d %f %f %f %f %f %f\n". The first three floats are the positions of each point in the target (in units of the lattice spacing), the next int is an index for the composition of that point, and the next six floats give the real and imaginary parts of each component of the electric field at that point (which should include the scattered electric field).
//							If input_mode != 1, then each line of this file should be of the form: "%f %f %f %f %f %f %f %f %f\n". The first three floats are the positions of each point in the target, and the next six floats are the real and imaginary parts of each component of the electric field at that point (which should include the scattered electric field).
//		-<output file name>: 	The location of the file to write the output into (existing files will be overwritten).
//								Each line of the output file is of the form: "%d %d %d %f %f %f\n", where the first three integers are the positions of each point in the target (in units of the lattice spacing), the next float is the temperature at that point, then the original heat source at that point, and finally the "effective heat source" at that point, i.e., the heat source including the thermal equivalent of "bound charges".

int main(int argc, char *argv[]) {
	if (argc != 5) {
		printf("Invalid syntax: Expected format is 'Lattice_Diffusion <Green's function value list> <parameter list> <input file name> <output file name>'\n");
		return 1;
        }else{
	        printf("You have the appropriate number of files specified\n");
	}

	char *dat_name = argv[1];
	char *par_name = argv[2];
	char *inp_name = argv[3];
	char *out_name = argv[4];

	FILE *dat_id, *par_id, *inp_id, *out_id;

	dat_id = fopen(dat_name, "r");
	par_id = fopen(par_name, "r");
	inp_id = fopen(inp_name, "r");
	out_id = fopen(out_name, "w");

	if (dat_id == NULL) {
		printf("ERROR: Invalid values file\n");
		return 0;
	}
	if (par_id == NULL) {
		printf("ERROR: Invalid param file\n");
		return 0;
	}
	if (inp_id == NULL) {
		printf("ERROR: Invalid input file\n");
		return 0;
	}
	if (out_id == NULL) {
		printf("ERROR: Invalid output file\n");
		return 0;
	}

	int i, j, k;
	
//	Reading in parameters. The parameter values are set in the "parameter list" file given as input to the program
	int num_k = 0; //Number of materials that make up the target, e.g., for a target composed of a homogeneous material set num_k = 1 in the "parameter list" file
	float *kappa; //kappa is the list of conductivities (the first element is the background conductivity)
	float n_m = 0.0, lambda = 0.0, I_0 = 0.0, unit = 0.0, x_plane = 0.0;//n_m: background medium's index of refraction, lambda is the illuminating light's wavelength, I_0 is the intensity of the illuminating beam, and unit is the unit of length in nm that's used in the calculations
	int d = 0, x_min = 0, x_max = 0, y_min = 0, y_max = 0, z_min = 0, z_max = 0; //d is the lattice spacing in units of "unit", the next six ints define the size of the grid on which calculations are done (the entire target should be included in the region defined by these values). 
	int input_mode = 1; //Set input_mode = 1 to require that the input file of target points also give an index for each point's composition, i.e., if a target point has an index of i, then that point has a thermal conductivity given by *(kappa + i). If input_mode != 1, then all target points are assumed to be the same material.
	printf("Initializing parameters\n");            
	fscanf(par_id, "num_k: %d\n", &num_k);
	kappa = (float *)malloc((num_k + 1)*sizeof(float));
	fscanf(par_id, "k_out: %f\n", kappa);
	for (i = 0; i < num_k; i++) {
	  fscanf(par_id, "k_in: %f\n", kappa + i + 1);
	  printf("%f\n", *(kappa+i));
	}
	fscanf(par_id, "k_sub: %f\n", kappa + num_k + 1);
	fscanf(par_id, "lambda: %f\nn_m: %f\n", &lambda, &n_m);
	fscanf(par_id, "I_0: %f\nunit: %f\n\n", &I_0, &unit);
	fscanf(par_id, "d: %d\nx_min: %d\nx_max: %d\ny_min: %d\ny_max: %d\nz_min: %d\nz_max: %d\n\n", &d, &x_min, &x_max, &y_min, &y_max, &z_min, &z_max);
	fscanf(par_id, "x_plane: %f\n", &x_plane);
	
	printf("%f %f %f %f %f %f %f\n", *kappa, *(kappa + 1), *(kappa + 2), *(kappa + 3) ,lambda*1e6, n_m, I_0);
	printf("%d %d %d %d %d %d %d\n", d, x_min, x_max, y_min, y_max, z_min, z_max);
	printf("%d\n", input_mode);
	
//      Calculating Substrate Ratios
	float upper = (*(kappa + num_k + 1) - *kappa)/(*(kappa + num_k + 1) + *kappa);
	float lower = (*(kappa)*2)/(*(kappa + num_k + 1) + *kappa);

	printf("%f %f\n",upper,lower);

//	Defining additional variables
	int num_x = (x_max - x_min)/d + 1;
	int num_y = (y_max - y_min)/d + 1;
	int num_z = (z_max - z_min)/d + 1;
		
//	Creating arrays of all possible x coordinates, y coordinates, z coordinates
	int *x_coords = (int *)malloc(num_x*sizeof(int));
	int *y_coords = (int *)malloc(num_y*sizeof(int));
	int *z_coords = (int *)malloc(num_z*sizeof(int));
	for (i = 0; i < num_x; i++) {
		x_coords[i] = x_min + i*d;
	}
	for (j = 0; j < num_y; j++) {
		y_coords[j] = y_min + j*d;
	}
	for (k = 0; k < num_z; k++) {
		z_coords[k] = z_min + k*d;
	}
	
//	Reading in lattice Green's function values
	int x_lim = 0, y_lim = 0, z_lim = 0;
	fscanf(dat_id, "%d %d %d\n", &x_lim, &y_lim, &z_lim);
	float *G_r = (float *)malloc(x_lim*y_lim*z_lim*sizeof(float));
	float dat_val = 0.0;
	while (!feof(dat_id)) {
		fscanf(dat_id, "%d %d %d %f\n", &i, &j, &k, &dat_val);
		*(G_r + i + x_lim*j + x_lim*y_lim*k) = dat_val;
	}
	
//	Declaring useful variables
	int N_init = 0; //Counts the number of target points in the input file
	int N_init_guess = 200000; //Initial guess as to how many target points are in the input file
	int *initial_reading = (int *)malloc(N_init_guess*sizeof(int)); //List of each target point
	int *material_ir = (int *)malloc(N_init_guess*sizeof(int)); //List of the corresponding material index of each target point
	float *Q_ir = (float *)malloc(N_init_guess*sizeof(float)); //List of the corresponding heat source applied to each target point
	int x, y, z;
	int l, index, m_val;
	size_t n, p;

//	Reading in input file
	float x_f, y_f, z_f, E_xr, E_xi, E_yr, E_yi, E_zr, E_zi, P_xr, P_xi, P_yr, P_yi, P_zr, P_zi;
	int x_tar_min = 0, x_tar_max = 0, y_tar_min = 0, y_tar_max = 0, z_tar_min = 0, z_tar_max = 0; //To speed up the rest of the calculation, we want to obtain the tightest bound on the region in which the target is contained. These bounds are given by these six ints.
	
	while (!feof(inp_id)) {
		if (input_mode == 1) {
		  fscanf(inp_id, "%f %f %f %d %f %f %f %f %f %f %f %f %f %f %f %f\n", &x_f, &y_f, &z_f, &m_val, &E_xr, &E_xi, &E_yr, &E_yi, &E_zr, &E_zi, &P_xr, &P_xi, &P_yr, &P_yi, &P_zr, &P_zi);
		} else {
			fscanf(inp_id, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &x_f, &y_f, &z_f, &E_xr, &E_xi, &E_yr, &E_yi, &E_zr, &E_zi, &P_xr, &P_xi, &P_yr, &P_yi, &P_zr, &P_zi);
			m_val = 1;
		}
		x = (x_f >= 0 ? (int)(x_f + 0.5) : (int)(x_f - 0.5));
		y = (y_f >= 0 ? (int)(y_f + 0.5) : (int)(y_f - 0.5));
		z = (z_f >= 0 ? (int)(z_f + 0.5) : (int)(z_f - 0.5));
		
		if (x < x_tar_min) x_tar_min = x; //We've found a point outside the existing bound on the region the target is in, so extend the region to contain this new point.
		if (x > x_tar_max) x_tar_max = x; //'' '' '' '' ''...
		if (y < y_tar_min) y_tar_min = y;
		if (y > y_tar_max) y_tar_max = y;
		if (z < z_tar_min) z_tar_min = z;
		if (z > z_tar_max) z_tar_max = z;
		
		i = (x - x_min)/d;
		j = (y - y_min)/d;
		k = (z - z_min)/d;
		*(initial_reading + N_init) = k + (num_z + 1)*(j + (num_y + 1)*i); //Each triplet giving a target point (i.e., (x, y, z)) is compressed into a single integer index, then stored
		*(material_ir + N_init) = m_val;
		*(Q_ir + N_init) = 4*M_PI*(2*M_PI*n_m/lambda)*I_0 * (unit*d) * (unit*d) * (unit*d) * (E_xr*P_xi - E_xi*P_xr + E_yr*P_yi - E_yi*P_yr + E_zr*P_zi - E_zi*P_zr) * (1e-18);
                printf("%f ",*(Q_ir + N_init));
                *(Q_ir + N_init) -= 8*M_PI/3.0*(2*M_PI*n_m/lambda)*(2*M_PI*n_m/lambda)*(2*M_PI*n_m/lambda)*(2*M_PI*n_m/lambda)*I_0*(unit*d)*(unit*d)*(unit*d)*(unit*d)*(unit*d)*(unit*d)*(P_xr*P_xi + P_yr*P_yi + P_zr*P_zi)*(1e-45);
                printf("%f\n",*(Q_ir + N_init));
		N_init++;
		if (N_init >= N_init_guess) { //If the arrays have run out of room, then expand the memory allocated to them
			initial_reading = (int *)realloc(initial_reading, 2*N_init_guess*sizeof(int));
			material_ir = (int *)realloc(material_ir, 2*N_init_guess*sizeof(int));
			Q_ir = (float *)realloc(Q_ir, 2*N_init_guess*sizeof(float));
			N_init_guess *= 2;
		}
	}
	
//	Sorting initial_reading (and material_ir and Q_ir) in order of increasing target point index, so that initial_reading contains entries in increasing order.
//	This makes it easier to search for specific values in the array later.
	float q_val;
	int node_num; //Node #, i.e., location in array, of minimum index
	for (n = 0; n < N_init; n++) {
		node_num = n;
		index = *(initial_reading + node_num);
		for (l = n + 1; l < N_init; l++) {
			if (*(initial_reading + l) < index) { //Found new lowest index, this becomes the new minimum
				node_num = l;
				index = *(initial_reading + node_num);
			}
		}
		m_val = *(material_ir + node_num);
		q_val = *(Q_ir + node_num);
		*(initial_reading + node_num) = *(initial_reading + n);
		*(material_ir + node_num) = *(material_ir + n);
		*(Q_ir + node_num) = *(Q_ir + n);
		*(initial_reading + n) = index; //Put lowest target point index in the next spot in array
		*(material_ir + n) = m_val;
		*(Q_ir + n) = q_val;
	}
	
//	Constructing a list of all points that need to have a modified thermal conductivity with at least one of their neighbors, i.e., points that are either inside the target or directly adjacent to points that are inside the target. This constitutes the "extended target": the original target plus a surrounding layer of background points.
//	Also constructing a list of points in the background medium at which we'd like to know the steady-state temperature.
	int N = N_init, N_out = 0;
	int N_guess = N_init_guess, N_out_guess = 5000;;
	int *target_points = (int *)malloc(N_guess*sizeof(int)); //List of indices for points that are in the extended target, i.e., either part of the target or have at least one neighbor that's part of the target
	int *outside_points = (int *)malloc(N_out_guess*sizeof(int));
	for (n = 0; n < N_init; n++) {
		*(target_points + n) = *(initial_reading + n); //All points from the initial reading of the input file are target points, so they should be in the list of extended target points
	}
	for (i = 0; i < num_x; i++) {
		x = x_coords[i];
		for (j = 0; j < num_y; j++) {
			y = y_coords[j];
			for (k = 0; k < num_z; k++) {
				z = z_coords[k];
				if (x_tar_min - 1 <= x && x <= x_tar_max + 1 && y_tar_min - 1 <= y && y <= y_tar_max + 1 && z_tar_min - 1 <= z && z <= z_tar_max + 1) { //We should check to see if the current point (given by x, y, z) is even close to the target. If not, there's no chance that this point is part of the extended target, so we shouldn't waste time checking for it.
					if (find_neighbor(initial_reading, N_init, k + (num_z + 1)*(j + (num_y + 1)*i)) < 0) { //If the current point (given by x, y, z) isn't in the list of target points from the input file...
						if (find_neighbor(initial_reading, N_init, k + (num_z + 1)*(j + (num_y + 1)*(i + 1))) >= 0 || //Adjacent point at (x + d, y, z)
							find_neighbor(initial_reading, N_init, k + (num_z + 1)*(j + (num_y + 1)*(i - 1))) >= 0 || //Adjacent point at (x - d, y, z)
							find_neighbor(initial_reading, N_init, k + (num_z + 1)*(j + 1 + (num_y + 1)*i)) >= 0 || //Adjacent point at (x, y + d, z)
							find_neighbor(initial_reading, N_init, k + (num_z + 1)*(j - 1 + (num_y + 1)*i)) >= 0 ||//Adjacent point at (x, y - d, z)
							find_neighbor(initial_reading, N_init, k + 1 + (num_z + 1)*(j + (num_y + 1)*i)) >= 0 ||
							find_neighbor(initial_reading, N_init, k - 1 + (num_z + 1)*(j + (num_y + 1)*i)) >= 0) { //And if at least one of the adjacent points is in the list of target points from the input file...
								*(target_points + N) = k + (num_z + 1)*(j + (num_y + 1)*i); //The current point is also part of the extended target
								*(material_ir + N) = 0; //The current point is technically still part of the background medium, so it has material index 0
								*(Q_ir + N) = 0.0;
								N++;
								if (N >= N_guess) { //Arrays have run out of memory, reallocate more for them
									target_points = (int *)realloc(target_points, 2*N_guess*sizeof(int));
									material_ir = (int *)realloc(material_ir, 2*N_guess*sizeof(int));
									Q_ir = (float *)realloc(Q_ir, 2*N_guess*sizeof(float));
									N_guess *= 2;
								}
						} else if (x >= x_min && x <= x_max) { //See if the current point, although far away from the target, is still one that we'd like to know the temperature at. This should be manually changed to get the temperature at different outside points.
							*(outside_points + N_out) = k + (num_z + 1)*(j + (num_y + 1)*i);
							N_out++;
							if (N_out >= N_out_guess) { //Array has run out of memory, reallocate more
								outside_points = (int *)realloc(outside_points, 2*N_out_guess*sizeof(int));
								N_out_guess *= 2;
							}
						}
					}
				} else if (x >= x_min && x <= x_max) { //See if the current point, although far away from the target, is still one that we'd like to know the temperature at. This should be manually changed to get the temperature at different outside points.
					*(outside_points + N_out) = k + (num_z + 1)*(j + (num_y + 1)*i);
					N_out++;
					if (N_out >= N_out_guess) { //Array has run out of memory, reallocate more
						outside_points = (int *)realloc(outside_points, 2*N_out_guess*sizeof(int));
						N_out_guess *= 2;
					}
				}
			}
		}
	}
	
//	Determining the composition of each neighbor to each node 
	int *neighbor_material_ir = (int *)calloc(6*N, sizeof(int)); //6xN array listing the material composition of each neighbor to each node, i.e., *(neighbor_material_ir + n + i*N) gives the material index of the (i+1)th neighbor of the nth extended target point
	int ne_index;
	for (n = 0; n < N; n++) {
		index = *(target_points + n);
		i = index/((num_z + 1)*(num_y + 1));
		j = (index % ((num_z + 1)*(num_y + 1)))/(num_z + 1);
		k = index % (num_z + 1);
		if ((ne_index = find_neighbor(initial_reading, N_init, k + (num_z + 1)*(j + (num_y + 1)*(i + 1)))) >= 0) *(neighbor_material_ir + n + 0*N) = *(material_ir + ne_index); //The neighbor at (x + d, y, z) is part of the target, so copy the neighbor's material index
		if ((ne_index = find_neighbor(initial_reading, N_init, k + (num_z + 1)*(j + (num_y + 1)*(i - 1)))) >= 0) *(neighbor_material_ir + n + 1*N) = *(material_ir + ne_index); //The neighbor at (x - d, y, z) is '' '' '' ''...
		if ((ne_index = find_neighbor(initial_reading, N_init, k + (num_z + 1)*(j + 1 + (num_y + 1)*i))) >= 0) *(neighbor_material_ir + n + 2*N) = *(material_ir + ne_index); //The neighbor at (x, y + d, z) is '' '' '' ''...
		if ((ne_index = find_neighbor(initial_reading, N_init, k + (num_z + 1)*(j - 1 + (num_y + 1)*i))) >= 0) *(neighbor_material_ir + n + 3*N) = *(material_ir + ne_index); //The neighbor at (x, y - d, z) is '' '' '' ''...
		if ((ne_index = find_neighbor(initial_reading, N_init, k + 1 + (num_z + 1)*(j + (num_y + 1)*i))) >= 0) *(neighbor_material_ir + n + 4*N) = *(material_ir + ne_index);
		if ((ne_index = find_neighbor(initial_reading, N_init, k - 1 + (num_z + 1)*(j + (num_y + 1)*i))) >= 0) *(neighbor_material_ir + n + 5*N) = *(material_ir + ne_index);
	}
	
//	Splitting the target points into boundary points and interior points: boundary points have at least one neighbor that has a different material composition, interior points have all neighbors being of the same composition as they are
	size_t N_in = 0, N_bound = 0; 
	int *interior_points = (int *)malloc(N*sizeof(int)); //List of interior points
	int *boundary_points = (int *)malloc(N*sizeof(int)); //List of boundary points
	int *material = (int *)malloc(N*sizeof(int)); //Final list of material indices (this will now be sorted so that all interior points come first)
	int *neighbor_material = (int *)malloc(6*N*sizeof(int)); //Final list of the material indices of each point's neighbors (also sorted so that all interior points come first)
	float *Q = (float *)malloc(N*sizeof(float)); //Final list of the heat source at each point (also sorted so that all interior points come first)
	for (n = 0; n < N; n++) {
		m_val = *(material_ir + n);
		if (m_val == *(neighbor_material_ir + n + 0*N) && m_val == *(neighbor_material_ir + n + 1*N) && m_val == *(neighbor_material_ir + n + 2*N) && m_val == *(neighbor_material_ir + n + 3*N) && m_val == *(neighbor_material_ir + n + 4*N) && m_val == *(neighbor_material_ir + n + 5*N)) { //If this point has the same composition as all its neighbors...
			*(interior_points + N_in) = *(target_points + n); //This point is an interior point
			*(material + N_in) = *(material_ir + n);
			*(neighbor_material + N_in + 0*N) = *(neighbor_material_ir + n + 0*N);
			*(neighbor_material + N_in + 1*N) = *(neighbor_material_ir + n + 1*N);
			*(neighbor_material + N_in + 2*N) = *(neighbor_material_ir + n + 2*N);
			*(neighbor_material + N_in + 3*N) = *(neighbor_material_ir + n + 3*N);
			*(neighbor_material + N_in + 4*N) = *(neighbor_material_ir + n + 4*N);
			*(neighbor_material + N_in + 5*N) = *(neighbor_material_ir + n + 5*N);
			*(Q + N_in) = *(Q_ir + n);
			N_in++;
		}
	}
	for (n = 0; n < N; n++) {
		m_val = *(material_ir + n);
		if (m_val != *(neighbor_material_ir + n + 0*N) || m_val != *(neighbor_material_ir + n + 1*N) || m_val != *(neighbor_material_ir + n + 2*N) || m_val != *(neighbor_material_ir + n + 3*N) || m_val != *(neighbor_material_ir + n + 4*N) || m_val != *(neighbor_material_ir + n + 5*N)) { //If this point's composition differs from at least one of its neighbors'...
			*(boundary_points + N_bound) = *(target_points + n); //This point is a boundary point
			*(material + N_in + N_bound) = *(material_ir + n);
			*(neighbor_material + N_in + N_bound + 0*N) = *(neighbor_material_ir + n + 0*N);
			*(neighbor_material + N_in + N_bound + 1*N) = *(neighbor_material_ir + n + 1*N);
			*(neighbor_material + N_in + N_bound + 2*N) = *(neighbor_material_ir + n + 2*N);
			*(neighbor_material + N_in + N_bound + 3*N) = *(neighbor_material_ir + n + 3*N);
			*(neighbor_material + N_in + N_bound + 4*N) = *(neighbor_material_ir + n + 4*N);
			*(neighbor_material + N_in + N_bound + 5*N) = *(neighbor_material_ir + n + 5*N);
			*(Q + N_in + N_bound) = *(Q_ir + n);
			N_bound++;
		}
	}
	
//	Constructing surface denominator matrix and volume "heat source" (see a description of the T-DDA method for an explanation of this part)
	float *D = (float *)calloc(N_bound*N_bound, sizeof(float)); //Denominator matrix
	float *Q_bound = (float *)calloc(N_bound, sizeof(float)); //"Heat source"
	float factor;
	int d_index, d_x, d_y, d_z, s_index, s_x, s_y, s_z;
	for (n = 0; n < N_bound; n++) {
		d_index = *(boundary_points + n);
		d_x = x_coords[d_index/((num_z + 1)*(num_y + 1))];
		d_y = y_coords[(d_index % ((num_z + 1)*(num_y + 1)))/(num_z + 1)];
		d_z = z_coords[d_index % (num_z + 1)];
		
		*(D + n + N_bound*n) = 1.0;
		for (p = 0; p < N_bound; p++) {
			s_index = *(boundary_points + p);
			s_x = x_coords[s_index/((num_z + 1)*(num_y + 1))];
			s_y = y_coords[(s_index % ((num_z + 1)*(num_y + 1)))/(num_z + 1)];
			s_z = z_coords[s_index % (num_z + 1)];
			

			if ( d_x - d < 0 ){

			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 0*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 0*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x + d - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 1*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 1*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x -d - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x - d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 2*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 2*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 3*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 3*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 4*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 4*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 5*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 5*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			}else if ( d_x < 0){

			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 0*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 0*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x + d - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 1*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 1*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x - d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 2*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 2*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 3*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 3*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 4*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 4*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 5*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 5*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));

			}else if( d_x + d < 0 ){

			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 0*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 0*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x + d - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 1*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 1*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x - d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 2*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 2*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 3*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 3*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 4*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 4*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 5*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 5*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));

			}else{

			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 0*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 0*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x + d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 1*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 1*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x - d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 2*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 2*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 3*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 3*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 4*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 4*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 5*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 5*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));

			}
		}
		
		*(Q_bound + n) = *(Q + N_in + n);
		for (l = 0; l < N_in; l++) {
			s_index = *(interior_points + l);
			s_x = x_coords[s_index/((num_z + 1)*(num_y + 1))];
			s_y = y_coords[(s_index % ((num_z + 1)*(num_y + 1)))/(num_z + 1)];
			s_z = z_coords[s_index % (num_z + 1)];
			
			if ( d_x - d < 0 ){

			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 0*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 0*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x + d - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + d + s_x) + x_lim*abs(d_y - s_y) 
			               + x_lim*y_lim*abs(d_z - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + 
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 1*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 1*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x -d - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x - d + s_x) + x_lim*abs(d_y - s_y) 
                                       + x_lim*y_lim*abs(d_z - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) +
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 2*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 2*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y + d - s_y) 
				       + x_lim*y_lim*abs(d_z - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + 
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 3*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 3*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - d - s_y) 
                                       + x_lim*y_lim*abs(d_z - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + 
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 4*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 4*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) 
                                       + x_lim*y_lim*abs(d_z + d - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) +
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 5*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 5*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) 
                                       + x_lim*y_lim*abs(d_z - d - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + 
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));

			}else if (d_x < 0){

			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 0*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 0*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x + d - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + d + s_x) + x_lim*abs(d_y - s_y) 
			               + x_lim*y_lim*abs(d_z - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + 
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 1*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 1*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x - d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
				       (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 2*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 2*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y + d - s_y) 
				       + x_lim*y_lim*abs(d_z - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + 
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 3*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 3*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - d - s_y) 
                                       + x_lim*y_lim*abs(d_z - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + 
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 4*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 4*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) 
                                       + x_lim*y_lim*abs(d_z + d - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) +
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 5*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 5*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) 
                                       + x_lim*y_lim*abs(d_z - d - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + 
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));

			
			}else if ( d_x + d < 0){

			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 0*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 0*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x + d - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + d + s_x) + x_lim*abs(d_y - s_y) 
                                       + x_lim*y_lim*abs(d_z - s_z)))) - (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 1*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 1*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x - d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
                                                                               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 2*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 2*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 3*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 3*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 4*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 4*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)))) - 
								               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 5*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 5*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)))) - 
								               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));

			}else{

			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 0*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 0*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x + d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 1*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 1*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x - d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
                                                                               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 2*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 2*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 3*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 3*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 4*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 4*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)))) - 
								               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 5*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 5*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)))) - 
								               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));

			}
		}
	}
	
//	Calculating "effective heat", i.e., bound charges
	int n_rhs = 1, info = 0;
	int *ipiv = (int *)malloc(N_bound*sizeof(int));
	float *Q_eff = (float *)malloc(N*sizeof(float));
	sgesv_(&N_bound, &n_rhs, D, &N_bound, ipiv, Q_bound, &N_bound, &info);
	free(D);
	for (n = 0; n < N; n++) {
		if (n < N_in) { //If the current point is an interior point
			*(Q_eff + n) = (*kappa)/(*(kappa + *(material + n)))*(*(Q + n)); //The heat source at this point is simply reduced by a factor of the target conductivity
		} else {
			*(Q_eff + n) = *(Q_bound + n - N_in); //The heat source is that computed above by solving the system of equations D*Q_eff = Q_bound
		}
	}
	
//	Writing header to output
//	fprintf(out_id, "%d %d %d %d %d %d %d\n", d, x_min, x_max, y_min, y_max, z_min, z_max);
//	fprintf(out_id, "%d %lu %lu %d\n", N, N_in, N_bound, N + N_out);
	
//	Calculating and writing out temperatures
	float *T = (float *)calloc(N + N_out, sizeof(float));
	for (n = 0; n < N + N_out; n++) { //Iterate through every point at which we want to determine the temperature
		if (n < N_in) {
			d_index = *(interior_points + n);
		} else if (n < N) {
			d_index = *(boundary_points + n - N_in);
		} else {
			d_index = *(outside_points + n - N);
		}
		d_x = x_coords[d_index/((num_z + 1)*(num_y + 1))];
		d_y = y_coords[(d_index % ((num_z + 1)*(num_y + 1)))/(num_z + 1)];
		d_z = z_coords[d_index % (num_z + 1)];
		
		for (l = 0; l < N; l++) { //Iterate through every point at which there is a heat source
			if (l < N_in) {
				s_index = *(interior_points + l);
			} else {
				s_index = *(boundary_points + l - N_in);
			}
			s_x = x_coords[s_index/((num_z + 1)*(num_y + 1))];
			s_y = y_coords[(s_index % ((num_z + 1)*(num_y + 1)))/(num_z + 1)];
			s_z = z_coords[s_index % (num_z + 1)];
			
			if ( d_x < 0) {
			  *(T + n) += 1/((*kappa)*unit*d)*( (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z))) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z))))*(*(Q_eff + l)); //This source point's contribution to the temperature at (d_x, d_y, d_z) is simply its effective heat multiplied by the Green's function value that takes (s_x, s_y, s_z) to (d_x, d_y, d_z)
			}else{			  
			  *(T + n) += 1/((*kappa)*unit*d)*( lower*(*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z))))*(*(Q_eff + l)); //This source point's contribution to the temperature at (d_x, d_y, d_z) is simply its effective heat multiplied by the Green's function value that takes (s_x, s_y, s_z) to (d_x, d_y, d_z)
			}



		}
		
		if (n < N) {
			fprintf(out_id, "%d %d %d %f %f %f\n", d_x, d_y, d_z, *(T + n), *(Q + n), *(Q_eff + n));
		} else {
			fprintf(out_id, "%d %d %d %f %f %f\n", d_x, d_y, d_z, *(T + n), 0.0, 0.0);
		}
	}

	fclose(dat_id);
	fclose(par_id);
	fclose(inp_id);	
	fclose(out_id);
	
	return 0;
}
