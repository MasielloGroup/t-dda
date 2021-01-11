//Determine whether 'target' is located in 'arr'
int find_neighbor(int *arr, int N, int target) {
	int imin = 0, imax = N - 1;
	int midpoint;
	while (imax >= imin) { //Check if 'arr' contains 'target'
		midpoint = (imax + imin)/2;
		if (*(arr + midpoint) == target) {
			return midpoint; //'target' is in 'arr'
		} else if (*(arr + midpoint) > target) {
			imax = midpoint - 1;
		} else {
			imin = midpoint + 1;
		}
	}
	
	return -1; //'target' is not in 'arr', so just return -1
}