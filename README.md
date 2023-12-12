# MACS2
MARSIS coherent surface clutter simulator

## SET-UP:
To use OpenMP set `#define OMP 1` in `setParallel.h`
To use MPI set `#define MPI 1` in `setParallel.h`

OpenMP and MPI can be used simultaneously.

To compile for big endian machines set `#define BIG_ENDIAN_COMP 1` in `marsislib.h` line 65

Select the complex variable error function algorithm by setting `MARSIS_erf_alg` in `marsislib.h` line 54

## COMPILATION:
Example using gcc:

### Serial
`gcc main.c marsislib.c faddeewa_w.c -lm -O3 -o macs`

### OpenMP:
`gcc main.c marsislib.c faddeewa_w.c -lm -O3 -fopenmp -o macs`

### MPI:
`mpicc main.c marsislib.c faddeewa_w.c -lm -O3 -o macs`

### MPI+OpenMP
`mpicc main.c marsislib.c faddeewa_w.c -lm -O3 -fopenmp -o macs`


## RUN:
`macs -m <path_to_MOLA_dir>  -i <path_to_intput_file> -o <output_file>`

Number of used nodes, number of OpenMP threads etc. defined by the scheduler job configuration file.

