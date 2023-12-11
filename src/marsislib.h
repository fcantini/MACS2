/* Physical and numerical constants */

#define C          299792458 // Velocity of light in vacuo in m/s
#define r_circle      100000 // Radius of the area contributing to scattering in m
#define r_sphere     3396000 // Mars mean radius in m
#define r_s_max      3417500 // Maximum radius of Mars
#define B              1.0e6 // MARSIS bandwidth in Hz
#define PRF           127.27 // MARSIS processing PRF in Hz
#define trx           350e-6 // MARSIS receiving window duration in s
#define FS             1.4e6 // MARSIS A/D converter sampling frequency in Hz
#define NS               512 // number of MARSIS echo samples
#define NBANDS			   2 // number of radar bands
#define RECORDBYTES		 160 //modified to comply the new sim format (9/4/2014)
#define EPS		  2.2204e-16
#define meno134	(-0.70710678118654746 + 0.70710678118654757 * I)

#define NFILTER			   3

#define COUNTS			   0
#define RADIUS			   1
#define TOPOGRAPHY		   2
#define DEF_SimDataFile "SIMDATA_SS3_TRK_CMP_EDR_2665.dat"	//Parameters input file
//#define SimDataFile "testMPIinput.dat"	//Parameters input file
//#define SimDataFile "input4punti.dat"	//Parameters input file
#define DEF_EdrSimFile  "SIM_SS3_TRK_CMP_EDR_2665.dat"		//Output file
#define DEF_MolaPath "../MOLA_data/"
#define comperffilterfile "comperffilter50.txt"


#define MARSIS_double	double
#define MARSIS_sqrt		sqrt
#define MARSIS_csqrt	csqrt
#define MARSIS_exp		exp
#define MARSIS_cexp		cexp
#define MARSIS_erf		erf
#define MARSIS_fabs		fabs
#define MARSIS_cos		cos
#define MARSIS_sin		sin
#define MARSIS_tan		tan
#define MARSIS_asin		asin
#define MARSIS_atan		atan
#define MARSIS_atan2	atan2
#define MARSIS_pow		pow
#define MARSIS_cpow		cpow
#define MARSIS_log		log
#define MARSIS_conj		conj
#define MARSIS_aton		atof
//#define MARSIS_cerf		erfz
#define	MARSIS_8		long
#define	MOLA_4			int
#define	MOLA_2			short


#define MARSIS_erf_alg	2 // 0-> Abramowitz-Stegun; 1-> Weideman; 2->Johnson 


#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288
#endif

#ifndef complex
#define complex _Complex
#endif

#define BIG_ENDIAN_COMP 0

//#define NPROC	2
#include "setParallel.h"
/******  STRUCTURES DEFINITION  *********************************/
typedef struct record_elem_struct{
	long int Nelem;
	long int Index;
}record_elem; 
/***************************************************************/

/***********************************************************************************/
/*          FUNCTIONS little_endian and big_endian                                 */
/*          Hardware endian detection                                              */
/*                                                                                 */
/*                           Last update 14/09/2012                                */
/***********************************************************************************/
int little_endian();
int big_endian();

/***********************************************************************************/
/*          FUNCTION readmarsissimrecords                                          */
/*          reads data length for a MARSIS orbit simulation                        */
/*                                                                                 */
/*                           Last update 13/09/2012                                */
/***********************************************************************************/
long int readmarsissimrecords(FILE *, int);


/***********************************************************************************/
/*          FUNCTION readmarsissimdata                                             */
/*          reads necessary data for a MARSIS orbit simulation                     */
/*                                                                                 */
/*                           Last update 09/04/2014                                */
/***********************************************************************************/
int readmarsissimdata(FILE *,
					  int,
					  long int,
					  double *,
					  double **,
					  double *,
					  double *,
					  double *,
					  double *,
					  double *,
					  double *,
					  double **,
					  double *,
					  double *,
					  double *,
					  double *,
					  double *,
					  double *,
					  double *,
					  double *,
					  double *);

/***********************************************************************************/
/*          FUNCTION sb_marsis_fread                                               */
/*          Swap bytes MARSIS readind file function                                */
/*                                                                                 */
/*                           Last update 30/11/2012                                */
/***********************************************************************************/
double sb_marsis_fread(FILE *fin);

/***********************************************************************************/
/*          FUNCTION sb_readmarsissimdata                                          */
/*          reads necessary data for a MARSIS orbit simulation                     */
/*                                                                                 */
/*                           Last update 09/04/2014                                */
/***********************************************************************************/
int sb_readmarsissimdata(FILE *fin,
						 int Nbands,
						 long int FileRecords,
						 double *ostline,
						 double **f0,
						 double *theta,
						 double *frameid,
						 double *scetfw,
						 double *scetff,
						 double *Vt,
						 double *Vr,
						 double **NA,
						 double *x0,
						 double *y0,
						 double *z0,
						 double *alt0,
						 double *Vx0,
						 double *Vy0,
						 double *Vz0,
						 double *lon_0,
						 double *lat_0);

/***********************************************************************************/
/*          FUNCTION sb_marsis_fwrite                                              */
/*          Swap bytes fwrite                                                      */
/*                                                                                 */
/*                           Last update 30/11/2012                                */
/***********************************************************************************/
int sb_marsis_fwrite(double *OutBuf, int ns, FILE *fout);

/***********************************************************************************/
/*          FUNCTION sb_writemarsissimdata                                         */
/*          Swap bytes Write simulation parameters in the output file              */
/*                                                                                 */
/*                           Last update 09/04/2014                                */
/***********************************************************************************/
int sb_writemarsissimdata(FILE *fout,
						  int Nbands,
						  long int index,					   
						  double *ostline,
						  double **f0,
						  double *theta,
						  double *frameid,
						  double *scetfw,
						  double *scetff,
						  double *Vt,
						  double *Vr,
						  double **NA,
						  double *x0,
						  double *y0,
						  double *z0,
						  double *alt0,
						  double *Vx0,
						  double *Vy0,
						  double *Vz0,
						  double *lon_0,
						  double *lat_0);

/***********************************************************************************/
/*          FUNCTION fftvars                                                       */
/*          Computes vectors containing the correct time and frequency values for  */
/*          FFT computation                                                        */
/*                                                                                 */
/*                           Last update 26/09/2012                                */
/***********************************************************************************/

int fftvars(MARSIS_double fs, long int Ns, MARSIS_double *t, MARSIS_double *f);

#if MPI == 1
/***********************************************************************************/
/*          FUNCTION writemarsissimdata_MPI                                        */
/*          Write simulation parameters in the output file                         */
/*                                                                                 */
/*                           Last update 09/04/2014                                */
/***********************************************************************************/
int writemarsissimdata_MPI(MPI_File MPI_fout,
						   MPI_Offset offset,
						   MPI_Request *request,
						   int Nbands,
						   long int index,					   
						   double *ostline,
						   double **f0,
						   double *theta,
						   double *frameid,
						   double *scetfw,
						   double *scetff,
						   double *Vt,
						   double *Vr,
						   double **NA,
						   double *x0,
						   double *y0,
						   double *z0,
						   double *alt0,
						   double *Vx0,
						   double *Vy0,
						   double *Vz0,
						   double *lon_0,
						   double *lat_0);

/***********************************************************************************/
/*          FUNCTION sb_marsis_MPI_File_iwrite                                     */
/*          Swap bytes MPI_File_iwrite                                             */
/*                                                                                 */
/*                           Last update 30/11/2012                                */
/***********************************************************************************/
int sb_marsis_MPI_File_iwrite(MPI_File MPI_fout, double *OutBuf, int ns, MPI_Request *request);

/***********************************************************************************/
/*          FUNCTION sb_writemarsissimdata_MPI                                     */
/*          Write simulation parameters in the output file                         */
/*                                                                                 */
/*                           Last update 09/04/2014                                */
/***********************************************************************************/
int sb_writemarsissimdata_MPI(MPI_File MPI_fout,
							  MPI_Offset offset,
							  MPI_Request *request,
							  int Nbands,
							  long int index,					   
							  double *ostline,
							  double **f0,
							  double *theta,
							  double *frameid,
							  double *scetfw,
							  double *scetff,
							  double *Vt,
							  double *Vr,
							  double **NA,
							  double *x0,
							  double *y0,
							  double *z0,
							  double *alt0,
							  double *Vx0,
							  double *Vy0,
							  double *Vz0,
							  double *lon_0,
							  double *lat_0);

#endif
/***********************************************************************************/
/*          FUNCTION writemarsissimdata                                            */
/*          Write simulation parameters in the output file                         */
/*                                                                                 */
/*                           Last update 09/04/2014                                */
/***********************************************************************************/
int writemarsissimdata(FILE *,
					   int ,
					   long int,
					   double *,
					   double **,
					   double *,
					   double *,
					   double *,
					   double *,
					   double *,
					   double *,
					   double **,
					   double *,
					   double *,
					   double *,
					   double *,
					   double *,
					   double *,
					   double *,
					   double *,
					   double *);

/***********************************************************************************/
/*          FUNCTION molagridsize                                                  */
/*          MOLA set size calculation                                              */
/*                                                                                 */
/*                           Last update 13/10/2012                                */
/***********************************************************************************/
int molagridsize(double lon_0, 
				 double lat_0, 
				 double R_circle, 
				 double R_sphere,  
				 MARSIS_double *lon_w, 
				 MARSIS_double *lon_e, 
				 MARSIS_double *lat_s, 
				 MARSIS_double *lat_n, 
				 long int *N_I, 
				 long int *N_J);

/***********************************************************************************/
/*          FUNCTION bestmolagrid                                                  */
/*          Best MOLA set selection                                                */
/*                                                                                 */
/*                           Last update 13/10/2012                                */
/***********************************************************************************/
int bestmolagrid(double lon_0, 
				 double lat_0, 
				 double R_circle, 
				 double R_sphere, 
				 MARSIS_double lon_w, 
				 MARSIS_double lon_e, 
				 MARSIS_double lat_s, 
				 MARSIS_double lat_n, 
				 MARSIS_double **mola, 
				 MARSIS_double **lon, 
				 MARSIS_double **lat,
				 long int N_i,
				 long int N_j);

/***********************************************************************************/
/*          FUNCTION lonlatbox                                                     */
/*          Used by bestmolagrid. Finds longitude and latitude                     */
/*          box enclosing a circle of a given radiuson a sphere                    */
/*                                                                                 */
/*                           Last update 26/09/2012                                */
/***********************************************************************************/
int lonlatbox(double lon_0, 
			  double lat_0, 
			  double R_circle, 
			  double R_sphere, 
			  MARSIS_double *lon_w, 
			  MARSIS_double *lon_e, 
			  MARSIS_double *lat_s, 
			  MARSIS_double *lat_n);

/***********************************************************************************/
/*         FUNCTION molagrid                                                       */
/*         Used by bestmolagrid. MOLA MEGDR values corresponding to a              */
/*         regularly-spaced grid of % planetocentric coordinates.                  */
/*                                                                                 */
/*                           Last update 16/11/2012                                */
/***********************************************************************************/
int molagrid(MARSIS_double lon_w, 
			 MARSIS_double lon_e, 
			 MARSIS_double lat_s, 
			 MARSIS_double lat_n, 
			 int type, 
			 double ppd, 
			 MARSIS_double **mola, 
			 MARSIS_double **lon, 
			 MARSIS_double **lat,
			 long int N_i,
			 long int N_j);

/***********************************************************************************/
/*         FUNCTION molapolargrid                                                  */
/*         Used by bestmolagrid.                                                   */
/*         extract a subset of the MOLA POLAR dataset                              */
/*                                                                                 */
/*                           Last update 16/11/2012                                */
/***********************************************************************************/
int molapolargrid(double lon_0, 
				  double lat_0, 
				  double delta_x,
				  double delta_y,
				  double type,
				  double ppd,
				  MARSIS_double **mola, 
				  MARSIS_double **lon, 
				  MARSIS_double **lat,
				  long int N_i,
				  long int N_j);

/***********************************************************************************/
/*         FUNCTION surfaspect                                                     */
/*         given a surface described by a set of points expressed in cartesian     */
/*         coordinates and a point external to the surface, computes the size,     */
/*         orientation and distance of each portion of the discretized surface     */
/*         w.r.t. the external point                                               */
/*                                                                                 */
/*                           Last update 26/09/2012                                */
/***********************************************************************************/
int surfaspect(MARSIS_double **X, 
			   MARSIS_double **Y, 
			   MARSIS_double **Z, 
			   double x0, 
			   double y0, 
			   double z0,
			   long int N_i,
			   long int N_j,
			   MARSIS_double **la, 
			   MARSIS_double **lb, 
			   MARSIS_double **Ux, 
			   MARSIS_double **Uy, 
			   MARSIS_double **Uz, 
			   MARSIS_double **R);

/***********************************************************************************/
/*         FUNCTION facet2nd                                                       */
/*         Determine the complex reflectivity of all facets                        */
/*                                                                                 */
/*                           Last update 26/09/2012                                */
/***********************************************************************************/
int facet2nd(MARSIS_double **la, 
			 MARSIS_double **lb, 
			 MARSIS_double **Ux, 
			 MARSIS_double **Uy, 
			 MARSIS_double **Uz, 
			 MARSIS_double **R, 
			 MARSIS_double f, 
			 MARSIS_double complex **E, 
			 long int N_i, 
			 long int N_j,
			 MARSIS_double L, 
			 MARSIS_double *coef, 
			 int ncoef);

/***********************************************************************************/
/*         FUNCTION erfz                                                           */
/*         Error function of complex numbers                                       */
/*                                                                                 */
/*                           Last update 26/09/2012                                */
/***********************************************************************************/
MARSIS_double complex erfz(MARSIS_double complex z);

/***********************************************************************************/
/*         FUNCTION parts                                                          */
/*         Used by erfz: Evaluates all partial functions E(z) to G(z).             */
/*                                                                                 */
/*                           Last update 26/09/2012                                */
/***********************************************************************************/
MARSIS_double complex parts(MARSIS_double Re, MARSIS_double Im);

/***********************************************************************************/
/*         FUNCTION marsisazproc                                                   */
/*         Reproduces the effect of synthetic aperture processing on signals       */
/*         coming from different directions                                        */
/*                                                                                 */
/*                           Last update 26/09/2012                                */
/***********************************************************************************/
int marsisazproc(MARSIS_double f, 
				 double Vr, 
				 double Vt, 
				 double theta_s, 
				 int m, 
				 double prf, 
				 MARSIS_double complex NA, 
				 MARSIS_double **X, 
				 MARSIS_double **Y,
				 MARSIS_double **Z, 
				 double x0, 
				 double y0, 
				 double z0, 
				 double Vx0, 
				 double Vy0, 
				 double Vz0,
				 MARSIS_double complex **G,
				 long int N_i,
				 long int N_j);

/***********************************************************************************/
/*         FUNCTION readcomperffilterlen                                           */
/*         Returns the number of lines in complex erf filter text file.            */
/*         Returns -1 if an error occours                                          */
/*                                                                                 */
/*                           Last update 05/10/2012                                */
/***********************************************************************************/
int readcomperffilterlen(FILE *fid);


/***********************************************************************************/
/*         FUNCTION readcomperffilter                                              */
/*         Returns complex erf filter coefficients.                                */
/*                                                                                 */
/*                           Last update 05/10/2012                                */
/***********************************************************************************/
int readcomperffilter(FILE *fid, MARSIS_double *L, MARSIS_double *coef, int ncoef);

/***********************************************************************************/
/*         FUNCTION cerf                                                           */
/*         Returns erf value with complex argument.                                */
/*                                                                                 */
/*                           Last update 05/10/2012                                */
/***********************************************************************************/
MARSIS_double complex cerf(MARSIS_double complex z, MARSIS_double L, MARSIS_double *coef, int ncoef);

/***********************************************************************************/
/*         FUNCTION struct_cmp_by_elements                                         */
/*         Comparing function for qsort.                                           */
/*                                                                                 */
/*                           Last update 14/10/2012                                */
/***********************************************************************************/
int struct_cmp_by_elements(const void *a, const void *b);

/***********************************************************************************/
/*         FUNCTION select_proc_elem                                               */
/*         Selects trajectory points to be simulated by each MPI process.          */
/*                                                                                 */
/*                           Last update 14/10/2012                                */
/***********************************************************************************/
int select_MPI_proc_elem_NO(record_elem *Relem, long int Relem_len, long int *ProcElem_MPI, long int *N_ProcElem_MPI,  int N_proc_MPI, int Proc_MPI);

/***********************************************************************************/
/*         FUNCTION select_proc_elem                                             */
/*         Selects trajectory points to be simulated by each MPI process.          */
/*                                                                                 */
/*                           Last update 16/10/2012                                */
/***********************************************************************************/
int select_MPI_proc_elem(record_elem *Relem, long int Relem_len, long int *ProcElem_MPI, long int *N_ProcElem_MPI,  int N_proc_MPI, int Proc_MPI);

/***********************************************************************************/
/*         FUNCTION Faddeeva_w                                                     */
/*         From http://ab-initio.mit.edu/wiki/index.php/Faddeeva_w                 */
/*                                                                                 */
/*         Slightly modified form C++ to C                                         */
/***********************************************************************************/
double complex Faddeeva_w(double complex z, double relerr);

double complex Faddeeva_erf(double complex z, double relerr);
