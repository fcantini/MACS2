/***********************************************************************************/
/*          PROGRAM MARSIS_coherent_simulator                                      */
/*                                                                                 */
/*                           Last update 09/04/2014                                */
/***********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <string.h>
#include "setParallel.h"
#if MPI == 1
#include <mpi.h>
#endif
#if OMP == 1
#include <omp.h>
#endif
#include "marsislib.h"

#define LOOP_COMP_OFF 0
#define MULTI_OUT_FILE 0

char *MolaPath;


int main (int argc, char *argv[]) {
	
	/*******************  VARIABLES DECLARATION   *******************************************************************/	
	char *SimDataFile;//[50] = "SIMDATA_SS3_TRK_CMP_EDR_2665.dat";
	char *EdrSimFile;//[50]  = "SIM_SS3_TRK_CMP_EDR_2665.dat";
	
	FILE *fin;
#if MARSIS_erf_alg == 1
	FILE *ffilt;
#endif
#if MPI == 1
	MPI_File MPI_fout;
	MPI_Request request; 
	MPI_Status status;
#else
	FILE *fout;
#endif
#if MULTI_OUT_FILE == 1
	FILE *PointFile;
	char PointFilename[100];
	char Tstring[20];
#endif
	double *ostline;
	double *f0[NBANDS];
	double *theta;
	double *frameid;
	double *scetfw;
	double *scetff;
	double *Vt;
	double *Vr;
	double *NA[NBANDS];
	double *x0;
	double *y0;
	double *z0;
	double *alt0;
	double *Vx0;
	double *Vy0;
	double *Vz0;
	double *lon_0;
	double *lat_0;
	long int ii;
	long int jj;
	long int kk;
	long int ll;
	long int NP;
	int mm, nn;
	long int FileRecords;
	MARSIS_double *t;
	MARSIS_double *f, *fband;
	int	   *iband;
	int		Nf;
	MARSIS_double **Z_mola;
	MARSIS_double **X_lon;
	MARSIS_double **Y_lat;
	MARSIS_double **la, **lb, **Ux, **Uy, **Uz, **R;
	MARSIS_double complex **E, **G;
	MARSIS_double *lon_w;
	MARSIS_double *lon_e;
	MARSIS_double *lat_s;
	MARSIS_double *lat_n;
	long int *N_I, *N_J;
	long int N_I_max, N_J_max;	
	double rcoselev,az,elev,deltat;
	MARSIS_double complex phase[NS];
//	MARSIS_double complex espectr[NS][NFILTER];
	MARSIS_double complex espectr;
//	MARSIS_double Respectr[NFILTER][NS];
//	MARSIS_double Iespectr[NFILTER][NS];
	double Respectr[NFILTER][NS];
	double Iespectr[NFILTER][NS]; // perché deve scrivere double	
	MARSIS_double *comperffilt;
	int comperffiltlen;
	MARSIS_double L;
	record_elem *Relem;
	long int *ProcElem_MPI;
	long int N_ProcElem_MPI;
//	double empty;
#if MPI == 1
	int RecordSize, RecHeadSize, RecDataSize;
#endif
	
	/* MPI VARIABLES */
	int commSize, commRank;
	/*****************/
#if OMP == 1
	double start, end;
#else
	clock_t start, end;
#endif
	double elapsed;
	
	/****************************************************************************************************************/

#if MPI == 1
	// Initialize MPI state
    MPI_Init(&argc, &argv);
	
    // Get our MPI node number and node count
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);	
#else
	commSize = 1;
	commRank = 0;
#endif	
	
#if OMP == 1
	start = omp_get_wtime();
#else
	start = clock();
#endif
	
	/***********************************/
//    if(big_endian())
/*        printf("big_endian\n");
    if(little_endian())
        printf("little_endian\n");*/
	/***********************************/
//	MPI_Finalize();
//	return 0;
	/****ENDIANNES CHECK ***************/
	if (big_endian() != BIG_ENDIAN_COMP)
	{
		printf("ENDIANNES consistency check: wrong endianness selection at compiling time.");
		return -1;				
	}
	/***********************************/	
	/****MOLA DATA SIZE CHECK **********/
	if (sizeof(MOLA_2)!=2){
		printf("MOLA data size consistency check:Reading size for function molagrid is %ld. Required 2. \n", sizeof(MOLA_2));
		return -1;		
	}
	if (sizeof(MOLA_4)!=4){
		printf("MOLA data size consistency check:Reading size for function molapolargrid is %ld. Required 4. \n", sizeof(MOLA_4));
		return -1;		
	}	
	/***********************************/
	/******  DEFAULT  INPUT MANAGEMENT *************/
	SimDataFile  = malloc((strlen(DEF_SimDataFile)+1) * sizeof(char));
	strcpy(SimDataFile,DEF_SimDataFile);

	EdrSimFile  = malloc((strlen(DEF_EdrSimFile)+1) * sizeof(char));
	strcpy(EdrSimFile,DEF_EdrSimFile);

	MolaPath  = malloc((strlen(DEF_MolaPath)+1) * sizeof(char));
	strcpy(MolaPath,DEF_MolaPath);
	
	/******  OPTIONAL INPUT MANAGEMENT *************/
	for (ii = 1; ii < argc; ii+=2)
	{
		if (strcmp(argv[ii],"-i") == 0)
		{
			free(SimDataFile);
			SimDataFile  = malloc((strlen(argv[ii+1])+1) * sizeof(char));
			strcpy(SimDataFile,argv[ii+1]);
		}
		if (strcmp(argv[ii],"-o") == 0)
		{
			free(EdrSimFile);
			EdrSimFile  = malloc((strlen(argv[ii+1])+1) * sizeof(char));
			strcpy(EdrSimFile,argv[ii+1]);
		}
		if (strcmp(argv[ii],"-m") == 0)
		{
			free(MolaPath);
			MolaPath  = malloc((strlen(argv[ii+1])+1) * sizeof(char));
			strcpy(MolaPath,argv[ii+1]);
		}		
	}
	/***********************************************/
	if (commRank == 0)
	{
		printf("Orbit data INPUT  file: %s\n", SimDataFile);
		printf("Simulation OUTPUT file: %s\n", EdrSimFile);
		printf("MOLA data path: %s\n", MolaPath);
	}
//	return 0;
	/******************************   READING Weideman complex erf filter coefficients  ********************************/
#if MARSIS_erf_alg	== 1
	ffilt = fopen(comperffilterfile, "r");
	if (ffilt == NULL){
		printf("ReadMarsisSimData:MissingFilterFile - Unable to open %s\n", comperffilterfile);
		return -1;
	}
	comperffiltlen = readcomperffilterlen(ffilt);
	comperffiltlen--;
	comperffilt = malloc((comperffiltlen)*sizeof(MARSIS_double));
	readcomperffilter(ffilt, &L, comperffilt, comperffiltlen);
	fclose(ffilt);
#else
	comperffilt = malloc(sizeof(MARSIS_double));
	comperffilt[0] = 0;
#endif
//	for(ii=0;ii<comperffiltlen;ii++)
//		printf("coef[%d]=%g\n",ii,comperffilt[ii]);
//	printf("L=%g\n",L);
//	return 0;
	/******************************************************************************************************************/	
	/******************************   READING instrument and geometric parameters  *********************************/
	/* Opening input file */
	fin = fopen(SimDataFile, "rb");
	if (fin == NULL){
		printf("ReadMarsisSimData:MissingInputFile - The required simulation data file %s could not be opened\n", SimDataFile);
		return -1;
	}	
//	fcheck = fopen("fcheckC.txt","wb");
	/* Retrieving the number of file records*/
	FileRecords = readmarsissimrecords(fin, RECORDBYTES);
	
	/* Memory allocation */
	ostline = malloc(		FileRecords * sizeof(double));
	for (jj=0; jj<NBANDS; jj++)
		f0[jj] = malloc(	FileRecords * sizeof(double));
	theta	= malloc(		FileRecords * sizeof(double));
	frameid = malloc(		FileRecords * sizeof(double));
	scetfw	= malloc(		FileRecords * sizeof(double));
	scetff	= malloc(		FileRecords * sizeof(double));
	Vt		= malloc(		FileRecords * sizeof(double));
	Vr		= malloc(		FileRecords * sizeof(double));
	for (jj=0; jj<NBANDS; jj++)
		NA[jj] = malloc(	FileRecords * sizeof(double));	
	x0		= malloc(		FileRecords * sizeof(double));
	y0		= malloc(		FileRecords * sizeof(double));
	z0		= malloc(		FileRecords * sizeof(double));
	alt0		= malloc(		FileRecords * sizeof(double));
	Vx0		= malloc(		FileRecords * sizeof(double));
	Vy0		= malloc(		FileRecords * sizeof(double));
	Vz0		= malloc(		FileRecords * sizeof(double));
	lon_0	= malloc(		FileRecords * sizeof(double));
	lat_0	= malloc(		FileRecords * sizeof(double));
	
	lon_w   = malloc(		FileRecords * sizeof(double)); 
	lon_e   = malloc(		FileRecords * sizeof(double));
	lat_s   = malloc(		FileRecords * sizeof(double));
	lat_n   = malloc(		FileRecords * sizeof(double));
	N_I     = malloc(		FileRecords * sizeof(long int));
	N_J     = malloc(		FileRecords * sizeof(long int));
	
	Relem		= malloc(		FileRecords * sizeof(record_elem));
	ProcElem_MPI= malloc(		FileRecords * sizeof(long int));

	/*Data reading*/
#if	BIG_ENDIAN_COMP	== 1
	sb_readmarsissimdata(fin, NBANDS, FileRecords, ostline, f0, theta, frameid, scetfw, scetff, Vt, Vr, NA, x0, y0, z0, alt0, Vx0, Vy0, Vz0, lon_0, lat_0);
//	printf("sb_readmarsissimdata \n");
#else
	readmarsissimdata(fin, NBANDS, FileRecords, ostline, f0, theta, frameid, scetfw, scetff, Vt, Vr, NA, x0, y0, z0, alt0, Vx0, Vy0, Vz0, lon_0, lat_0);
#endif
	/*Closing input file*/
	fclose(fin);
	/***************************************************************************************************************/

	/*****************************************************************/
	/* Preparing quantities used in radar computations**************************************************************/
	t = malloc(NS*sizeof(MARSIS_double));
	f = malloc(NS*sizeof(MARSIS_double));
	fftvars(FS, NS, t, f);


	iband = malloc(NS*sizeof(int));
	Nf=0;
	for (ii=0; ii<NS; ii++)
		if (fabs(f[ii]) <= B/2)
			iband[Nf++] = ii;
		
	fband = malloc(Nf*sizeof(MARSIS_double));
	for (ii=0; ii<Nf; ii++)
		fband[ii]=f[iband[ii]];
	/***************************************************************************************************************/
	/************** Extracting MOLA grid size ****************/
	N_ProcElem_MPI = 0;
	for (ii=0; ii<FileRecords; ii++)		
	{
		molagridsize(lon_0[ii], 
					 lat_0[ii], 
					 r_circle, 
					 r_sphere, 
					 &lon_w[ii], 
					 &lon_e[ii], 
					 &lat_s[ii], 
					 &lat_n[ii], 
					 &N_I[ii], 
					 &N_J[ii]);

		Relem[ii].Index = ii;
		Relem[ii].Nelem = N_I[ii]*N_J[ii];
	/*********************************************************/	
		ProcElem_MPI[ii] = 0;
	}
	/************** Selecting orbital points for current MPI process ****************/
#if MPI == 1	
	qsort(Relem, FileRecords, sizeof(struct record_elem_struct), struct_cmp_by_elements);
	select_MPI_proc_elem(Relem, FileRecords, ProcElem_MPI, &N_ProcElem_MPI, commSize, commRank);
#else
	for (ii=0; ii<FileRecords; ii++)
		ProcElem_MPI[ii] = ii;
	N_ProcElem_MPI = FileRecords;
	printf("Number of points %ld \n", FileRecords);
#endif
	/*******************************************************************************/	
	
	/***** Calculating current process max N_I, N_J  for matrix allocation *********/	
	N_I_max = 0;
	N_J_max = 0;
	
	for (NP=0; NP<N_ProcElem_MPI; NP++)
	{
		ii = ProcElem_MPI[NP]; /* ii = ith point of the orbit */
		if (N_I[ii] > N_I_max)
			N_I_max = N_I[ii];
		 
		if (N_J[ii] > N_J_max)
			N_J_max = N_J[ii];
	}
	/*******************************************************************************/		

	/*********************************************************/
	/* MOLA Memory allocation **********************/
	Z_mola = malloc(N_I_max * sizeof(MARSIS_double *));
	X_lon  = malloc(N_I_max * sizeof(MARSIS_double *));
	Y_lat  = malloc(N_I_max * sizeof(MARSIS_double *));
	for(jj = 0; jj < N_I_max; jj++)
	{
		Z_mola[jj] = malloc(N_J_max * sizeof(MARSIS_double));
		X_lon[jj]  = malloc(N_J_max * sizeof(MARSIS_double));
		Y_lat[jj]  = malloc(N_J_max * sizeof(MARSIS_double));
	}
	/***********************************************/
	/* SURFASPECT Memory allocation *****************/
	la = malloc(N_I_max * sizeof(double *));
	lb = malloc(N_I_max * sizeof(double *));
	Ux = malloc(N_I_max * sizeof(double *));
	Uy = malloc(N_I_max * sizeof(double *));
	Uz = malloc(N_I_max * sizeof(double *));
	R  = malloc(N_I_max * sizeof(double *));
	for(jj = 0; jj < N_I_max; jj++)
	{
		la[jj] = malloc(N_J_max * sizeof(double));
		lb[jj] = malloc(N_J_max * sizeof(double));
		Ux[jj] = malloc(N_J_max * sizeof(double));
		Uy[jj] = malloc(N_J_max * sizeof(double));
		Uz[jj] = malloc(N_J_max * sizeof(double));
		R[jj]  = malloc(N_J_max * sizeof(double));
	}
	/***********************************************/
	/* FACET2ND MARSISAZPROC Memory allocation *****	
	E = malloc(N_I_max * sizeof(double complex *));
	G = malloc(N_I_max * sizeof(double complex *));
	for(jj = 0; jj < N_I_max; jj++)
	{
		E[jj] = malloc(N_J_max * sizeof(double complex));
		G[jj] = malloc(N_J_max * sizeof(double complex));
	}*/	
	/***********************************************/	
#if MPI == 1
//	MPI_fout = NULL; /* To avoid warning with intel cmpiler */
	MPI_File_open(MPI_COMM_WORLD, EdrSimFile, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &MPI_fout);
//	MPI_File_set_view(MPI_fout, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
	MPI_File_set_atomicity(MPI_fout, 1);
	
	RecDataSize =sizeof(double)*(2*NS*NFILTER);
//	RecHeadSize= sizeof(double)*(NBANDS*2+15);
	RecHeadSize= RECORDBYTES;
	RecordSize = NBANDS*RecDataSize+RecHeadSize;
	MPI_File_preallocate(MPI_fout, FileRecords*RecordSize);
#else
	fout = fopen(EdrSimFile, "wb");
#endif

	/***************************************************************************************************************/
	/*  STARTING MAIN LOOP: for each echo to be simulated                                                          */
	/***************************************************************************************************************/
//	for (ii=0; ii<FileRecords; ii++)
//	for (NP=0; NP<1; NP++)	
	for (NP=0; NP<N_ProcElem_MPI; NP++)
	{
		ii = ProcElem_MPI[NP]; /* ii = ith point of the orbit */
		/*Data writing*/
#if MPI == 1
	#if BIG_ENDIAN_COMP == 1
		sb_writemarsissimdata_MPI(MPI_fout, ii*RecordSize, &request, NBANDS, ii, ostline, f0, theta, frameid, scetfw, scetff, Vt, Vr, NA, x0, y0, z0, alt0, Vx0, Vy0, Vz0, lon_0, lat_0);	
//		printf("sb_writemarsissimdata_MPI \n");
	#else
		writemarsissimdata_MPI(MPI_fout, ii*RecordSize, &request, NBANDS, ii, ostline, f0, theta, frameid, scetfw, scetff, Vt, Vr, NA, x0, y0, z0, alt0, Vx0, Vy0, Vz0, lon_0, lat_0);	
	#endif		
#else
	#if BIG_ENDIAN_COMP == 1
		sb_writemarsissimdata(fout, NBANDS, ii, ostline, f0, theta, frameid, scetfw, scetff, Vt, Vr, NA, x0, y0, z0, alt0, Vx0, Vy0, Vz0, lon_0, lat_0);	
//		printf("sb_writemarsissimdata \n");
	#else	
		printf("Processing point %ld ", NP);
		writemarsissimdata(fout, NBANDS, ii, ostline, f0, theta, frameid, scetfw, scetff, Vt, Vr, NA, x0, y0, z0, alt0, Vx0, Vy0, Vz0, lon_0, lat_0);	
	#endif
#endif	
#if MULTI_OUT_FILE == 1
		strcpy(PointFilename,"");

		sprintf(Tstring, "%ld", ii);
		strcat(PointFilename,Tstring);

		strcat(PointFilename, "_");
		strcat(PointFilename, EdrSimFile);
		
		PointFile = fopen(PointFilename, "wb");
	#if BIG_ENDIAN_COMP == 1
		sb_writemarsissimdata(PointFile, NBANDS, 0, ostline, f0, theta, frameid, scetfw, scetff, Vt, Vr, NA, x0, y0, z0, alt0, Vx0, Vy0, Vz0, lon_0, lat_0);	
	#else		
		writemarsissimdata(PointFile, NBANDS, 0, ostline, f0, theta, frameid, scetfw, scetff, Vt, Vr, NA, x0, y0, z0, alt0, Vx0, Vy0, Vz0, lon_0, lat_0);	
	#endif
#endif
//		printf("lon=%f lat=%f\n", lon_0[ii], lat_0[ii]);
//		fclose(fout);
//		return 0;
#if LOOP_COMP_OFF == 0		
		/************** Extracting MOLA topography ***************/
//		molagridsize(lon_0[ii], lat_0[ii], r_circle, r_sphere, &lon_w, &lon_e, &lat_s, &lat_n, &N_I, &N_J);
//		fprintf (fcheck,"%ld %ld %ld\n",ii+1, N_I, N_J);
		bestmolagrid(lon_0[ii],lat_0[ii], r_circle, r_sphere, lon_w[ii], lon_e[ii], lat_s[ii], lat_n[ii], Z_mola, X_lon, Y_lat, N_I[ii], N_J[ii]);
		
//		printf("ii=%d",ii);
		/************** END Extracting MOLA topography ***********/
/*						int yy,zz;
		 for(yy=0; yy<N_I[ii]; yy++)
		 {
		 for(zz=0; zz<N_J[ii]; zz++)
		 {
		 printf("%8.20f\t",Z_mola[yy][zz]);
		 }
		 printf("\n");
		 }
		 return 0;*/
		
		// DEBUG
#endif
		/************** Inizializing  espectr ***************/
		for (jj=0; jj<NFILTER; jj++)
			for (kk=0; kk<NS; kk++)
				Respectr[jj][kk] = Iespectr[jj][kk] = 0;
		/****************************************************/		
				
#if LOOP_COMP_OFF == 0
		/******** Compute position and orientation of the observed surface ******/
		/* spherical to cartesian coordinate */

		for(kk=0; kk<N_I[ii]; kk++)
		{
			for(ll=0; ll<N_J[ii]; ll++)
			{
				az = M_PI/180.0*X_lon[kk][ll];
				elev   = M_PI/180.0*Y_lat[kk][ll];
				rcoselev = Z_mola[kk][ll] * MARSIS_cos(elev);				

				X_lon[kk][ll] = rcoselev * MARSIS_cos(az);
				Y_lat[kk][ll] = rcoselev * MARSIS_sin(az);
				Z_mola[kk][ll] = Z_mola[kk][ll] * MARSIS_sin(elev);
			}
		}
		/*************************************/		
		surfaspect(X_lon, Y_lat, Z_mola, x0[ii], y0[ii], z0[ii], N_I[ii], N_J[ii], la, lb, Ux, Uy, Uz, R);
		// DEBUG
	
		/*  Compute two-way travel time difference between the current spacecraft
		 position and a reference distance from Mars equal to the highest value
		 of the Martian radius, to be used in aligning echoes to a common time
		 reference		*/
		
		deltat = 2 * (MARSIS_sqrt( x0[ii]*x0[ii] + y0[ii]*y0[ii] + z0[ii]*z0[ii] ) - r_s_max ) / C;
#endif	
#if OMP == 1	
#pragma omp parallel private (E, G, jj, kk, phase)
{
//			printf("omp threads = %d\n", omp_get_num_threads());
#endif
			/* FACET2ND MARSISAZPROC Memory allocation *****/	
			E = malloc(N_I_max * sizeof(double complex *));
			G = malloc(N_I_max * sizeof(double complex *));
			for(jj = 0; jj < N_I_max; jj++)
			{
				E[jj] = malloc(N_J_max * sizeof(double complex));
				G[jj] = malloc(N_J_max * sizeof(double complex));
			}	
			/***********************************************/	
			
		/***************************************************************************************************************/
		/*  STARTING BANDS LOOP: for each band to be simulated                                                         */
		/***************************************************************************************************************/
		for (jj=0; jj<NBANDS; jj++)
		{
#if LOOP_COMP_OFF == 0			
			for(kk=0; kk<NS; kk++)
				phase[kk] =  MARSIS_cexp( 2 * I * M_PI * deltat * (f[kk] + f0[jj][ii]));
#endif			
//			printf("cabs: %g done\n",cabs(1+2*I));
//			printf("test: %g\n", erf(-INFINITY));
			/***************************************************************************************************************/
			/*  STARTING IN BANDS FREQUENCIES LOOP: for each frequency in the band to be simulated                         */
			/***************************************************************************************************************/			
#if OMP == 1
#pragma omp for private (ll, mm, nn, espectr)			
#endif
			for (kk=0; kk<Nf; kk++)
			{
//				printf("Band=%ld/%d Freq=%ld/%d \n",jj+1,NBANDS,kk+1,Nf);
//				start = clock();
#if LOOP_COMP_OFF == 0				
				facet2nd( la, lb, Ux, Uy, Uz, R, fband[kk] + f0[jj][ii], E, N_I[ii], N_J[ii], L, comperffilt, comperffiltlen);		
#endif				
//				end = clock();
//				elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
//				printf("Elapsed time (facet2nd): %g \n", elapsed);

/*				FILE *fmt;
				fmt = fopen("MatrixE_C_abram.dat","wt");
				for(jj = 0; jj < N_I[ii]; jj++)
				{ 
					for(kk = 0; kk < N_J[ii]; kk++)	
					{
						fprintf(fmt,"row=%d col=%d %lg %lg",jj+1,kk+1,creal(E[jj][kk]),cimag(E[jj][kk]));
						fprintf(fmt,"\n");
					}
				}
				fclose(fmt);
				return 0;
*/				
				/***************************************************************************************************************/
				/*  STARTING FILTERS LOOP: for each filter                                                                     */
				/***************************************************************************************************************/			
				for (ll=0; ll<NFILTER; ll++)
				{
//					start = clock();
#if LOOP_COMP_OFF == 0					
					marsisazproc(fband[kk] + f0[jj][ii], 
								 Vr[ii], 
								 Vt[ii], 
								 theta[ii], 
								 ll-1, 
								 PRF, 
								 (MARSIS_double complex)(NA[jj][ii]), 
								 X_lon, 
								 Y_lat, 
								 Z_mola, 
								 x0[ii], 
								 y0[ii], 
								 z0[ii], 
								 Vx0[ii], 
								 Vy0[ii], 
								 Vz0[ii],
								 G,
								 N_I[ii],
								 N_J[ii]);
#endif					
//					end = clock();
//					elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
//					printf("Elapsed time (marsisazproc): %g \n", elapsed);					
/*					FILE *fmt;
					fmt = fopen("MatrixG_C.dat","wt");
					for(jj = 0; jj < N_I[ii]; jj++)
					{ 
						for(kk = 0; kk < N_J[ii]; kk++)	
						{
							fprintf(fmt,"row=%d col=%d %lg %lg",jj+1,kk+1,creal(G[jj][kk]),cimag(G[jj][kk]));
							fprintf(fmt,"\n");
						}
					}
					fclose(fmt);
					return 0;
*/					
#if LOOP_COMP_OFF == 0
					espectr=0;
					for(mm=0; mm<N_I[ii]; mm++)
						for(nn=0; nn<N_J[ii]; nn++)
							espectr += E[mm][nn]*cabs(G[mm][nn]);
					espectr *= phase[iband[kk]];
//					printf("%8.20f\t",		creal(espectr));							
//					printf("%8.20f * i\n",	cimag(espectr));		
					
					Respectr[ll][iband[kk]] = (double)(creal(espectr));
					Iespectr[ll][iband[kk]] = (double)(cimag(espectr));
#endif
#if LOOP_COMP_OFF == 1
					Respectr[ll][iband[kk]] = (double)( ii+1);
					Iespectr[ll][iband[kk]] = (double)(-ii-1);
#endif					
				}
				/***************************************************************************************************************/
				/*  END OF FILTERS LOOP: for each filter                                                                       */
				/***************************************************************************************************************/							
			}
			/***************************************************************************************************************/
			/*  END OF IN BANDS FREQUENCIES LOOP: for each frequency in the band to be simulated                           */
			/***************************************************************************************************************/	
//#if MPI == 1			
//			MPI_File_seek(MPI_fout, ii*RecordSize+RecHeadSize+jj*RecDataSize, MPI_SEEK_SET);
//#endif
#if OMP == 1
#pragma omp single
{
#endif			
			for (ll=0; ll<NFILTER; ll++)
			{
#if MPI == 1
	#if BIG_ENDIAN_COMP == 1
				sb_marsis_MPI_File_iwrite(MPI_fout, Respectr[ll], NS, &request);
	#else
				MPI_File_iwrite(MPI_fout, Respectr[ll], NS, MPI_DOUBLE, &request);
	#endif
				MPI_Wait( &request, &status );
	#if BIG_ENDIAN_COMP == 1
				sb_marsis_MPI_File_iwrite(MPI_fout, Iespectr[ll], NS, &request);				
	#else				
				MPI_File_iwrite(MPI_fout, Iespectr[ll], NS, MPI_DOUBLE, &request);
	#endif
				MPI_Wait( &request, &status );
#else
	#if BIG_ENDIAN_COMP == 1
				sb_marsis_fwrite(Respectr[ll],NS,fout);
				sb_marsis_fwrite(Iespectr[ll],NS,fout);
	#else					
				fwrite(Respectr[ll],sizeof(double),NS,fout);
				fwrite(Iespectr[ll],sizeof(double),NS,fout);
	#endif				
#endif
#if MULTI_OUT_FILE == 1
	#if BIG_ENDIAN_COMP == 1
				sb_marsis_fwrite(Respectr[ll],NS,PointFile);
				sb_marsis_fwrite(Iespectr[ll],NS,PointFile);
	#else					
				fwrite(Respectr[ll],sizeof(double),NS,PointFile);
				fwrite(Iespectr[ll],sizeof(double),NS,PointFile);
	#endif
#endif				
			}
#if OMP == 1	
}
#endif
		}
		/***************************************************************************************************************/
		/*  END OF BANDS LOOP: for each band to be simulated                                                           */
		/***************************************************************************************************************/	

		/* FACET2ND MARSISAZPROC Memory allocation *****/
		for(jj = 0; jj < N_I_max; jj++)
		{
			free(E[jj]);
			free(G[jj]);
		}	
		free(E);
		free(G);
		/***********************************************/
	
#if OMP == 1	
}		
#endif
#if MULTI_OUT_FILE == 1
		fclose(PointFile);
#endif				
	}
	/***************************************************************************************************************/
	/*  END OF THE MAIN LOOP: for each echo to be simulated                                                        */
	/***************************************************************************************************************/
//	fclose(fcheck);	
	
	/* MOLA Memory deallocation ********************/
	for(jj = 0; jj < N_I_max; jj++)
	{
		free(Z_mola[jj]);
		free(X_lon[jj]);
		free(Y_lat[jj]);
	}
	free(Z_mola);
	free(X_lon);
	free(Y_lat);
	/***********************************************/
	/* SURFASPECT Memory deallocation **************/
	for(jj = 0; jj < N_I_max; jj++)
	{
		free(la[jj]);
		free(lb[jj]);
		free(Ux[jj]);
		free(Uy[jj]);
		free(Uz[jj]);
		free(R[jj]);
	}
	free(la);
	free(lb);
	free(Ux);
	free(Uy);
	free(Uz);
	free(R);
	/***********************************************/
	/* FACET2ND MARSISAZPROC Memory allocation *****
	for(jj = 0; jj < N_I_max; jj++)
	{
		free(E[jj]);
		free(G[jj]);
	}	
	free(E);
	free(G);*/
	/***********************************************/
	/***** Memory deallocation *********************/
	free(comperffilt);
	free(t);
	free(f);
	free(iband);
	free(fband);
	free(ostline);
	for (jj=0; jj<NBANDS; jj++)
	{
		free(f0[jj]);
		free(NA[jj]);
	}
	free(theta);
	free(frameid);
	free(scetfw);
	free(scetff);
	free(Vt);
	free(Vr);
	free(x0);
	free(y0);
	free(z0);
	free(alt0);
	free(Vx0);
	free(Vy0);
	free(Vz0);
	free(lon_0);
	free(lat_0);
	free(lon_w); 
	free(lon_e);
	free(lat_s);
	free(lat_n);
	free(N_I);
	free(N_J);
	free(Relem);
	free(ProcElem_MPI);
	
	free(SimDataFile);
	free(EdrSimFile);
	free(MolaPath);
	/***********************************************/	
#if OMP == 1
	end = omp_get_wtime();
	elapsed = ((double) (end - start));
#else
	end = clock();
	elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
#endif
	
    
#if MPI == 1
	printf("Process Id: %d - Elapsed time: %g \n", commRank, elapsed);			
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&MPI_fout);
#else
    printf("Elapsed time: %g \n", elapsed);			
	fclose(fout);
#endif
	
#if MPI == 1	
	MPI_Finalize();
#endif	
    return 0;
}
