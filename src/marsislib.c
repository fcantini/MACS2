
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "setParallel.h"
#if MPI == 1
#include <mpi.h>
#endif
#include "marsislib.h"

extern char *MolaPath;

/***********************************************************************************/
/*          FUNCTIONS little_endian and big_endian                                 */
/*          Hardware endian detection                                              */
/*                                                                                 */
/*                           Last update 14/09/2012                                */
/***********************************************************************************/
int little_endian(){
    int x = 1;
    return *(char*)&x;
}
int big_endian(){
    return !little_endian();
}

/***********************************************************************************/
/*          FUNCTION readmarsissimrecords                                          */
/*          reads data length a MARSIS orbit simulation                            */
/*                                                                                 */
/*                           Last update 13/09/2012                                */
/***********************************************************************************/

long int readmarsissimrecords(FILE *fin, int RecordBytes )
{
	long int FileBytes, FileRecords;
	
	fseek( fin, 0, SEEK_END);
	FileBytes = ftell(fin);
	FileRecords = FileBytes / RecordBytes;
	
	if ((FileBytes % RecordBytes) != 0)
	{
		fclose(fin);
		fprintf(stderr, "ReadMarsisSimData:FractionalNumberOfRecords - The simulation data file contains a non integer number of records.");
		return -1;
	}
	return FileRecords;
}

/***********************************************************************************/
/*          FUNCTION readmarsissimdata                                             */
/*          reads necessary data for a MARSIS orbit simulation                     */
/*                                                                                 */
/*                           Last update 09/04/2014                                */
/***********************************************************************************/
int readmarsissimdata(FILE *fin,
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
					  double *lat_0)

{
	long int ii;
	int jj;	
		
	fseek( fin, 0, SEEK_SET);
	for (ii=0; ii<FileRecords; ii++)
	{
		fread(&ostline[ii], sizeof(double), 1, fin);
		for (jj=0; jj<Nbands; jj++)
			fread(&f0[jj][ii], sizeof(double), 1, fin);
		fread(&theta[ii],	sizeof(double), 1, fin);
		fread(&frameid[ii], sizeof(double), 1, fin);
		fread(&scetfw[ii],	sizeof(double), 1, fin);
		fread(&scetff[ii],	sizeof(double), 1, fin);
		fread(&Vt[ii],		sizeof(double), 1, fin);
		fread(&Vr[ii],		sizeof(double), 1, fin);
		for (jj=0; jj<Nbands; jj++)
			fread(&NA[jj][ii], sizeof(double), 1, fin);
		fread(&x0[ii],		sizeof(double), 1, fin);
		fread(&y0[ii],		sizeof(double), 1, fin);
		fread(&z0[ii],		sizeof(double), 1, fin);
		fread(&alt0[ii],	sizeof(double), 1, fin);
		fread(&Vx0[ii],		sizeof(double), 1, fin);
		fread(&Vy0[ii],		sizeof(double), 1, fin);
		fread(&Vz0[ii],		sizeof(double), 1, fin);
		fread(&lon_0[ii],	sizeof(double), 1, fin);
		fread(&lat_0[ii],	sizeof(double), 1, fin);
	}

	return 0;
}

/***********************************************************************************/
/*          FUNCTION fftvars                                                       */
/*          Computes vectors containing the correct time and frequency values for  */
/*          FFT computation                                                        */
/*                                                                                 */
/*                           Last update 26/09/2012                                */
/***********************************************************************************/

int fftvars(MARSIS_double fs, long int Ns, MARSIS_double *t, MARSIS_double *f)
{
	MARSIS_double dt, df;
	int ii;
	
	if (fs <= 0)
	{
		fprintf(stderr, "Fftvars:NegativeSamplingFrequency - The sampling frequency for the FFT vector must be a positive number.");
		return -1;
	}
		
	if (Ns <= 0)
	{
		fprintf(stderr, "FFTvars:NonIntegerNumberOfSamples - The number of samples for the FFT is not a positive integer.");
		return -1;
	}
	
	
	dt = 1/fs;
	for (ii=0; ii<Ns; ii++)
		t[ii] = dt*ii;
	
	if (Ns%2 == 0)
	{
		df = fs/Ns;
		for(ii=0; ii<Ns/2; ii++)
			f[ii] = df*(ii);
		for(ii=Ns/2; ii<Ns; ii++)
			f[ii] = df*( ii-Ns);		
	}
	else
	{
		df = fs/(Ns-1);
		for(ii=0; ii<(Ns+1)/2; ii++)
			f[ii] = df*(ii);
		for(ii=(Ns+1)/2; ii<Ns; ii++)
			f[ii] = df*(ii-Ns+1);		
	}
    
	return 0;
}
#if BIG_ENDIAN_COMP == 1
/***********************************************************************************/
/*          FUNCTION sb_marsis_fread                                               */
/*          Swap bytes MARSIS readind file function                                */
/*                                                                                 */
/*                           Last update 30/11/2012                                */
/***********************************************************************************/
double sb_marsis_fread(FILE *fin)
{
	union 
	{
		double d;
		MARSIS_8 l;
	} marsis_data;
	double t;
	
	fread(&t, sizeof(double), 1, fin);
	marsis_data.d = t;
	marsis_data.l = htole64(marsis_data.l);
	
	return marsis_data.d;
}


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
						 double *lat_0)

{
	long int ii;
	int jj;	
	
	fseek( fin, 0, SEEK_SET);
	for (ii=0; ii<FileRecords; ii++)
	{
		ostline[ii] = sb_marsis_fread(fin);
		
		for (jj=0; jj<Nbands; jj++)
			f0[jj][ii] = sb_marsis_fread(fin);
		theta[ii] = sb_marsis_fread(fin);
		frameid[ii] = sb_marsis_fread(fin);
		scetfw[ii] = sb_marsis_fread(fin);
		scetff[ii] = sb_marsis_fread(fin);
		Vt[ii] = sb_marsis_fread(fin);
		Vr[ii] = sb_marsis_fread(fin);
		for (jj=0; jj<Nbands; jj++)
			NA[jj][ii] = sb_marsis_fread(fin);
		x0[ii] = sb_marsis_fread(fin);
		y0[ii] = sb_marsis_fread(fin);
		z0[ii] = sb_marsis_fread(fin);
		alt0[ii] = sb_marsis_fread(fin);
		Vx0[ii] = sb_marsis_fread(fin);
		Vy0[ii] = sb_marsis_fread(fin);
		Vz0[ii] = sb_marsis_fread(fin);
		lon_0[ii] = sb_marsis_fread(fin);
		lat_0[ii] = sb_marsis_fread(fin);
	}
	
	return 0;
}

/***********************************************************************************/
/*          FUNCTION sb_marsis_fwrite                                              */
/*          Swap bytes fwrite                                                      */
/*                                                                                 */
/*                           Last update 30/11/2012                                */
/***********************************************************************************/
int sb_marsis_fwrite(double *OutBuf, int ns, FILE *fout)
{
	union 
	{
		double d;
		MARSIS_8 l;
	} marsis_data;
	
	double SwapBuf[ns];	
	int ii;
	
	for (ii = 0; ii < ns; ii++)
	{
		marsis_data.d = OutBuf[ii];
		marsis_data.l = htole64(marsis_data.l);
		SwapBuf[ii] = marsis_data.d;
	}

	fwrite(SwapBuf,sizeof(double),ns,fout);
	
	return 0;
}

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
					   double *lat_0)
{
	int jj;
	sb_marsis_fwrite(&ostline[index], 1, fout);
	for (jj=0; jj<Nbands; jj++)
		sb_marsis_fwrite(&f0[jj][index], 1, fout);
	sb_marsis_fwrite(&theta[index], 1, fout);
	sb_marsis_fwrite(&frameid[index], 1, fout);
	sb_marsis_fwrite(&scetfw[index], 1, fout);
	sb_marsis_fwrite(&scetff[index], 1, fout);
	sb_marsis_fwrite(&Vt[index], 1, fout);
	sb_marsis_fwrite(&Vr[index], 1, fout);
	for (jj=0; jj<Nbands; jj++)
		sb_marsis_fwrite(&NA[jj][index], 1, fout);
	sb_marsis_fwrite(&x0[index], 1, fout);
	sb_marsis_fwrite(&y0[index], 1, fout);
	sb_marsis_fwrite(&z0[index], 1, fout);
	sb_marsis_fwrite(&alt0[index], 1, fout);
	sb_marsis_fwrite(&Vx0[index], 1, fout);
	sb_marsis_fwrite(&Vy0[index], 1, fout);
	sb_marsis_fwrite(&Vz0[index], 1, fout);
	sb_marsis_fwrite(&lon_0[index], 1, fout);
	sb_marsis_fwrite(&lat_0[index], 1, fout);
	
	return 0;
}

#endif
#if MPI == 1
#if BIG_ENDIAN_COMP == 1
/***********************************************************************************/
/*          FUNCTION sb_marsis_MPI_File_iwrite                                     */
/*          Swap bytes MPI_File_iwrite                                             */
/*                                                                                 */
/*                           Last update 30/11/2012                                */
/***********************************************************************************/
int sb_marsis_MPI_File_iwrite(MPI_File MPI_fout, double *OutBuf, int ns, MPI_Request *request)
{
	union 
	{
		double d;
		MARSIS_8 l;
	} marsis_data;
	
	double SwapBuf[ns];
	int ii;
	
	for (ii = 0; ii < ns; ii++)
	{
		marsis_data.d = OutBuf[ii];
		marsis_data.l = htole64(marsis_data.l);
		SwapBuf[ii] = marsis_data.d;
	}
	
	MPI_File_iwrite(MPI_fout, SwapBuf, ns, MPI_DOUBLE, request);
	
	return 0;
}
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
							  double *lat_0)
{
	int jj;
	MPI_Status status;
	
	MPI_File_seek(MPI_fout, offset, MPI_SEEK_SET);
	
	sb_marsis_MPI_File_iwrite(MPI_fout, &ostline[index], 1, request);
	MPI_Wait( request, &status );
	for (jj=0; jj<Nbands; jj++)
	{
		sb_marsis_MPI_File_iwrite(MPI_fout, &f0[jj][index],	1, request);
		MPI_Wait( request, &status );
	}
	sb_marsis_MPI_File_iwrite(MPI_fout, &theta[index],	1, request);
	MPI_Wait( request, &status );
	sb_marsis_MPI_File_iwrite(MPI_fout, &frameid[index],1, request);
	MPI_Wait( request, &status );
	sb_marsis_MPI_File_iwrite(MPI_fout, &scetfw[index],	1, request);
	MPI_Wait( request, &status );
	sb_marsis_MPI_File_iwrite(MPI_fout, &scetff[index],	1, request);
	MPI_Wait( request, &status );
	sb_marsis_MPI_File_iwrite(MPI_fout, &Vt[index],		1, request);
	MPI_Wait( request, &status );
	sb_marsis_MPI_File_iwrite(MPI_fout, &Vr[index],		1, request);
	MPI_Wait( request, &status );
	for (jj=0; jj<Nbands; jj++)
	{
		sb_marsis_MPI_File_iwrite(MPI_fout, &NA[jj][index],	1, request);
		MPI_Wait( request, &status );
	}
	sb_marsis_MPI_File_iwrite(MPI_fout, &x0[index],		1, request);
	MPI_Wait( request, &status );
	sb_marsis_MPI_File_iwrite(MPI_fout, &y0[index],		1, request);
	MPI_Wait( request, &status );
	sb_marsis_MPI_File_iwrite(MPI_fout, &z0[index],		1, request);
	MPI_Wait( request, &status );
	sb_marsis_MPI_File_iwrite(MPI_fout, &alt0[index],	1, request);
	MPI_Wait( request, &status );
	sb_marsis_MPI_File_iwrite(MPI_fout, &Vx0[index],	1, request);
	MPI_Wait( request, &status );
	sb_marsis_MPI_File_iwrite(MPI_fout, &Vy0[index],	1, request);
	MPI_Wait( request, &status );
	sb_marsis_MPI_File_iwrite(MPI_fout, &Vz0[index],	1, request);
	MPI_Wait( request, &status );
	sb_marsis_MPI_File_iwrite(MPI_fout, &lon_0[index],	1, request);
	MPI_Wait( request, &status );
	sb_marsis_MPI_File_iwrite(MPI_fout, &lat_0[index],	1, request);
	MPI_Wait( request, &status );
	
	return 0;
}

#endif
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
						   double *lat_0)
{
	int jj;
	MPI_Status status;
	
	MPI_File_seek(MPI_fout, offset, MPI_SEEK_SET);
	
	MPI_File_iwrite(MPI_fout, &ostline[index],	1, MPI_DOUBLE, request);
	MPI_Wait( request, &status );
	for (jj=0; jj<Nbands; jj++)
	{
		MPI_File_iwrite(MPI_fout, &f0[jj][index],	1, MPI_DOUBLE, request);
		MPI_Wait( request, &status );
	}
	MPI_File_iwrite(MPI_fout, &theta[index],	1, MPI_DOUBLE, request);
	MPI_Wait( request, &status );
	MPI_File_iwrite(MPI_fout, &frameid[index],	1, MPI_DOUBLE, request);
	MPI_Wait( request, &status );
	MPI_File_iwrite(MPI_fout, &scetfw[index],	1, MPI_DOUBLE, request);
	MPI_Wait( request, &status );
	MPI_File_iwrite(MPI_fout, &scetff[index],	1, MPI_DOUBLE, request);
	MPI_Wait( request, &status );
	MPI_File_iwrite(MPI_fout, &Vt[index],		1, MPI_DOUBLE, request);
	MPI_Wait( request, &status );
	MPI_File_iwrite(MPI_fout, &Vr[index],		1, MPI_DOUBLE, request);
	MPI_Wait( request, &status );
	for (jj=0; jj<Nbands; jj++)
	{
		MPI_File_iwrite(MPI_fout, &NA[jj][index],	1, MPI_DOUBLE, request);		
		MPI_Wait( request, &status );
	}
	MPI_File_iwrite(MPI_fout, &x0[index],		1, MPI_DOUBLE, request);
	MPI_Wait( request, &status );
	MPI_File_iwrite(MPI_fout, &y0[index],		1, MPI_DOUBLE, request);
	MPI_Wait( request, &status );
	MPI_File_iwrite(MPI_fout, &z0[index],		1, MPI_DOUBLE, request);
	MPI_Wait( request, &status );
	MPI_File_iwrite(MPI_fout, &alt0[index],		1, MPI_DOUBLE, request);
	MPI_Wait( request, &status );
	MPI_File_iwrite(MPI_fout, &Vx0[index],		1, MPI_DOUBLE, request);
	MPI_Wait( request, &status );
	MPI_File_iwrite(MPI_fout, &Vy0[index],		1, MPI_DOUBLE, request);
	MPI_Wait( request, &status );
	MPI_File_iwrite(MPI_fout, &Vz0[index],		1, MPI_DOUBLE, request);
	MPI_Wait( request, &status );
	MPI_File_iwrite(MPI_fout, &lon_0[index],	1, MPI_DOUBLE, request);
	MPI_Wait( request, &status );
	MPI_File_iwrite(MPI_fout, &lat_0[index],	1, MPI_DOUBLE, request);
	MPI_Wait( request, &status );
	
	return 0;
}
#endif
/***********************************************************************************/
/*          FUNCTION writemarsissimdata                                            */
/*          Write simulation parameters in the output file                         */
/*                                                                                 */
/*                           Last update 09/04/2014                                */
/***********************************************************************************/
int writemarsissimdata(FILE *fout,
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
					  double *lat_0)
{
	int jj;
	fwrite(&ostline[index], sizeof(double), 1, fout);
	for (jj=0; jj<Nbands; jj++)
		fwrite(&f0[jj][index], sizeof(double), 1, fout);
	fwrite(&theta[index],	sizeof(double), 1, fout);
	fwrite(&frameid[index], sizeof(double), 1, fout);
	fwrite(&scetfw[index],	sizeof(double), 1, fout);
	fwrite(&scetff[index],	sizeof(double), 1, fout);
	fwrite(&Vt[index],		sizeof(double), 1, fout);
	fwrite(&Vr[index],		sizeof(double), 1, fout);
	for (jj=0; jj<Nbands; jj++)
		fwrite(&NA[jj][index], sizeof(double), 1, fout);		
	fwrite(&x0[index],		sizeof(double), 1, fout);
	fwrite(&y0[index],		sizeof(double), 1, fout);
	fwrite(&z0[index],		sizeof(double), 1, fout);
	fwrite(&alt0[index],		sizeof(double), 1, fout);
	fwrite(&Vx0[index],		sizeof(double), 1, fout);
	fwrite(&Vy0[index],		sizeof(double), 1, fout);
	fwrite(&Vz0[index],		sizeof(double), 1, fout);
	fwrite(&lon_0[index],	sizeof(double), 1, fout);
	fwrite(&lat_0[index],	sizeof(double), 1, fout);
	
	return 0;
}

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
				 long int *N_J)
{
	int ppd;
	double n;
	long int ii, n_rows, n_cols, i_min, i_max, j_min, j_max;
	MARSIS_double lat, lon, lat_max;
	MARSIS_double x_0, y_0, x_min, y_min, x_max, y_max;
	MARSIS_double r_0, delta;
	
	if (fabs(lat_0) < 55) /* non sono implementati tutti i casi ma solo ppd=128 (chiamata di bestmolagrid) */
	{
		lonlatbox( lon_0, lat_0, R_circle, R_sphere, lon_w, lon_e, lat_s, lat_n);
		ppd = 128;
		lat_max     = 88.0;
		n_rows      = lat_max * 2 * ppd;
		n_cols      =         360 * ppd;
		
		*N_I = 0;
		for (ii=0; ii<n_rows; ii++)
		{
			lat = lat_max - (1.0 / (MARSIS_double)(ppd)) * (MARSIS_double)(ii+1) + (1.0 / (MARSIS_double)( 2 * ppd ));
			if ((lat >= *lat_s) && (lat <= *lat_n))
				(*N_I)++;
		}
		
		*N_J = 0;
		if (*lon_w <= *lon_e)
		{
			for (ii=0; ii<n_cols; ii++)
			{
				lon = (1.0 / (MARSIS_double)(ppd)) * (MARSIS_double)(ii+1) - (1.0 / (MARSIS_double)( 2 * ppd ));
				if ((lon >= *lon_w) && (lon <= *lon_e))
					(*N_J)++;
			}			
		}
		else
		{
			for (ii=0; ii<n_cols; ii++)
			{
				lon = (1.0 / (MARSIS_double)(ppd)) * (MARSIS_double)(ii+1) - (1.0 / (MARSIS_double)( 2 * ppd ));
				if ((lon >= *lon_w))
					(*N_J)++;
				if ((lon <= *lon_e))
					(*N_J)++;				
			}						
		}
	}
	else
	{
		if (fabs(lat_0) < 70) /* Implementato solo metodo RADIUS */
		{
			ppd = 128;
			n   = 10240;
		}
		else
		{
			if (fabs(lat_0) < 80)
			{
				ppd = 256;
				n   = 11520;
			}
			else
			{
				ppd = 512;
				n   = 12288;
			}
		}
		
		if (lat_0>=0)
		{
			r_0 = 360 / M_PI * MARSIS_tan( ( M_PI/2 - M_PI / 180 * lat_0 ) / 2 );
			x_0 = r_0 * MARSIS_cos( M_PI / 180 * lon_0 );
		}
		else
		{
			r_0 = 360 / M_PI * MARSIS_tan( ( M_PI/2 + M_PI / 180 * lat_0 ) / 2 );
			x_0 = -r_0 * MARSIS_cos( M_PI / 180 * lon_0 );
		}
		
		y_0 = r_0 * MARSIS_sin( M_PI / 180 * lon_0 );
	
		delta = R_circle / R_sphere * 180/M_PI;
	
		x_min = x_0 - delta;
		x_max = x_0 + delta;
		y_min = y_0 - delta;
		y_max = y_0 + delta;
		
		i_min = round( ppd * x_min + n/2 + 0.5 );
		i_max = round( ppd * x_max + n/2 + 0.5 );
		j_min = round( ppd * y_min + n/2 + 0.5 );
		j_max = round( ppd * y_max + n/2 + 0.5 );
		
		*N_I = i_max-i_min+1;
		*N_J = j_max-j_min+1;
	}

	return 0;
}


/***********************************************************************************/
/*          FUNCTION bestmolagrid                                                  */
/*          Best MOLA set selection                                                */
/*                                                                                 */
/*                           Last update 26/09/2012                                */
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
				 long int N_j)
{
	double theta_circle;
	int ppd;
//	printf("lat_0 =%g\n",lat_0);
	if (fabs(lat_0) < 55)
	{
//		lonlatbox( lon_0, lat_0, R_circle, R_sphere, &lon_w, &lon_e, &lat_s, &lat_n); //trasferito in molagridsize
		molagrid( lon_w, lon_e, lat_s, lat_n, 1, 128.0, mola, lon, lat, N_i, N_j);// quella matlab restituisce 2 vettori. Devono essere  matrici
//		printf("molagrid\n");
	}
	else
	{
		if (fabs(lat_0) < 70)
			ppd = 128;
		else
		{
			if (fabs(lat_0) < 80)
				ppd = 256;
			else
				ppd = 512;
		}
		theta_circle = R_circle / R_sphere * 180/M_PI;
		molapolargrid( lon_0, lat_0, theta_circle, theta_circle, 1, (double)(ppd), mola, lon, lat, N_i, N_j );
//		printf("molapolargrid\n");		
		
	}
	
	
	return 0;
}

/***********************************************************************************/
/*          FUNCTION lonlatbox                                                     */
/*          Used by bestmolagrid. Finds longitude and latitude                     */
/*          box enclosing a circle of a given radius on a sphere                   */
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
			  MARSIS_double *lat_n)
{
	MARSIS_double delta_lat, delta_lon, a, c;
	
	if ((lon_0 < -180) || (lon_0 > 360))
	{
		fprintf(stderr, "Lonlatbox:LongitudeOutOfRange - Input longitude is %g°.\n Allowed longitude range is [ 0°, 360° ].", lon_0 );
		return -1;
	}	
	if (lon_0 < 0)
	{
		fprintf(stderr, "Lonlatbox:NegativeLongitude - Input longitude is %g°.\n Allowed longitude range is [ 0°, 360° ].",lon_0 );
		lon_0 = lon_0 + 360;
	}
	
	if ((lat_0 < -90) || (lat_0 > 90))
	{
		fprintf(stderr, "Lonlatbox:LatitudeOutOfRange - Input latitude is %g°.\n Allowed latitude range is [ -90°, 90° ].",lat_0 );
		return -1;
	}
	/* Pathologies */
	if (R_circle/R_sphere >= M_PI/2)
	{	
		*lon_w =   0;
		*lon_e = 360;
		*lat_s = -90;
		*lat_n =  90;
		return 0;
	}
	
	a = R_circle / R_sphere;
	delta_lat = 180/M_PI * a;
	
	*lat_s = lat_0 - delta_lat;
	*lat_n = lat_0 + delta_lat;
	
	if (*lat_n >= 90)
	{	
		*lon_w =   0;
		*lon_e = 360;
		*lat_n =  90;
		return 0;
	}
	
	if (*lat_s <= -90)
	{	
		*lon_w =   0;
		*lon_e = 360;
		*lat_s = -90;
		return 0;
	}

	/* use the sine theorem: sin A / sin a = sin C / sin c */
	
	c = M_PI/180 * (90 - fabs(lat_0));
	
	delta_lon = 180/M_PI * MARSIS_asin( MARSIS_sin( a ) / MARSIS_sin( c ) );
	
	*lon_w = lon_0 - delta_lon;
	
	if (*lon_w < 0)
		*lon_w = *lon_w + 360;
		
	*lon_e = lon_0 + delta_lon;
	
	if (*lon_e > 360)
		*lon_e = *lon_e - 360;
	
	return 0;
}

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
			 long int N_j)
{
	char id[2];
	size_t prec;
	int samplebytes;
	MARSIS_double valueoffset;
	FILE *fid[4][4];
//	char filename[100];
	char *filename;	
//	char precision[20];
	MOLA_2 mola_t[N_j];
	
    MARSIS_double lat_max;
	MARSIS_double lat_l, lon_l;
    long int n_rows;
    long int n_cols;
    int i_tiles;
    int j_tiles;
    long int n_rows_tile;
    long int n_cols_tile;
	long int ii,jj,kk,ll;
	long int I_g[N_i], J_g[N_j], I_t[N_i], J_t[N_j];
	long int I_wt, J_wt, J_wt_l, I_wt_l;
	long int offset, skip;
	long int mola_I, mola_J;
	
	filename = malloc((strlen(MolaPath)+30) * sizeof(char));
	
	if (lat_s > lat_n)
	{
		fprintf(stderr,"Molagrid:SwappedLatitudes - Minimum latitude should be smaller than maximum latitude.\nMinimum longitude is %g°, maximum longitude is %g°",lat_s, lat_n );
		return -1;
	}	
	if( (lon_w < 0) || (lon_w > 360) || (lon_e < 0) || (lon_e > 360) || (lat_s  < -90) || (lat_n  >  90))
	{
		fprintf(stderr,"Molagrid:CoordinatesOutOfRange - Allowed longitude range is [ 0°, 360° ], allowed latitude range is [ -90°, 90° ].\nInput longitude range is [ %g°, %g° ], input latitude range is [ %g°, %g° ].",lon_w, lon_e, lat_s, lat_n );
		return -1;
	}
	
	/* Al momento implementanto solo type RADIUS */
    strcpy(id,"r");
//    strcpy(prec,"int16");
	prec = sizeof(MOLA_2);
    samplebytes = 2;
    valueoffset = 3396000;	
	/*********************************************/
	
	/* Al momento implementanto solo ppd=128 *****/
	strcpy(filename,"");
	strcat(filename, MolaPath);
	strcat(filename, "meg");
	strcat(filename, id);
	strcat(filename, "88n000hb.img");
	fid[0][0] = fopen(filename, "rb");
	strcpy(filename,"");
	strcat(filename, MolaPath);
	strcat(filename, "meg");
	strcat(filename, id);
	strcat(filename, "88n090hb.img");	
    fid[0][1] = fopen(filename, "rb");
	strcpy(filename,"");
	strcat(filename, MolaPath);
	strcat(filename, "meg");
	strcat(filename, id);
	strcat(filename, "88n180hb.img");	
    fid[0][2] = fopen(filename, "rb");
	strcpy(filename,"");
	strcat(filename, MolaPath);
	strcat(filename, "meg");
	strcat(filename, id);
	strcat(filename, "88n270hb.img");	
    fid[0][3] = fopen(filename, "rb");
	strcpy(filename,"");
	strcat(filename, MolaPath);
	strcat(filename, "meg");
	strcat(filename, id);
	strcat(filename, "44n000hb.img");	
    fid[1][0] = fopen(filename, "rb");
	strcpy(filename,"");
	strcat(filename, MolaPath);
	strcat(filename, "meg");
	strcat(filename, id);
	strcat(filename, "44n090hb.img");	
    fid[1][1] = fopen(filename, "rb");
	strcpy(filename,"");
	strcat(filename, MolaPath);
	strcat(filename, "meg");
	strcat(filename, id);
	strcat(filename, "44n180hb.img");	
    fid[1][2] = fopen(filename, "rb");
	strcpy(filename,"");
	strcat(filename, MolaPath);
	strcat(filename, "meg");
	strcat(filename, id);
	strcat(filename, "44n270hb.img");	
    fid[1][3] = fopen(filename, "rb");
	strcpy(filename,"");
	strcat(filename, MolaPath);
	strcat(filename, "meg");
	strcat(filename, id);
	strcat(filename, "00n000hb.img");	
    fid[2][0] = fopen(filename, "rb");
	strcpy(filename,"");
	strcat(filename, MolaPath);
	strcat(filename, "meg");
	strcat(filename, id);
	strcat(filename, "00n090hb.img");	
    fid[2][1] = fopen(filename, "rb");
	strcpy(filename,"");
	strcat(filename, MolaPath);
	strcat(filename, "meg");
	strcat(filename, id);
	strcat(filename, "00n180hb.img");	
    fid[2][2] = fopen(filename, "rb");
	strcpy(filename,"");
	strcat(filename, MolaPath);
	strcat(filename, "meg");
	strcat(filename, id);
	strcat(filename, "00n270hb.img");	
    fid[2][3] = fopen(filename, "rb");
	strcpy(filename,"");
	strcat(filename, MolaPath);
	strcat(filename, "meg");
	strcat(filename, id);
	strcat(filename, "44s000hb.img");	
    fid[3][0] = fopen(filename, "rb");
	strcpy(filename,"");
	strcat(filename, MolaPath);
	strcat(filename, "meg");
	strcat(filename, id);
	strcat(filename, "44s090hb.img");	
    fid[3][1] = fopen(filename, "rb");
	strcpy(filename,"");
	strcat(filename, MolaPath);
	strcat(filename, "meg");
	strcat(filename, id);
	strcat(filename, "44s180hb.img");	
    fid[3][2] = fopen(filename, "rb");
	strcpy(filename,"");
	strcat(filename, MolaPath);
	strcat(filename, "meg");
	strcat(filename, id);
	strcat(filename, "44s270hb.img");	
    fid[3][3] = fopen(filename, "rb");
	
    lat_max     = 88;
    n_rows      = lat_max * 2 * ppd;
    n_cols      =         360 * ppd;
    i_tiles     = 4;
    j_tiles     = 4;
    n_rows_tile = n_rows / i_tiles;
    n_cols_tile = n_cols / j_tiles;
	/*********************************************/
	
	/* A check is made to verify that input files have been opened correctly. */
		
	for (ii=0; ii<4; ii++)
		for (jj=0; jj<4; jj++)
			if (fid[ii][jj]==NULL)
			{
				fprintf(stderr, "Molagrid:MissingInputFile - One of the required MOLA EGDR files is either missing or corrupted.\n");
				return -1;
			}
	/* The latitude range is checked again and, if necessary, shrunk to fit
	   within the one covered by the selected MOLA MEGDR data set.	*/
	if ((lat_s > lat_max) || (lat_n < -lat_max))
	{
		fprintf(stderr, "Molagrid:LatitudeRangeOutOfBounds - Input latitude range is [ %g°, %g° ].\nAllowed latitude range for the requested resolution is [ -%g°, %g° ].", lat_s, lat_n, lat_max, lat_max );
		return -1;
	}	
	else if (lat_n > lat_max)
	{
		lat_n = lat_max;
		fprintf(stderr, "Molagrid:MaxLatitudeExceeded - Northernmost latitude exceeds the one covered by the requested data set.\nMaximum latitude for which data exist is 88°.");
	}
	else if (lat_s < -lat_max)
	{
		lat_s = -lat_max;
		fprintf(stderr, "Molagrid:MinLatitudeExceeded Southernmost latitude exceeds the one covered by the requested data set.\nMinimum latitude for which data exist is -88°.");
	}

	/* Definition of the output data topographic grid
	
	 The correspondence between the latitude range and a range of row indexes
	 in the MOLA EGDR grid is determined, remembering that latitudes and row
	 indexes in the MOLA EGDR grid are ordered in the opposite sense. Also, the
	 vector containing the latitude corresponding to the rows of the MOLA EGDR
	 output grid is computed. */
	kk=0;
	for (ii=0; ii<n_rows; ii++)
	{	
		lat_l = lat_max - (1.0 / (MARSIS_double)(ppd)) * (MARSIS_double)(ii+1) + (1.0 / (MARSIS_double)( 2 * ppd ));
		if ((lat_l >= lat_s) && (lat_l <= lat_n))
		{
			for(jj=0; jj<N_j; jj++)
				lat[kk][jj] = lat_l;
			I_g[kk] = ii;
			I_t[kk] = I_g[kk]/n_rows_tile + 1;
			kk++;			
		}
	}
	
	kk=0;
	if (lon_w <= lon_e)
	{
		for (ii=0; ii<n_cols; ii++)
		{
			lon_l = (1.0 / (MARSIS_double)(ppd)) * (MARSIS_double)(ii+1) - (1.0 / (MARSIS_double)( 2 * ppd ));
			if ((lon_l >= lon_w) && (lon_l <= lon_e))
			{
				for(jj=0; jj<N_i; jj++)
					lon[jj][kk] = lon_l;
//				J_t = floor ( ( J_g - 1 ) ./ n_cols_tile ) + 1; versione MATLAB
				J_g[kk] = ii;
				J_t[kk] = J_g[kk]/n_cols_tile + 1;
				kk++;
			}
			
		}			
	}
	else
	{
		for (ii=0; ii<n_cols; ii++)
		{
			lon_l = (1.0 / (MARSIS_double)(ppd)) * (MARSIS_double)(ii+1) - (1.0 / (MARSIS_double)( 2 * ppd ));
			if ((lon_l >= lon_w))
			{
				for(jj=0; jj<N_i; jj++)
					lon[jj][kk] = lon_l;
				J_g[kk] = ii;
				J_t[kk] = J_g[kk]/n_cols_tile + 1;
				kk++;				
			}			
			if ((lon_l <= lon_e))
			{
				for(jj=0; jj<N_i; jj++)
					lon[jj][kk] = lon_l;
				J_g[kk] = ii;
				J_t[kk] = J_g[kk]/n_cols_tile + 1;
				kk++;				
			}			
		}						
	}
	
	for (jj = 0; jj<j_tiles; jj++)
	{
		J_wt_l=0;
		for(kk=0; kk<N_j; kk++)
			if (J_t[kk] == (jj+1))
			{
				if (J_wt_l==0)
					mola_J = kk;
				J_wt_l++;	
			}
		if (J_wt_l>0)
			for (ii = 0; ii<i_tiles; ii++)
			{
				I_wt_l=0;
				for(kk=0; kk<N_i; kk++)
					if (I_t[kk] == (ii+1))
					{
						if (I_wt_l==0)
							mola_I = kk;
						I_wt_l++;
					}
				if (I_wt_l>0)
				{
					I_wt = I_g[mola_I] - ( I_t[mola_I] - 1 ) * n_rows_tile;
					J_wt = J_g[mola_J] - ( J_t[mola_J] - 1 ) * n_cols_tile;
					offset = samplebytes * ( ( I_wt ) * n_cols_tile + J_wt );
					skip   = samplebytes * ( n_cols_tile - J_wt_l );
					fseek( fid[ii][jj], offset, SEEK_SET);
					for (ll=0; ll<I_wt_l; ll++)
//					for (ll=0; ll<2; ll++)						
					{
						fread(&mola_t, prec, J_wt_l, fid[ii][jj]);
						fseek( fid[ii][jj], skip, SEEK_CUR);						
						for (kk=0; kk<J_wt_l; kk++)
						{
#if BIG_ENDIAN_COMP == 0							
							mola_t[kk]=(mola_t[kk] << 8) | ((mola_t[kk] >> 8) & 0xFF); //Reversing byte order
#endif
							mola[mola_I+ll][mola_J+kk] = (MARSIS_double)(mola_t[kk]+valueoffset);
						}
					}
				}
			}
	}
	for (ii=0; ii<4; ii++)
		for (jj=0; jj<4; jj++)
			fclose(fid[ii][jj]);
	
	free(filename);
	
	return 0;
}

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
				  long int N_j)
{
    char	id[2];
	char	ppd_s[4];
    size_t	prec;
    int		samplebytes;
    double	scalingfact;
    MARSIS_double valueoffset;
	long int ii, jj;
	long int offset, skip;
	long int i_min, i_max, j_min, j_max;
	long int i_0, j_0;
	MARSIS_double x_0, y_0, x_min, y_min, x_max, y_max;
	MARSIS_double r_0, r, x, x2, y;
	MOLA_4 mola_t[N_j];
	
	union 
	{
		float f;
		MOLA_4 l;
//		char c[sizeof(float)];
	} mola_u;

//	char mola_t[N_j][4];	
	
	long int n;
//	char filename[100];
	char *filename;
	
	
	//FILE *fid[2];
	FILE *fid;
	
	filename = malloc((strlen(MolaPath)+30) * sizeof(char));
	
//	printf("ppd = %g\n",ppd);
	if ((lon_0 < 0) || (lon_0 > 360) || (lat_0  < -90) || (lat_0  >  90))
	{
		fprintf(stderr, "MolaPolarGrid:CoordinatesOutOfRange - Allowed longitude range is [ 0°, 360° ], allowed latitude range is [ -90°, 90° ].\nInput center longitude is %g°, input center latitude range is %g°.'", lon_0, lat_0 );
		return -1;
	}
	/* Al momento implementanto solo type RADIUS */
	strcpy(id,"r");
    prec        = sizeof(long);
    samplebytes = 4;
    scalingfact = 1;
    valueoffset = 3396000;
	/*********************************************/
		
	if ((int)(ppd) == 128)
	{
		n = 10240;
		strcpy(ppd_s, "128");
	}
	else if ((int)(ppd) == 256)
	{
		n = 11520;
		strcpy(ppd_s, "256");
	}		
	else if ((int)(ppd) == 512)
	{
		n = 12288;
		strcpy(ppd_s, "512");
	}		
	else
	{
		fprintf(stderr, "MolaPolarGrid:ResolutionNotFound - Permitted values for input variable _ppd_ are 128, 256 and 512.\n'Input value is %g.", ppd);
		return -1;
	}
	/*	
	sprintf(filename,"meg%s_n_%d.img",id,ppd);
	fid[0] = fopen(filename, "rb");
	sprintf(filename,"meg%s_s_%d.img",id,ppd);
	fid[1] = fopen(filename, "rb");
	*/
	strcpy(filename,"");
	strcat(filename, MolaPath);
	strcat(filename, "meg");
	strcat(filename, id);
	
	if (lat_0 > 0)
		strcat(filename, "_n_");
//		sprintf(filename,"meg%s_n_%d.img",id,ppd);
	else
		strcat(filename, "_s_");		
//		sprintf(filename,"meg%s_s_%d.img",id,ppd);

	strcat(filename, ppd_s);
	strcat(filename, ".img");
	
//	printf("%s\n",filename);
//	printf("lon_0= %g lat_0= %g\n",lon_0, lat_0);
	
	fid = fopen(filename, "rt");
	
	
//	for (ii=0; ii<2; ii++)
		if (fid==NULL)
		{
			fprintf(stderr, "MolaPolarGrid:MissingInputFile - One of the required MOLA EGDR files is either missing or corrupted.\n");
			return -1;
		}
	
	/* The correspondence between the latitude range and a range of row indexes
	   in the MOLA EGDR grid is determined, using formulas presented in
	   http://pds-geosciences.wustl.edu/mgs/mgs-m-mola-5-megdr-l3-v1/mgsl_300x/catalog/dsmap_polar.cat
	   Because polar MOLA EGDR data sets do not cover the entire range of
	   latitudes, indexes for such latitudes will be ignored.	*/
	
	if (lat_0>=0)
	{
		r_0 = 360 / M_PI * MARSIS_tan( ( M_PI/2 - M_PI / 180 * lat_0 ) / 2 );
		x_0 = r_0 * cos( M_PI / 180 * lon_0 );
	}
	else
	{
		r_0 = 360 / M_PI * MARSIS_tan( ( M_PI/2 + M_PI / 180 * lat_0 ) / 2 );
		x_0 = -r_0 * MARSIS_cos( M_PI / 180 * lon_0 );
	}
	
	y_0 = r_0 * sin( M_PI / 180 * lon_0 );
	
	i_0 = round( ppd * x_0 + n/2 + 0.5 );
	j_0 = round( ppd * y_0 + n/2 + 0.5 );
	
	if ((lat_0 == 0) || (i_0 < 1) || (i_0 > n) || (j_0 < 1) || (j_0 > n))
	{
		fprintf(stderr, "MolaPolarGrid:AreaCenterOutsideMap - The input central point is outside the area covered by the requested data set.\nThe coordinates of the central point are [ %g°E, %g°N ].", lon_0, lat_0 );
		return -1;
	}

	x_min = x_0 - delta_x;
	x_max = x_0 + delta_x;
	y_min = y_0 - delta_y;
	y_max = y_0 + delta_y;
	
	i_min = round( ppd * x_min + n/2 + 0.5 );
	i_max = round( ppd * x_max + n/2 + 0.5 );
	j_min = round( ppd * y_min + n/2 + 0.5 );
	j_max = round( ppd * y_max + n/2 + 0.5 );
	
	/* Index values out of range are brought back within range */
	
	if ((i_min < 1) || (i_min > n) || (i_max < 1) || (i_max > n) || (j_min < 1) || (j_min > n) || (j_max < 1) || (j_max > n))
		fprintf(stderr, "MolaPolarGrid:PartOfAreaOutsideMap - Part of the input zone is outside the area covered by the requested data set.");
	
		
	if (i_min < 1)
		i_min = 1;
	if (i_max < 1)
		i_max = 1;
	if (j_min < 1)
		j_min = 1;
	if (j_max < 1)
		j_max = 1;
	
	if (i_min > n)
		i_min = n;
	if (i_max > n)
		i_max = n;
	if (j_min > n)
		j_min = n;
	if (j_max > n)
		j_max = n;
	
	/* Longitudes and latitudes corresponding to index values are computed using
	   formulas derived from the inversion of those presented in
	   http://pds-geosciences.wustl.edu/mgs/mgs-m-mola-5-megdr-l3-v1/mgsl_300x/catalog/dsmap_polar.cat */	
	if (lat_0 > 0)
	{
		for (ii = 0; ii<N_i; ii++)
		{
			x = ((MARSIS_double)(ii+i_min) -n/2.0 -0.5)/ppd;
			x2 = x*x;
			for (jj = 0; jj<N_j; jj++)
			{
				y = ((MARSIS_double)(jj+j_min) -n/2.0 -0.5)/ppd;
				r = MARSIS_sqrt( x2 + y*y );
				lon[ii][jj] = MARSIS_atan2( y, x ) * 180.0/M_PI;
				if (lon[ii][jj] < 0)
					lon[ii][jj]+=360;
				lat[ii][jj] =  90.0 - 2.0 * MARSIS_atan( r * M_PI/360.0 ) * 180.0/M_PI;				
			}
		}
	}
	else
	{
		for (ii = 0; ii<N_i; ii++)
		{
			x = ((MARSIS_double)(ii+i_min) -n/2.0 -0.5)/ppd;
			x2 = x*x;
			for (jj = 0; jj<N_j; jj++)
			{
				y = ((MARSIS_double)(jj+j_min) -n/2.0 -0.5)/ppd;
				r = MARSIS_sqrt( x2 + y*y );
				lon[ii][jj] = 180.0 - MARSIS_atan2( y, x ) * 180.0/M_PI;
				if (lon[ii][jj] < 0)
					lon[ii][jj]+=360;
				lat[ii][jj] =  -90.0 + 2.0 * MARSIS_atan( r * M_PI/360.0 ) * 180.0/M_PI;				
			}		
		}
	}
	
	/* Reading of data from input file(s) */
	
	/*The file position indicator is set at the beginning of the section of
	  data to be retrieved, remembering that each MOLA EGDR sample is
	  _samplebytes_ bytes long. Data are to be read one row at the time,
	  retrieving the correct number of samples of each MOLA EGDR tile row, and
	  skipping the portion of the row which falls outside the area of interest. */
	
	
	offset    = samplebytes * ( ( i_min - 1 ) * n + j_min - 1 );

	skip      = samplebytes * ( n - N_j  );
	fseek( fid, offset, SEEK_SET);
//	printf("pos %ld \n",ftell(fid));
	
	for (ii = 0; ii < N_i; ii++)
	{
		fread(mola_t, samplebytes, N_j, fid);	
		fseek(fid, skip, SEEK_CUR);
		
		for (jj=0; jj<N_j; jj++)
		{
			
			mola_u.l = mola_t[jj];
#if BIG_ENDIAN_COMP == 0			
			mola_u.l = ((mola_u.l>>24)&0xff) | ((mola_u.l<<8)&0xff0000) | ((mola_u.l>>8)&0xff00) | ((mola_u.l<<24)&0xff000000);
#endif			
			mola[ii][jj] = scalingfact * (mola_u.f+valueoffset);
		}		
	}
	
	fclose(fid);
	free(filename);
	
	return 0;
}


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
			   MARSIS_double **R)
{
	long int ii,jj,ii_u,ii_d,jj_l,jj_r;
	MARSIS_double Xa, Xb, Xn, Ya, Yb, Yn, Za, Zb, Zn, XR, YR, ZR, n; 
	
	for (ii=0; ii<N_i; ii++)
		for (jj=0; jj<N_j; jj++)
			 {
				 jj_r = (jj<N_j-1) ? jj+1 : jj;
				 jj_l = (jj>0) ? jj-1 : 0;
				 
				 Xa = X[ii][jj_r] - X[ii][jj_l];
				 Ya = Y[ii][jj_r] - Y[ii][jj_l];
				 Za = Z[ii][jj_r] - Z[ii][jj_l];

				 ii_d = (ii<N_i-1) ? ii+1 : ii;
				 ii_u = (ii>0) ? ii-1 : 0;

				 Xb = X[ii_u][jj] - X[ii_d][jj];
				 Yb = Y[ii_u][jj] - Y[ii_d][jj];
				 Zb = Z[ii_u][jj] - Z[ii_d][jj];				 
				 
				 Xn = Ya * Zb - Za * Yb; 
				 Yn = Za * Xb - Xa * Zb;
				 Zn = Xa * Yb - Ya * Xb;
				 
				 XR = x0 - X[ii][jj];
				 YR = y0 - Y[ii][jj];
				 ZR = z0 - Z[ii][jj];
				 
				 la[ii][jj] = MARSIS_sqrt( Xa*Xa + Ya*Ya + Za*Za );
				 Xa = Xa / la[ii][jj];
				 Ya = Ya / la[ii][jj];
				 Za = Za / la[ii][jj];
				 
				 lb[ii][jj] = MARSIS_sqrt( Xb*Xb + Yb*Yb + Zb*Zb );
				 Xb = Xb / lb[ii][jj];
				 Yb = Yb / lb[ii][jj];
				 Zb = Zb / lb[ii][jj];
				 
				 n = MARSIS_sqrt( Xn*Xn + Yn*Yn + Zn*Zn );
				 Xn = Xn / n;
				 Yn = Yn / n;
				 Zn = Zn / n;
				 
				 R[ii][jj] = MARSIS_sqrt( XR*XR + YR*YR + ZR*ZR );
				 XR = XR / R[ii][jj];
				 YR = YR / R[ii][jj];
				 ZR = ZR / R[ii][jj];				 
				 
				 la[ii][jj] = la[ii][jj] / 2;
				 lb[ii][jj] = lb[ii][jj] / 2;				
				 
				 Ux[ii][jj] = XR * Xa + YR * Ya + ZR * Za;
				 Uy[ii][jj] = XR * Xb + YR * Yb + ZR * Zb;
				 Uz[ii][jj] = XR * Xn + YR * Yn + ZR * Zn;				 
			 }
	
	return 0;
}


/***********************************************************************************/
/*         FUNCTION facet2nd                                                       */
/*         Determine the complex reflectivity of all facets                        */
/*                                                                                 */
/*                           Last update 27/09/2012                                */
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
			 int ncoef)
{
	MARSIS_double k, k4, cos_alpha_x, cos_alpha_y, a, b, sqrtR;
	MARSIS_double complex sqrtk/*, meno134*/;
	MARSIS_double complex Mk, kI, I_M_PI;
	int ii, jj;
	
	k = 2.0 * M_PI * f / C; // angular wavenumber, 2 pi / lambda
	k4 = -4.0*k;
	kI = -I*k;
//	meno134 = MARSIS_cpow(-1.0,0.75); 
//	meno134 = -0.70710678118654746 + 0.70710678118654757 * I;
	sqrtk = MARSIS_csqrt(-k);
	Mk = meno134*sqrtk;
	I_M_PI =  I * M_PI;
//	printf("Entering facet2nd\n");
//			for (jj=0; jj<N_j; jj++)
	for (ii=0; ii<N_i; ii++)
		for (jj=0; jj<N_j; jj++)
		{
//			printf("row=%d col=%d \r",ii+1,jj+1);
			cos_alpha_x = Uz[ii][jj] / MARSIS_sqrt(Ux[ii][jj]*Ux[ii][jj] + Uz[ii][jj]*Uz[ii][jj] );
			cos_alpha_y = Uz[ii][jj] / MARSIS_sqrt(Uy[ii][jj]*Uy[ii][jj] + Uz[ii][jj]*Uz[ii][jj] );
			a = Ux[ii][jj]/Uz[ii][jj];
			b = Uy[ii][jj]/Uz[ii][jj];
			
			sqrtR = MARSIS_sqrt(R[ii][jj]);
			
#if MARSIS_erf_alg == 0
			E[ii][jj] = 
			( I_M_PI ) / ( k4 * R[ii][jj] ) * MARSIS_cexp( kI * R[ii][jj] * ( 2 - a*a - b*b ) ) * 
			( erfz(Mk * ( 2.0 * a * R[ii][jj] + la[ii][jj] * cos_alpha_x ) / ( 2.0 * sqrtR ) ) +
			 erfz(Mk * sqrtR * ( -a + la[ii][jj] * cos_alpha_x / ( 2.0 * R[ii][jj] ) ) ) ) *
			( erfz( Mk * ( 2.0 * b * R[ii][jj] + lb[ii][jj] * cos_alpha_y ) / ( 2.0 * sqrtR ) ) +
			 erfz( Mk * sqrtR * ( -b + lb[ii][jj] * cos_alpha_y / ( 2.0 * R[ii][jj] ) ) ) );
#endif
			
#if MARSIS_erf_alg == 1
			E[ii][jj] = 
			( I_M_PI ) / ( k4 * R[ii][jj] ) * MARSIS_cexp( kI * R[ii][jj] * ( 2 - a*a - b*b ) ) * 
			( cerf(Mk * ( 2.0 * a * R[ii][jj] + la[ii][jj] * cos_alpha_x ) / ( 2.0 * sqrtR ), L, coef, ncoef) +
			 cerf(Mk * sqrtR * ( -a + la[ii][jj] * cos_alpha_x / ( 2.0 * R[ii][jj] ) ), L, coef, ncoef )) *
			( cerf( Mk * ( 2.0 * b * R[ii][jj] + lb[ii][jj] * cos_alpha_y ) / ( 2.0 * sqrtR ), L, coef, ncoef ) +
			 cerf( Mk * sqrtR * ( -b + lb[ii][jj] * cos_alpha_y / ( 2.0 * R[ii][jj] ) ), L, coef, ncoef ) );
#endif
			
#if MARSIS_erf_alg == 2			
			E[ii][jj] =
			( I_M_PI ) / ( k4 * R[ii][jj] ) * MARSIS_cexp( kI * R[ii][jj] * ( 2 - a*a - b*b ) ) * 
			( Faddeeva_erf(Mk * ( 2.0 * a * R[ii][jj] + la[ii][jj] * cos_alpha_x ) / ( 2.0 * sqrtR ),0 ) +
			 Faddeeva_erf(Mk * sqrtR * ( -a + la[ii][jj] * cos_alpha_x / ( 2.0 * R[ii][jj] ) ),0 ) ) *
			( Faddeeva_erf( Mk * ( 2.0 * b * R[ii][jj] + lb[ii][jj] * cos_alpha_y ) / ( 2.0 * sqrtR ),0 ) +
			 Faddeeva_erf( Mk * sqrtR * ( -b + lb[ii][jj] * cos_alpha_y / ( 2.0 * R[ii][jj] ) ),0 ) );
#endif			
			
		}
	return 0;
}

/***********************************************************************************/
/*         FUNCTION erfz                                                           */
/*         Error function of complex numbers                                       */
/*                                                                                 */
/*                           Last update 26/09/2012                                */
/***********************************************************************************/
MARSIS_double complex erfz(MARSIS_double complex z)
{
	MARSIS_double complex e, g;
	MARSIS_double Re, Im; 
	
	Re=creal(z);
	Im=cimag(z);
//printf("-> %g + %gi\n",Re,Im);
	if (Im==0)
		return MARSIS_erf(Re);
	
	if ((Re==NAN) || (Im==NAN))
		return NAN;
	// siamo rimasti alla riga 41. Capire bene cosa fanno le righe 45 e 46
	e = MARSIS_erf(Re);
	g = parts(Re, MARSIS_fabs(Im));
	if (Im<0)
		g = MARSIS_conj(g);
	e+=g;

	return e;
}


/***********************************************************************************/
/*         FUNCTION parts                                                          */
/*         Used by erfz: Evaluates all partial functions E(z) to G(z).             */
/*                                                                                 */
/*                           Last update 10/10/2012                                */
/***********************************************************************************/
MARSIS_double complex parts(MARSIS_double Re, MARSIS_double Im)
{
#define N 13
	
	MARSIS_double R2, R3, F, H, Hr, Hi, G, Gr, Gi, n1, n2;
	MARSIS_double complex e2iRI, e, E;
	int nn, mm, M;
	
	R2 = Re*Re;
	e2iRI = MARSIS_cexp(-2.0*Re*Im*I);
	if (Re==0)
		E = 0;
	else
		E=(1-e2iRI)/(2.0*M_PI*Re);
	
	F=0;
	Hr=0;
	Hi=0;
	
	for (nn=1; nn<=N; nn++) //Riga 71 erfz.m (Può essere calcolato a partire dalla precisione del formato numerico)
	{
		H=(MARSIS_double)(nn*nn)/4.0;
		H=MARSIS_exp(-H)/(H + R2);
		F=F + H;
		H=H*MARSIS_exp(-nn*Im);
		Hi=Hi + (MARSIS_double)(nn)/2.0*H;
		Hr=Hr + H;		
	}
	
	e=MARSIS_exp(-R2)*(E + Re*F/M_PI - e2iRI*(Re*Hr+Hi*I)/(2.0*M_PI));
	R3=R2 + 1.837877066409345/*MARSIS_log(2.0*M_PI)*/;
	Gr=0;
	Gi=0;
	M=2*Im;
	nn=(int)(floor(M - N));
	if (nn<1)
		nn=1;
	
	M=(int)(ceil(M + N - nn));

	for(mm=0; mm<=M; mm++)
	{
		n1=(MARSIS_double)(nn)/2.0;
		n2=n1*n1;
		G=MARSIS_exp((MARSIS_double)(nn)*Im - n2 - R3 - MARSIS_log(n2 + R2));
		Gi=Gi - n1*G;
		Gr=Gr + G;
		nn++;
	}
	
	e -= e2iRI*(Re*Gr+Gi*I);
	if (Re==0)
		e=e + (Im/M_PI);
	
	return e;
}


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
				 long int N_j)
{
	MARSIS_double lambda, X_diff, Y_diff, Z_diff, R_diff, V_rel, dt, dt_corr;
	MARSIS_double complex dphi;
	int ii,jj;
	
	lambda = C/f;
	for (ii=0; ii<N_i; ii++)
		for (jj=0; jj<N_j; jj++)
		{
//			start = clock();	
			/* Relative velocity between spacecraft and discrete ground points	*/
			X_diff = x0 - X[ii][jj];
			Y_diff = y0 - Y[ii][jj];
			Z_diff = z0 - Z[ii][jj];
			
			R_diff = MARSIS_sqrt( X_diff*X_diff + Y_diff*Y_diff + Z_diff*Z_diff);
			
			V_rel = (Vx0 * X_diff + Vy0 * Y_diff + Vz0 * Z_diff) / R_diff;
			
			/* Phase difference	*/
			dt = -2.0 * V_rel / (prf * C);
			dt_corr = -( 2.0 * Vr / prf + 2.0 * Vt * theta_s / prf  + m * lambda / NA ) / C;
			
			/* Synthetic aperture pattern */
			dphi = MARSIS_cexp( 2*I * M_PI * ( dt - dt_corr ) * f );
			
			if(dphi == 1)
				dphi = MARSIS_cexp( 2*I * M_PI * EPS ); // to avoid singularities

				
			G[ii][jj] = dphi * ( 1.0 - MARSIS_cpow(dphi,NA)) / ( 1.0 - dphi );

						
		}
	return 0;
}

/***********************************************************************************/
/*         FUNCTION readcomperffilterlen                                           */
/*         Returns the number of lines in complex erf filter text file.            */
/*         Returns -1 if an error occours                                          */
/*                                                                                 */
/*                           Last update 05/10/2012                                */
/***********************************************************************************/
int readcomperffilterlen(FILE *fid)
{
	/* Returns the number of lines in a text file. Returns -1 if an error occours*/
	{
		char line[40];
		int len = 0;
		int cur_pos = ftell(fid);
		if (fseek(fid, 0, SEEK_SET))
			return -1;
		
		while (!feof(fid))
		{
			len++;
			fgets(line, 40, fid);
		}
		if (fseek(fid, cur_pos, SEEK_SET))
			return -1;	
		
		if (strchr(line,'\n')!=0)
			len--;
		
		return len;
	}
}

/***********************************************************************************/
/*         FUNCTION readcomperffilter                                              */
/*         Returns complex erf filter coefficients.                                */
/*                                                                                 */
/*                           Last update 05/10/2012                                */
/***********************************************************************************/
int readcomperffilter(FILE *fid, MARSIS_double *L, MARSIS_double *coef, int Ncoef)
{
	char line[40];
	int ll = Ncoef-1;
	
	//	int cur_pos = ftell(fid);
	if (fseek(fid, 0, SEEK_SET))
		return -1;
	
	fgets(line, 40, fid);
	*L = MARSIS_aton(line);
	
	for(ll=Ncoef-1; ll>=0; ll--)
	{
		fgets(line, 40, fid);
		coef[ll] = MARSIS_aton(line);
	}
	
	return 1;
}

/***********************************************************************************/
/*         FUNCTION cerf                                                           */
/*         Returns erf value with complex argument.                                */
/*                                                                                 */
/*                           Last update 05/10/2012                                */
/***********************************************************************************/
MARSIS_double complex cerf(MARSIS_double complex z, MARSIS_double L, MARSIS_double *coef, int ncoef)
{
	MARSIS_double complex Z,w,a,Zacc,z1,z2;
	MARSIS_double k;
//	MARSIS_double Re, Im; 
	int ii;
//	if ((fabs(creal(z))<.0001) || (fabs(cimag(z))<.0001))
//		printf("%lg %lg \n",fabs(creal(z)),fabs(cimag(z))); 
//	return 0;
	if(creal(z)>0)
	{
		k = -1;
		z1 =  I*z;		
	}
	else
	{
		k = 1;
		z1 =  -I*z;		
	}
//	Re=creal(z);
//	Im=cimag(z);
	
//	if (Im==0)
//		return MARSIS_erf(Re);
	
//	if ((Re==NAN) || (Im==NAN))
//		return NAN;

	z2 = (-I)*z;
	a = (L-I*z1);
	Zacc = Z = (L+I*z1)/a;
	
	w = coef[0];
	for(ii=1; ii<ncoef; ii++)
	{
		w += coef[ii]*Zacc;
		Zacc *= Z;
	}
	
	w = 2*w/(a*a)+(0.564189583547756)/a;
	
	return k*(w*MARSIS_cexp(z2*z2)-1);
}

/***********************************************************************************/
/*         FUNCTION struct_cmp_by_elements                                         */
/*         Comparing function for qsort.                                           */
/*                                                                                 */
/*                           Last update 13/10/2012                                */
/***********************************************************************************/
int struct_cmp_by_elements(const void *a, const void *b)
{ 
    struct record_elem_struct *ia = (struct record_elem_struct *)a;
    struct record_elem_struct *ib = (struct record_elem_struct *)b;
    return (int)(ib->Nelem - ia->Nelem);
}	

/***********************************************************************************/
/*         FUNCTION select_proc_elem                                               */
/*         Selects trajectory points to be simulated by each MPI process.          */
/*                                                                                 */
/*                           Last update 16/10/2012                                */
/***********************************************************************************/
int select_MPI_proc_elem_NO(record_elem *Relem, long int Relem_len, long int *ProcElem_MPI, long int *N_ProcElem_MPI, int N_proc_MPI, int Proc_MPI)
{
	long int acc, ii, jj, kk, N_elem_tot;
	long int selected[Relem_len];
	float elempgroup;
	
	N_elem_tot = 0;
	for (ii=0; ii<Relem_len; ii++)
	{
		selected[ii] = 0;
		N_elem_tot+=Relem[ii].Nelem;
	}
		
	elempgroup = (float)(N_elem_tot)/(float)(N_proc_MPI);
	ii = 0;		
	while (ii<=Proc_MPI)
	{
		acc=0;
		kk=0;
		for (jj=0; jj<Relem_len; jj++)
		{
			if ((acc+Relem[jj].Nelem)<elempgroup)
			{
				if (selected[jj]==0)
				{
					selected[jj]=1;
					acc=acc+Relem[jj].Nelem;
					if (ii==Proc_MPI)
						ProcElem_MPI[kk++]=Relem[jj].Index;
				}
			}        
		}
		if (ii==Proc_MPI)
			*N_ProcElem_MPI = kk;
		ii=ii+1;
	}
	
	return 0;
}


/***********************************************************************************/
/*         FUNCTION select_proc_elem                                               */
/*         Selects trajectory points to be simulated by each MPI process.          */
/*                                                                                 */
/*                           Last update 15/10/2012                                */
/***********************************************************************************/
int select_MPI_proc_elem(record_elem *Relem, long int Relem_len, long int *ProcElem_MPI, long int *N_ProcElem_MPI,  int N_proc_MPI, int Proc_MPI)
{
	long int ii, jj, kk, ll, N_elem_tot, start_lap, min, I_min;
	long int selected[Relem_len];
	long int elemgroupN[N_proc_MPI];
	float elempgroup;
	
	N_elem_tot = 0;
	for (ii=0; ii<Relem_len; ii++)
	{
		selected[ii] = 0;
		N_elem_tot+=Relem[ii].Nelem;
	}
	elempgroup = (float)(N_elem_tot)/(float)(N_proc_MPI);
	
	for (ii=0; ii<N_proc_MPI; ii++)
	{
		elemgroupN[ii] = 0;
	}	
	
	ii = 0; //Processi
	jj = 0;	//punti traiettoria
	kk = 0;
	start_lap = -1;
	while (jj<Relem_len)
	{
		if (elemgroupN[ii]+Relem[jj].Nelem<elempgroup)
		{
			elemgroupN[ii]+=Relem[jj].Nelem;
			if (ii==Proc_MPI)
				ProcElem_MPI[kk++]=Relem[jj].Index;
			selected[jj]=1;
			jj++;
			start_lap=-1;
		}
		else
		{
			if (start_lap == -1)
				start_lap=ii;
		}
		
		ii++;
		if (ii>=N_proc_MPI)
			ii=0;
		if (ii==start_lap)
		{
			min = elemgroupN[0];
			I_min = 0;
			for (ll=1; ll<N_proc_MPI; ll++)
				if (elemgroupN[ll]<min)
				{
					min = elemgroupN[ll];
					I_min = ll;
				}
			elemgroupN[I_min]+=Relem[jj].Nelem;
			if (I_min==Proc_MPI)
				ProcElem_MPI[kk++]=Relem[jj].Index;			
			selected[jj]=1;
			jj++;
		}
	}	
	*N_ProcElem_MPI = kk;
//	for(ii=0;ii<N_proc_MPI;ii++)
//		printf("Elementi[%d]: %d\n",ii,elemgroupN[ii]);
	return 1;
}

/***********************************************************************************/
/*         FUNCTION Faddeeva_w                                                     */
/*         From http://ab-initio.mit.edu/wiki/index.php/Faddeeva_w                 */
/*                                                                                 */
/*         Slightly modified form C++ to C                                         */
/***********************************************************************************/

double complex Faddeeva_erf(double complex z, double relerr)
{
	MARSIS_double complex z1,z2;
	MARSIS_double k;
	
	if(creal(z)>0)
	{
		k = -1;
		z1 =  I*z;		
	}
	else
	{
		k = 1;
		z1 =  -I*z;		
	}
	
	z2 = (-I)*z;
	
	return k*(Faddeeva_w(z1,relerr)*MARSIS_cexp(z2*z2)-1);
}

