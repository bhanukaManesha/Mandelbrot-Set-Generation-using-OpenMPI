//////////////////////////////////////////////////////////////////////////////////////
// mandelbrot.c program: Mandelbort Set Fractal (Color Serial Code Implementation).
// --------------------------------
//  1. Draws Mandelbrot set for Fc(z) = z*z +c
//  using Mandelbrot algorithm ( boolean escape time )
//	This code is modified from the original version as available at:
//	http://rosettacode.org/wiki/Mandelbrot_set#PPM_non_interactive
// -------------------------------         
// 2. Technique of creating ppm file is  based on the code of Claudio Rocchini
// http://en.wikipedia.org/wiki/Image:Color_complex_plot.jpg
// create 24 bit color graphic file ,  portable pixmap file = PPM 
// see http://en.wikipedia.org/wiki/Portable_pixmap
// to see the file use external application ( graphic viewer)
//////////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

// Main program
int main(int argc, char *argv[])
 {

	//  MPI 
	int numtasks, rank;
    int shared_buffer;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Status stat;

	/* screen ( integer) coordinate */
	int iX,iY;
	const int iXmax = 8000; // default
	const int iYmax = 8000; // default

	/* world ( double) coordinate = parameter plane*/
	double Cx, Cy;
	const double CxMin = -2.5;
	const double CxMax = 1.5;
	const double CyMin = -2.0;
	const double CyMax = 2.0;

	/* */
	double PixelWidth = (CxMax - CxMin)/iXmax;
	double PixelHeight = (CyMax - CyMin)/iYmax;

	/* color component ( R or G or B) is coded from 0 to 255 */
	/* it is 24 bit color RGB file */
	const int MaxColorComponentValue = 255; 

	// RGB color array
	static unsigned char color[3];

	// Row RGB array
	static unsigned char row_color[iXmax * 3];

	/* Z = Zx + Zy*i;	Z0 = 0 */
	double Zx, Zy;
	double Zx2, Zy2; /* Zx2 = Zx*Zx;  Zy2 = Zy*Zy  */
	/*  */
	int Iteration;
	const int IterationMax = 1000; // default

	/* bail-out value , radius of circle ;  */
	const double EscapeRadius = 400;
	double ER2 = EscapeRadius * EscapeRadius;
	
	/* Clock information */
	clock_t start, end;
	double cpu_time_used;

	if (rank == 0){
		FILE * fp;
		char *filename = "Mandelbrot.ppm";
		char *comment = "# ";	/* comment should start with # */

		/*create new file,give it a name and open it in binary mode  */
		fp = fopen(filename, "wb"); /* b -  binary mode */

		/*write ASCII header to the file (PPM file format)*/
		fprintf(fp,"P6\n %s\n %d\n %d\n %d\n", comment, iXmax, iYmax, MaxColorComponentValue);

		printf("File: %s successfully opened for writing.\n", filename);
		printf("Computing Mandelbrot Set. Please wait...\n");
		
		// Get current clock time.
		start = clock();


		/* compute and write image data bytes to the file */
		for(iY = 0; iY < iYmax; iY++)
		{	
			// printf("Started Row %i.\n",iY);
			Cy = CyMin + (iY * PixelHeight);
			if (fabs(Cy) < (PixelHeight / 2))
			{
				Cy = 0.0; /* Main antenna */
			}

			MPI_Send(&Cy, sizeof(double), MPI_DOUBLE, 1 + (iY % numtasks - 1), 0, MPI_COMM_WORLD);
			MPI_Recv(&row_color, sizeof(char) * iXmax * 3, MPI_UNSIGNED_CHAR, 1 + (iY % numtasks - 1), 0, MPI_COMM_WORLD, &stat);

			/* write color to the file */
			fwrite(row_color, 1, sizeof(char) * iXmax * 3 , fp);

			// printf("Completed Row %i.\n",iY);

		}
		Cy = 5;
		MPI_Send(&Cy, sizeof(int), MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
		MPI_Send(&Cy, sizeof(int), MPI_DOUBLE, 2, 0, MPI_COMM_WORLD);
		MPI_Send(&Cy, sizeof(int), MPI_DOUBLE, 3, 0, MPI_COMM_WORLD);
		
		// Get the clock current time again
		// Subtract end from start to get the CPU time used.
		end = clock();
		cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

		fclose(fp);


		printf("Completed Computing Mandelbrot Set.\n");
		printf("File: %s successfully closed.\n", filename);
		printf("Mandelbrot computational process time: %lf\n", cpu_time_used);
		MPI_Finalize();
		return 0;


	}else{

		while (1){
			MPI_Recv(&Cy, sizeof(double), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &stat);

			if (Cy == 5){
				break;
			}

			for(iX = 0; iX < iXmax; iX++)
			{
				Cx = CxMin + (iX * PixelWidth);
				/* initial value of orbit = critical point Z= 0 */
				Zx = 0.0;
				Zy = 0.0;
				Zx2 = Zx * Zx;
				Zy2 = Zy * Zy;

				/* */
				for(Iteration = 0; Iteration < IterationMax && ((Zx2 + Zy2) < ER2); Iteration++)
				{
					Zy = (2 * Zx * Zy) + Cy;
					Zx = Zx2 - Zy2 + Cx;
					Zx2 = Zx * Zx;
					Zy2 = Zy * Zy;
				};

				/* compute  pixel color (24 bit = 3 bytes) */
				if (Iteration == IterationMax)
				{
					// Point within the set. Mark it as black
					row_color[iX * 3] = 0;
					row_color[iX * 3 + 1] = 0;
					row_color[iX * 3 + 2] = 0;
				}
				else 
				{
					// Point outside the set. Mark it as white
					double c = 3*log((double)Iteration)/log((double)(IterationMax) - 1.0);
					if (c < 1)
					{
						row_color[iX * 3] = 0;
						row_color[iX * 3 + 1] = 0;
						row_color[iX * 3 + 2] = 255*c;
					}
					else if (c < 2)
					{
						row_color[iX * 3] = 0;
						row_color[iX * 3 + 1] = 255*(c-1);
						row_color[iX * 3 + 2] = 255;
					}
					else
					{
						row_color[iX * 3] = 255*(c-2);
						row_color[iX * 3 + 1] = 255;
						row_color[iX * 3 + 2] = 255;
					}
				}

			}

			MPI_Send(&row_color, sizeof(char) * iXmax * 3, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);

		}

		MPI_Finalize();
		
		



	}

 }