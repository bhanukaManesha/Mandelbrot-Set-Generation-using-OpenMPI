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
#include <string.h>
#include <mpi.h>

#define iXmax 8000 // default
#define iYmax 8000 // default
#define sizeOfppm iXmax * iYmax * 3

// Main program
int main(int argc, char *argv[])
 {
	/* Clock information */
	double start, start_comp, end, end_comp;
	//double cpu_time_used;

	// Get current clock time.
	start = MPI_Wtime();

	//  MPI 
	int numtasks, rank;
    int shared_buffer;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Status stat;

	/* screen ( integer) coordinate */
	int iX,iY;

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
	const int IterationMax = 2000; // default

	/* bail-out value , radius of circle ;  */
	const double EscapeRadius = 400;
	double ER2 = EscapeRadius * EscapeRadius;

	static unsigned char ppm[sizeOfppm];

	int startRow, endRow, remainder;
	int offset;

    startRow = (iYmax / numtasks) * rank;
    endRow = (iYmax / numtasks) * (rank + 1) - 1;

	remainder = iYmax % numtasks;

    if (rank == numtasks -1 ){
        endRow = endRow + remainder;
    }

	printf("Rank %i :: %i :: %i\n", rank, startRow, endRow);

	// printf("rank : %i ->  start : %i -> end %i\n", rank,startRow, endRow);

	unsigned char* p_segment = (unsigned char*) malloc(((iYmax/numtasks) + remainder) * iXmax * 3);

	// Calculate CY Values ?? Can optimize
	double Cy_Ar[sizeof(double) * iYmax];
	/* compute and write image data bytes to the file */
	for(iY = 0; iY < iYmax; iY++)
	{	
		Cy_Ar[iY] = CyMin + (iY * PixelHeight);
		if (fabs(Cy_Ar[iY]) < (PixelHeight / 2))
		{
			Cy_Ar[iY] = 0.0; /* Main antenna */
		}

	}

	// Calculate CY Values ?? Can optimize
	double Cx_A[sizeof(double) * iXmax];
	for(iX = 0; iX < iXmax; iX++){
		Cx_A[iX] = CxMin + (iX * PixelWidth);
	}


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
		start_comp = MPI_Wtime();
		int row = 0;
		for(iY = startRow; iY <= endRow; iY++){
			offset = row * iXmax * 3;
			for(iX = 0; iX < iXmax; iX++)
				{
						/* initial value of orbit = critical point Z= 0 */
						Zx = 0.0;
						Zy = 0.0;
						Zx2 = Zx * Zx;
						Zy2 = Zy * Zy;
						/* */
						for(Iteration = 0; Iteration < IterationMax && ((Zx2 + Zy2) < ER2); Iteration++)
						{
							Zy = (2 * Zx * Zy) + Cy_Ar[iY];
							Zx = Zx2 - Zy2 + Cx_A[iX];
							Zx2 = Zx * Zx;
							Zy2 = Zy * Zy;
						};

						/* compute  pixel color (24 bit = 3 bytes) */
						if (Iteration == IterationMax)
						{
							// Point within the set. Mark it as black
							ppm[offset + iX * 3] = 0;
							ppm[offset + iX * 3 + 1] = 0;
							ppm[offset + iX * 3+ 2] = 0;

						}
						else 
						{
							// Point outside the set. Mark it as white
							double c = 3*log((double)Iteration)/log((double)(IterationMax) - 1.0);
							// printf("%f\n",c);
							if (c < 1)
							{
								ppm[offset + iX * 3] = 0;
								ppm[offset + iX * 3 + 1] = 0;
								ppm[offset + iX * 3 + 2] = 255*c;
								// printf("%c\n",rowStruct.row_color[iX * 3 + 2]);
							}
							else if (c < 2)
							{
								ppm[offset + iX * 3] = 0;
								ppm[offset + iX * 3 + 1] = 255*(c-1);
								ppm[offset + iX * 3 + 2] = 255;
							}
							else
							{
								ppm[offset + iX * 3] = 255*(c-2);
								ppm[offset + iX * 3 + 1] = 255;
								ppm[offset + iX * 3 + 2] = 255;
							}
					}
				}
				row++;
		}



		int count = 1;
		int startPoint;
		while (count <= numtasks - 1)
		{	
			startPoint = (iYmax / numtasks) * count;

			printf("Start ::::: %i -> remainder : %i\n", startPoint, remainder);
			MPI_Recv(&ppm[startPoint * iXmax * 3],((iYmax/numtasks) + remainder)* iXmax * 3, MPI_UNSIGNED_CHAR, count , 0, MPI_COMM_WORLD, &stat);
			count++;
		}

		end_comp = MPI_Wtime();
		printf("Mandelbrot computational process time: %f\n", end_comp - start_comp);
		printf("Completed Computing Mandelbrot Set.\n");

		/* write to the file */
		fwrite(ppm, 1, iYmax * iXmax * 3 , fp);
		fclose(fp);

		end = MPI_Wtime();
		printf("File: %s successfully closed.\n", filename);
		printf("Mandelbrot total process time: %f\n", end - start);
		// Get the clock current time again
		// Subtract end from start to get the CPU time used.

		MPI_Finalize();
		return 0;



	}else{
		int row = 0;
		for(iY = startRow; iY <= endRow; iY++){
			offset = row * iXmax * 3;
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
						Zy = (2 * Zx * Zy) + Cy_Ar[iY];
						Zx = Zx2 - Zy2 + Cx;
						Zx2 = Zx * Zx;
						Zy2 = Zy * Zy;
					};

					/* compute  pixel color (24 bit = 3 bytes) */
					if (Iteration == IterationMax)
					{
						// Point within the set. Mark it as black
						p_segment[offset + iX * 3] = 0;
						p_segment[offset + iX * 3 + 1] = 0;
						p_segment[offset + iX * 3 + 2] = 0;
					}
					else 
					{
						// Point outside the set. Mark it as white
						double c = 3*log((double)Iteration)/log((double)(IterationMax) - 1.0);
						// printf("%f\n",c);
						if (c < 1)
						{
							p_segment[offset + iX * 3] = 0;
							p_segment[offset+ iX * 3 + 1] = 0;
							p_segment[offset+ iX * 3 + 2] = 255*c;
							// printf("%c\n",rowStruct.row_color[iX * 3 + 2]);
						}
						else if (c < 2)
						{
							p_segment[offset + iX * 3] = 0;
							p_segment[offset + iX * 3 + 1] = 255*(c-1);
							p_segment[offset + iX * 3 + 2] = 255;
						}
						else
						{
							p_segment[offset + iX * 3] = 255*(c-2);
							p_segment[offset+ iX * 3 + 1] = 255;
							p_segment[offset + iX * 3 + 2] = 255;
						}
					}
			}
			row++;
		}
		MPI_Send(p_segment, ((iYmax / numtasks) +  remainder) * iXmax  * 3 , MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
		// printf("Rank %i -> Done\n", rank);
		MPI_Finalize();
		}

		free(p_segment);
		
		



	}
