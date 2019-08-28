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
#include <string.h>

#define iXmax 8000 // default
#define iYmax 8000 // default
#define divisionFactor 0.45


#define CxMin -2.5
#define CxMax 1.5
#define CyMin -2.0
#define CyMax 2.0


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

	/* world ( double) coordinate = parameter plane*/
	double Cx, Cy;

	/* */
	double PixelWidth = (CxMax - CxMin)/iXmax;
	double PixelHeight = (CyMax - CyMin)/iYmax;

	/* color component ( R or G or B) is coded from 0 to 255 */
	/* it is 24 bit color RGB file */
	const int MaxColorComponentValue = 255; 

	/* Z = Zx + Zy*i;	Z0 = 0 */
	double Zx, Zy;
	double Zx2, Zy2; /* Zx2 = Zx*Zx;  Zy2 = Zy*Zy  */
	/*  */
	int Iteration;
	const int IterationMax = 1000; // default

	/* bail-out value , radius of circle ;  */
	const double EscapeRadius = 400;
	double ER2 = EscapeRadius * EscapeRadius;

	// Pack Buffer
	int position;
	const int row_packsize = sizeof(unsigned char) * iYmax * 3 + sizeof(int);
	char row_packbuf[row_packsize];

	int rowLength = row_packsize - sizeof(int);

	int row_i;
	unsigned char row_color[iXmax * 3];

	/* Clock information */
	double start, end;

	int lowerI = 0;
	int upperI = 0;
	if (numtasks > 2){
		upperI = (iYmax * divisionFactor);
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
		start = MPI_Wtime();

		static unsigned char ppm[iXmax * 3 * iYmax];
		int stopCount = 0;

		// Resetting all
		for(int i=1; i < numtasks; i++){
			MPI_Recv(row_packbuf, row_packsize, MPI_PACKED, i , 0, MPI_COMM_WORLD, &stat);
			if (i % 2 == 0) {
				MPI_Send(&lowerI, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			}else{
				MPI_Send(&upperI, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			}
		}

		while (1){
			position = 0;
			MPI_Recv(row_packbuf, row_packsize, MPI_PACKED, MPI_ANY_SOURCE , 0, MPI_COMM_WORLD, &stat);
 			MPI_Unpack(row_packbuf, row_packsize, &position, &row_i, 1, MPI_INT, MPI_COMM_WORLD);
        	MPI_Unpack(row_packbuf, row_packsize, &position, row_color, rowLength, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);
        	// printf("Row : %i -> Rank : %i\n", row_i,stat.MPI_SOURCE);
		
			memcpy(&ppm[row_i * iXmax * 3],row_color, iXmax * 3 );
			// memcpy(&ppm[(iXmax - row_i - 1) * iXmax * 3],row_color, iXmax * 3 );
				
			// Assigning task to the top half of the cores
			if (stat.MPI_SOURCE % 2 == 0) {
				lowerI += 1;
				if (lowerI >= iYmax * divisionFactor){
					if (upperI >= iYmax){
						// printf("Rank : %i\n", stat.MPI_SOURCE);
						stopCount += 1;
						// printf("Stop : %i\n", rank);
						lowerI = iYmax + 1 ;
						MPI_Send(&lowerI, 1, MPI_INT, stat.MPI_SOURCE, 0, MPI_COMM_WORLD);
						if (stopCount == numtasks - 1){
							break;
						}
						continue;
					}else{
						lowerI = upperI + 1;
					}
				}
				MPI_Send(&lowerI, 1, MPI_INT, stat.MPI_SOURCE, 0, MPI_COMM_WORLD);

			}else{
				upperI += 1;
				if (upperI >= iYmax ){
					
					stopCount += 1;
					// printf("Stop : %i\n", rank);
					upperI = iYmax + 1;
					MPI_Send(&upperI, 1, MPI_INT, stat.MPI_SOURCE, 0, MPI_COMM_WORLD);
					if (stopCount == numtasks - 1){
						break;
					}
					continue;
				}
				MPI_Send(&upperI, 1, MPI_INT, stat.MPI_SOURCE, 0, MPI_COMM_WORLD);

			}
		}
		
		/* write color to the file */
		fwrite(ppm, 1, iYmax * iXmax * 3 , fp);

		end = MPI_Wtime(); 

		fclose(fp);

		printf("Completed Computing Mandelbrot Set.\n");
		printf("File: %s successfully closed.\n", filename);
		printf( "Mandelbrot computational process time %f\n", end - start ); 
		// printf("Mandelbrot computational process time: %lf\n", cpu_time_used);
		MPI_Finalize();
		return 0;


	}else{	

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

		while (1){
			
			position = 0;
			MPI_Pack( &row_i, 1, MPI_INT, row_packbuf, row_packsize, &position, MPI_COMM_WORLD );
			MPI_Pack(row_color,  rowLength, MPI_UNSIGNED_CHAR, row_packbuf, row_packsize, &position, MPI_COMM_WORLD );
			MPI_Send(row_packbuf, position, MPI_PACKED, 0, 0, MPI_COMM_WORLD);

			MPI_Recv(&row_i, 1, MPI_INT, 0 , 0, MPI_COMM_WORLD, &stat);

			if (row_i > iYmax){
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
					Zy = (2 * Zx * Zy) + Cy_Ar[row_i];
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
					// printf("%f\n",c);
					if (c < 1)
					{
						row_color[iX * 3] = 0;
						row_color[iX * 3 + 1] = 0;
						row_color[iX * 3 + 2] = 255*c;
						// printf("%c\n",rowStruct.row_color[iX * 3 + 2]);
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

		}

		MPI_Finalize();

	}

 }