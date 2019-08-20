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
	const int c_packsize = sizeof(double) + sizeof(int);
	const int row_packsize = sizeof(unsigned char) * iYmax * 3 + sizeof(int);
    char c_packbuf[c_packsize];
	char row_packbuf[row_packsize];

	int rowLength = sizeof(unsigned char) * iYmax * 3;

	int row_i;
	unsigned char row_color[iXmax * 3];

	
	/* Clock information */
	double start, end;

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
		double Cy_Ar[sizeof(double) * iYmax];
		
		const float DIVISION_FACTOR = 0.83;
		int lowerI = 0;
		int upperI = (iYmax / 2 * DIVISION_FACTOR);
		int stopCount = 0;

		/* compute and write image data bytes to the file */
		for(iY = 0; iY < iYmax / 2; iY++)
		{	
			// printf("Started Row %i.\n",iY);
			Cy_Ar[iY] = CyMin + (iY * PixelHeight);
			if (fabs(Cy_Ar[iY]) < (PixelHeight / 2))
			{
				Cy_Ar[iY] = 0.0; /* Main antenna */
			}

		}
		
		// Resetting all
		for(int i=1; i<numtasks; i++){
			position = 0;
			MPI_Recv(row_packbuf, row_packsize, MPI_PACKED, i , 0, MPI_COMM_WORLD, &stat);
			if (i < numtasks / 2) {
				MPI_Pack( &lowerI, 1, MPI_INT, c_packbuf, c_packsize, &position, MPI_COMM_WORLD );
				MPI_Pack( &Cy_Ar[lowerI], 1 , MPI_DOUBLE, c_packbuf, c_packsize, &position, MPI_COMM_WORLD );
			}else{
				MPI_Pack( &upperI, 1, MPI_INT, c_packbuf, c_packsize, &position, MPI_COMM_WORLD );
				MPI_Pack( &Cy_Ar[upperI],  1, MPI_DOUBLE, c_packbuf, c_packsize, &position, MPI_COMM_WORLD );
			}
			
			MPI_Send(c_packbuf, position, MPI_PACKED, i, 0, MPI_COMM_WORLD);
		}

		while (1){
			position = 0;
			MPI_Recv(row_packbuf, row_packsize, MPI_PACKED, MPI_ANY_SOURCE , 0, MPI_COMM_WORLD, &stat);
 			MPI_Unpack(row_packbuf, row_packsize, &position, &row_i, 1, MPI_INT, MPI_COMM_WORLD);
        	MPI_Unpack(row_packbuf, row_packsize, &position, row_color, rowLength, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);
        
		
			for (int t = 0;t < (iXmax) * 3;t++){ 

				ppm[ (row_i * iXmax * 3) + t ] = row_color[t];
				ppm[ ((iXmax - row_i - 1) * iXmax * 3) + t ] = row_color[t];
				
			}
			// Assigning task to the top half of the cores
			position = 0;
			if (stat.MPI_SOURCE < numtasks / 2) {
				lowerI += 1;
				row_i = lowerI;
				if (lowerI >= (iYmax / 2) * DIVISION_FACTOR){
					stopCount += 1;
					row_i = iYmax + 100000 ;
					MPI_Pack( &row_i, 1, MPI_INT, c_packbuf, c_packsize, &position, MPI_COMM_WORLD );
					MPI_Send(c_packbuf, position, MPI_PACKED, stat.MPI_SOURCE, 0, MPI_COMM_WORLD);

					if (stopCount == numtasks - 1){
						break;
					}
					continue;
				}
				MPI_Pack( &row_i, 1, MPI_INT, c_packbuf, c_packsize, &position, MPI_COMM_WORLD );
				MPI_Pack( &Cy_Ar[lowerI], 1 , MPI_DOUBLE, c_packbuf, row_packsize, &position, MPI_COMM_WORLD );
				MPI_Send(c_packbuf, position, MPI_PACKED, stat.MPI_SOURCE, 0, MPI_COMM_WORLD);

			}else{
				upperI += 1;
				row_i = upperI;
				
				if (upperI >= (iYmax / 2 )){
					stopCount += 1;
					row_i = iYmax + 10000;
					MPI_Pack( &row_i, 1, MPI_INT, c_packbuf, c_packsize, &position, MPI_COMM_WORLD );
					MPI_Send(c_packbuf, position, MPI_PACKED, stat.MPI_SOURCE, 0, MPI_COMM_WORLD);
					if (stopCount == numtasks - 1){
						break;
					}
					continue;
				}
				MPI_Pack( &upperI, 1, MPI_INT, c_packbuf, c_packsize, &position, MPI_COMM_WORLD );
				MPI_Pack( &Cy_Ar[upperI],  1, MPI_DOUBLE, c_packbuf, row_packsize, &position, MPI_COMM_WORLD );
				MPI_Send(c_packbuf, position, MPI_PACKED, stat.MPI_SOURCE, 0, MPI_COMM_WORLD);

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
		
		double c;
		
		while (1){
				
			position = 0;
			MPI_Pack( &row_i, 1, MPI_INT, row_packbuf, row_packsize, &position, MPI_COMM_WORLD );
			MPI_Pack(row_color,  rowLength, MPI_UNSIGNED_CHAR, row_packbuf, row_packsize, &position, MPI_COMM_WORLD );
			MPI_Send(row_packbuf, position, MPI_PACKED, 0, 0, MPI_COMM_WORLD);

			position = 0;			
			MPI_Recv(c_packbuf, c_packsize, MPI_PACKED, 0 , 0, MPI_COMM_WORLD, &stat);
 			MPI_Unpack(c_packbuf, c_packsize, &position, &row_i, 1, MPI_INT, MPI_COMM_WORLD);
        	MPI_Unpack(c_packbuf, c_packsize, &position, &c, 1, MPI_DOUBLE, MPI_COMM_WORLD);

			if (row_i > iYmax / 2){
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
					Zy = (2 * Zx * Zy) + c;
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