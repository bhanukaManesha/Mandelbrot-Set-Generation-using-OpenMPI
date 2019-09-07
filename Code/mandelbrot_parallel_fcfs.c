//////////////////////////////////////////////////////////////////////////////////////
// Bhanuka Manesha Samarasekera Vitharana Gamage
// 28993373
// bsam00002@student.monash.edu
// Parallel Computing - Assignment 1
// mandelbrot.c program: Mandelbort Set Fractal (First Come First Served Row Based Partition Scheme Implementation).
//
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

// Including the libraries for the preprocessor
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <string.h>

// Define the macro for the iXmax
#define iXmax 8000 // default
// Define the macro for the iYmax
#define iYmax 8000 // default
// Define the macro for the size of the final mandelbrot image
#define sizeOfppm iXmax * iYmax * 3


// Main program
int main(int argc, char *argv[])
 {
	// Defining the variables used for calculating the time
	double start, start_comp, end, end_comp;

	// Start the clock for calculating the total time
	start = MPI_Wtime();

	// Defining the variables to store the number of tasks and rank
	int numtasks, rank;

	// Initialize MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Status stat;

	// If core count less than 1, then abort
	if (numtasks < 2){
		printf("The task count should be more than 1\n");
		MPI_Finalize();
		return 0;
	}

	/* screen ( integer) coordinate */
	int iX,iY;

	/* world ( double) coordinate = parameter plane*/
	double Cx, Cy;
	const double CxMin = -2.5;
	const double CxMax = 1.5;
	const double CyMin = -2.0;
	const double CyMax = 2.0;

	/* Store the pixel width and height */
	double PixelWidth = (CxMax - CxMin)/iXmax;
	double PixelHeight = (CyMax - CyMin)/iYmax;

	/* color component ( R or G or B) is coded from 0 to 255 */
	/* it is 24 bit color RGB file */
	const int MaxColorComponentValue = 255; 

	// Define the variables for the Zx, Zy, Zx2 and Zy2 values
	double Zx, Zy;
	double Zx2, Zy2;
	
	// Define the variables to keep track of the Iteration and Iteration Max
	int Iteration;
	const int IterationMax = 2000; // default

	/* bail-out value , radius of circle ;  */
	const double EscapeRadius = 400;
	double ER2 = EscapeRadius * EscapeRadius;

	// Pointer to be used by the MPI_PACK
	int position;

	// Calculating the size of each row with the row indetifier
	const int row_packsize = sizeof(unsigned char) * iYmax * 3 + sizeof(int);

	// Pack Buffer to pack and send each row with the identifier
	char row_packbuf[row_packsize];

	// Calculating the length of each row without identifier
	int rowLength = row_packsize - sizeof(int);

	// Interger to keep track of the row
	int row_i;

	// Buffer to store each row without the identifier
	unsigned char row_color[iXmax * 3];

	// Integer that keeps track of the completed rows
	int rowCount = 0;

	// Code segment to be executed by the master node
	if (rank == 0){

		// Define the pointer to the file
		FILE * fp;
		char *filename = "Mandelbrot.ppm";
		char *comment = "# ";	/* comment should start with # */

		/*create new file,give it a name and open it in binary mode  */
		fp = fopen(filename, "wb"); /* b -  binary mode */

		/*write ASCII header to the file (PPM file format)*/
		fprintf(fp,"P6\n %s\n %d\n %d\n %d\n", comment, iXmax, iYmax, MaxColorComponentValue);
		
		// Print to console
		printf("File: %s successfully opened for writing.\n", filename);
		printf("Computing Mandelbrot Set. Please wait...\n");

		// Array to store the final image in memory for writing
		static unsigned char ppm[sizeOfppm];

		// Integer to keep track of the tasks completed so far
		int stopCount = 0;

		// Get current clock time.
		start_comp = MPI_Wtime();

		// Loop to get the initial message from the slaves
		for(int i=1; i < numtasks; i++){
			// Recieve the junk variable
			MPI_Recv(row_packbuf, row_packsize, MPI_PACKED, i , 0, MPI_COMM_WORLD, &stat);
			// Send the row needed to be done by the task
			MPI_Send(&rowCount, 1, MPI_INT, i, 0, MPI_COMM_WORLD);	
		}

		// Run the loop until all the rows are computed
		while (1){

			// Reset the buffer position to 0
			position = 0;

			// Receive the completed task from the slave
			MPI_Recv(row_packbuf, row_packsize, MPI_PACKED, MPI_ANY_SOURCE , 0, MPI_COMM_WORLD, &stat);
 			// Unpack the row ID to get the row that was completed
			MPI_Unpack(row_packbuf, row_packsize, &position, &row_i, 1, MPI_INT, MPI_COMM_WORLD);
        	// Unpack the row
			MPI_Unpack(row_packbuf, row_packsize, &position, row_color, rowLength, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);
        		
			// Copy the row to the memory at the correct location using the row id
			memcpy(&ppm[row_i * iXmax * 3],row_color, iXmax * 3 );
				
			// Assigning task to the top half of the cores
			rowCount += 1;

			// If the row count more than the iYmax, then all the rows are processed
			if (rowCount >= iYmax){
					// Increase the stop count
					stopCount += 1;
					// Add the row number greater than iYmax
					rowCount = iYmax + 10 ;
					// Send the row number so the task can abort
					MPI_Send(&rowCount, 1, MPI_INT, stat.MPI_SOURCE, 0, MPI_COMM_WORLD);

					// Check if all the tasks are completed
					if (stopCount == numtasks - 1){
						// If yes, then break out of the loop
						break;
					}
					// Else, continue until everyone is done
					continue;
			}
			// Send the row number to the next task
			MPI_Send(&rowCount, 1, MPI_INT, stat.MPI_SOURCE, 0, MPI_COMM_WORLD);

		}

		// Get the time in order to calculate the total computation time
		end_comp = MPI_Wtime(); 

		// Print to console
		printf( "Mandelbrot computational process time %f\n", end_comp - start_comp ); 
		printf("Completed Computing Mandelbrot Set.\n");
		
		/* write color to the file */
		fwrite(ppm, 1, iYmax * iXmax * 3 , fp);

		// Close the file
		fclose(fp);

		// Call MPI Finalize to end the MPI execution
		MPI_Finalize();

		// Get the time to calculate time for the parallel code
		end = MPI_Wtime(); 

		// Print to console
		printf("File: %s successfully closed.\n", filename);
		printf( "Mandelbrot total process time %f\n", end - start ); 

		return 0;


	}else{	

		// Calculate CY Values
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

		// Calculate Cx Values
		double Cx_Ar[sizeof(double) * iXmax];
		for(iX = 0; iX < iXmax; iX++){
			Cx_Ar[iX] = CxMin + (iX * PixelWidth);
		}

		// Unitl no more reamining tasks
		while (1){
			
			// Set the position of MPI_PACK to 0
			position = 0;

			// Pack the row number to be sent at the front
			MPI_Pack( &row_i, 1, MPI_INT, row_packbuf, row_packsize, &position, MPI_COMM_WORLD );
			// Pack the row with all the color pixels
			MPI_Pack(row_color,  rowLength, MPI_UNSIGNED_CHAR, row_packbuf, row_packsize, &position, MPI_COMM_WORLD );
			// Send the row buffer to the master
			MPI_Send(row_packbuf, position, MPI_PACKED, 0, 0, MPI_COMM_WORLD);

			// Recieve the next row from the master
			MPI_Recv(&row_i, 1, MPI_INT, 0 , 0, MPI_COMM_WORLD, &stat);

			// If the row number is greater than iYmax
			if (row_i > iYmax){
				// Break from the loop
				break;
			}
			
			// Else, calculate the row
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
					Zy = (2 * Zx * Zy) + Cy_Ar[row_i];
					Zx = Zx2 - Zy2 + Cx_Ar[iX];
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

		}
		// Call MPI Finalize to end the MPI execution
		MPI_Finalize();
		return 0;
	}

 }