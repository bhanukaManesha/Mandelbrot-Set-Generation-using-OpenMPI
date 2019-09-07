//////////////////////////////////////////////////////////////////////////////////////
// Bhanuka Manesha Samarasekera Vitharana Gamage
// 28993373
// bsam00002@student.monash.edu
// Parallel Computing - Assignment 1
// mandelbrot.c program: Mandelbort Set Fractal (Naive Row Based Partition Scheme Implementation).
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
#include <string.h>
#include <mpi.h>

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

	// Array to store the final image in memory for writing
	static unsigned char ppm[sizeOfppm];

	// Varaibles to store the start and end row for each task
	int startRow, endRow, remainder;

	// Varaible to store the off offset
	int offset;


	// Calculate the start and end values for each process
    startRow = (iYmax / numtasks) * rank;
    endRow = (iYmax / numtasks) * (rank + 1) - 1;

	// Calculate the remainder if the process cannot be equally divided
	remainder = iYmax % numtasks;

	// Allocate the remainder to the final process
    if (rank == numtasks -1 ){
        endRow = endRow + remainder;
    }

	// Initialze a dynamic array to store the rows calcuated by each task
	unsigned char* p_segment = (unsigned char*) malloc(((iYmax/numtasks) + remainder) * iXmax * 3);

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
		
		// Start the clock to calculate the time taken for generation of the Mandelbrot Set
		start_comp = MPI_Wtime();

		// Initialze variable to keep track of the row 
		int row = 0;
		// Loop though all the iY values until the end point for each process is reached
		for(iY = startRow; iY <= endRow; iY++){
			// Calculate the offset for each row
			offset = row * iXmax * 3;
			// Loop though all the iX values until iXmax
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
							Zx = Zx2 - Zy2 + Cx_Ar[iX];
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
				// Increment the row
				row++;
		}
		
		// Set the counter to 1 since rank 0 is done
		int count = 1;
		// Define a variable to store the start point
		int startPoint;

		// Loop until all the data is recieved
		while (count <= numtasks - 1)
		{	
			// Determine the start point
			startPoint = (iYmax / numtasks) * count;
			// Get the generated data from the slaves and set them straight into the final matrix
			MPI_Recv(&ppm[startPoint * iXmax * 3],((iYmax/numtasks) + remainder)* iXmax * 3, MPI_UNSIGNED_CHAR, count , 0, MPI_COMM_WORLD, &stat);
			// Increment the count to signify that one process is done
			count++;
		}

		// Get the time in order to calculate the total computation time
		end_comp = MPI_Wtime();

		// Print to console
		printf("Mandelbrot computational process time: %f\n", end_comp - start_comp);
		printf("Completed Computing Mandelbrot Set.\n");

		/* write to the file */
		fwrite(ppm, 1, iYmax * iXmax * 3 , fp);
		// Close the file
		fclose(fp);

		// Call MPI Finalize to end the MPI execution
		MPI_Finalize();

		// Get the time to calculate time for the parallel code
		end = MPI_Wtime();
		
		// Print to console
		printf("File: %s successfully closed.\n", filename);
		printf("Mandelbrot total process time: %f\n", end - start);
		return 0;
	}
	// Code segment to be executed by the slave nodes
	else
	{	
		// Initialize to varaible to keep track of the rows in the p_segment
		int row = 0;
		// Loop though all the iY values until the end point for each process is reached
		for(iY = startRow; iY <= endRow; iY++){
			// Calculate the offset for each row
			offset = row * iXmax * 3;
			// Loop though all the iX values until iXmax
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
						Zx = Zx2 - Zy2 + Cx_Ar[iX];
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
						
						if (c < 1)
						{
							p_segment[offset + iX * 3] = 0;
							p_segment[offset+ iX * 3 + 1] = 0;
							p_segment[offset+ iX * 3 + 2] = 255*c;
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
				// Increment the row count
				row++;
			}
		// Send the generated rows from each row back to the master node
		MPI_Send(p_segment, ((iYmax / numtasks) +  remainder) * iXmax  * 3 , MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
		// Call MPI Finalize to end the MPI execution
		MPI_Finalize();
		}

	// Finally free up the memory used by the dynamic array
	free(p_segment);
	return 0;
		



	}
