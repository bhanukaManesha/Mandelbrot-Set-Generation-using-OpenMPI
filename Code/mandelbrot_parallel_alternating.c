//////////////////////////////////////////////////////////////////////////////////////
// Bhanuka Manesha Samarasekera Vitharana Gamage
// 28993373
// bsam00002@student.monash.edu
// Parallel Computing - Assignment 1
// mandelbrot.c program: Mandelbort Set Fractal (Alternating Row Based Partition Scheme Implementation).
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

	// Define the variable to calculate the last row number that divides
	// the image in to evenly distributed sections
	int balancedStop = iYmax - (iYmax % numtasks);
	// Get the row to start of the remiander section
	int remainderStart = balancedStop + 1;

	// Define a dynamic array to store the equak number of tasks
	unsigned char* p_segment = (unsigned char*) malloc((balancedStop/numtasks) * iXmax * sizeof(unsigned char) * 3);
	// Define an array to store all the values from the nodes before assigning to the main ppm array
	unsigned char* r_segment = (unsigned char*) malloc(numtasks * (balancedStop/numtasks) * iXmax * sizeof(unsigned char) * 3);

	// Calculate CY Values
	double Cy_Ar[sizeof(double) * iYmax];
	/* compute and write image data bytes to the file */
	for(iY = rank; iY < iYmax; iY++)
	{	
		Cy_Ar[iY] = CyMin + (iY * PixelHeight);
		if (fabs(Cy_Ar[iY]) < (PixelHeight / 2))
		{
			Cy_Ar[iY] = 0.0; /* Main antenna */
		}

	}

	// Calculate CY Values
	double Cx_Ar[sizeof(double) * iXmax];
	for(iX = 0; iX < iXmax; iX++){
		Cx_Ar[iX] = CxMin + (iX * PixelWidth);
	}
		
	// Code segment to be executed by the master node
	if (rank == 0){

		// Define the pointer to the file
		static unsigned char ppm[sizeOfppm];
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

		// Set the iY value to the rank, in order to start from row 0
		iY = rank;

		// Initialize the variable to calculate the offset
		int offset;
		
		// Set the j to zero, to be used in the case of remainder
		int j = 0;

		// Loop from start till iYmax
		while (iY < iYmax){
			// Calculate the offset
			offset = iY * iXmax * 3;

			// Loop though each value in th row
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
			// If rank zero has completed all the rows allocated for it
			if (iY >= balancedStop){
				// Then add the iY value for the remainder part
				iY = remainderStart + j;
				// Increment the j
				j++;
			}else{
				// if rank zero is not done with the allocated part, then increase by the number of tasks
				iY+= numtasks;
				
			}
		}

		// Variable to keep track of the completed number of tasks
		int count = 0;

		// Calculate the offset for each segment
		int r_segment_offset = (balancedStop/numtasks)* iXmax * 3;
		
		// Define an aray of pointers to to store the pointers to each of the segments in the r_segment array
		unsigned char** receive_order;

		// Generate a dynamic array to store the pointers to each of the segments in the r_segment array
		receive_order = malloc(numtasks * sizeof(unsigned char*));

		// Loop untill all the tasks sends back the data
		while (count <= numtasks - 2)
		{
			// Get all the rows from the slave nodes
			MPI_Recv(r_segment + r_segment_offset * count,(balancedStop/numtasks)* iXmax * 3, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE , 0, MPI_COMM_WORLD, &stat);
			// Set the pointer to the recive_order[rank] index, in order to quickly access the 
			// segment recieved by each node
			receive_order[stat.MPI_SOURCE] = r_segment + r_segment_offset * count;
			// Increase the counter by 1
			count++;
		}

		// Set the row number as 0
		int rowNumber = 0;

		// Create a char pointer to point to the row to be copied
		unsigned char* pointerTo;

		// Write the pointer number to store the process
		int process;

		// Loop though the ppm image until the remainder
		while (rowNumber < balancedStop){
			// If the row number is not done by rank 0
			if (rowNumber % numtasks != 0){
				// Get the process number which did the row
				process = rowNumber % numtasks;

				// Get the pointer to the last row
				pointerTo = receive_order[process];

				// Copy the value from the r_segment array using the pointer to the ppm array
				memcpy(&ppm[rowNumber * iXmax * 3], pointerTo, iXmax * 3);

				// Get the offset to the specific row number, and point the next pointer 
				// to the next row
				receive_order[process] = pointerTo + iXmax * 3;

			}
			// Increment the row number
			rowNumber++;
		}	

		// Get the time in order to calculate the total computation time
		end_comp = MPI_Wtime();

		// Print to console
		printf("Mandelbrot computational process time: %f\n", end_comp - start_comp);
		printf("Completed Computing Mandelbrot Set.\n");
	
		/* write color to the file */
		fwrite(ppm, 1, iYmax * iXmax * 3 , fp);
		// Close the file
		fclose(fp);

		// Free the memory used by recieve order pointer array
		free(receive_order);

		// Get the time to calculate time for the parallel code
		end = MPI_Wtime();
		
		// Print to console
		printf("File: %s successfully closed.\n", filename);
		printf("Mandelbrot total process time: %f\n", end - start);

	}else{
			// Get current clock time.
			int x = MPI_Wtime();
			// Set the start point as the rank
			iY = rank;
			// Set the row number to 0
			int row = 0;
			// Initialize variable to store the offset
			int offset;

			// Loop unitl the balanced stop is reached
			while (iY < balancedStop){
				// Calculate the offset for each row
				offset = row * iXmax * 3;

				// Loop though all the items in the row
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
				// Increase the row count
				row++;
				// Increment by the number of tasks
				iY+= numtasks;
				
			}
			// Send all the calculated rows
			MPI_Send(p_segment, (balancedStop/numtasks) * iXmax  * 3, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
			
	}	
		// Call MPI Finalize to end the MPI execution
		MPI_Finalize();

		// Finally free up the memory used by the dynamic array
		free(p_segment);
		free(r_segment);

		return 0;
		

 }