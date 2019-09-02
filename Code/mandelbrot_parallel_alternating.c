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


#define CxMin -2.5
#define CxMax 1.5
#define CyMin -2.0
#define CyMax 2.0


// Main program
int main(int argc, char *argv[])
 {
	/* Clock information */
	double start, start_comp, end, end_comp;
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
	const int IterationMax = 2000; // default

	/* bail-out value , radius of circle ;  */
	const double EscapeRadius = 400;
	double ER2 = EscapeRadius * EscapeRadius;

	int balancedStop = iYmax - (iYmax % numtasks);
	int remainderStart = balancedStop + 1;

	unsigned char* p_segment = (unsigned char*) malloc((balancedStop/numtasks) * iXmax * sizeof(unsigned char) * 3);
	unsigned char* r_segment = (unsigned char*) malloc(numtasks * (balancedStop/numtasks) * iXmax * sizeof(unsigned char) * 3);

	// Calculate CY Values ?? Can optimize
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

	// Calculate CY Values ?? Can optimize
	double Cx_A[sizeof(double) * iXmax];
	for(iX = 0; iX < iXmax; iX++){
		Cx_A[iX] = CxMin + (iX * PixelWidth);
	}
		

	if (rank == 0){


		static unsigned char ppm[iXmax * 3 * iYmax];
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

		iY = rank;
		int offset;
		int j = 0;
		while (iY < iYmax){
			offset = iY * iXmax * 3;
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
			
			if (iY >= balancedStop){
				iY = remainderStart + j;
				j++;
			}else{
				iY+= numtasks;
				
			}
				// printf("Rank : %i -> iY : %i\n", rank, iY);
		}
		int count = 0;
		int r_segment_offset = (balancedStop/numtasks)* iXmax * 3;
		
		unsigned char** receive_order;
		receive_order = malloc(numtasks * sizeof(unsigned char*));

		while (count <= numtasks - 2)
		{
			MPI_Recv(r_segment + r_segment_offset * count,(balancedStop/numtasks)* iXmax * 3, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE , 0, MPI_COMM_WORLD, &stat);
			receive_order[stat.MPI_SOURCE] = r_segment + r_segment_offset * count;
			count++;
		}

		int rowNumber = 0;
		unsigned char* pointerTo;
		int process;

		while (rowNumber < balancedStop){
			if (rowNumber % numtasks != 0){
				// Get the process number which did the row
				process = rowNumber % numtasks;

				// Get the pointer to the last row
				pointerTo = receive_order[process];

				memcpy(&ppm[rowNumber * iXmax * 3], pointerTo, iXmax * 3);

				// Get the offset to the specific row number
				receive_order[process] = pointerTo + iXmax * 3;

			}
			rowNumber++;
		}	


		end_comp = MPI_Wtime();
		printf("Mandelbrot computational process time: %f\n", end_comp - start_comp);
		printf("Completed Computing Mandelbrot Set.\n");
	
		/* write color to the file */
		fwrite(ppm, 1, iYmax * iXmax * 3 , fp);
		fclose(fp);

		free(receive_order);
		end = MPI_Wtime();
		
		printf("File: %s successfully closed.\n", filename);
		printf("Mandelbrot total process time: %f\n", end - start);
		
		

	}else{
			// Get current clock time.
			int x = MPI_Wtime();
			iY = rank;
			int row = 0;
			int offset;
			while (iY < balancedStop){
				// printf("Rank : %i -> iY : %i\n", rank, iY);
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
				iY+= numtasks;
				
			}
			MPI_Send(p_segment, (balancedStop/numtasks) * iXmax  * 3, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
			
	}	

		MPI_Finalize();
		free(p_segment);
		free(r_segment);
		return 0;
		

 }