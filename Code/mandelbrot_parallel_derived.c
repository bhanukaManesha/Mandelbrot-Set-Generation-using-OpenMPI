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
	// static unsigned char color[3];	

	/* Z = Zx + Zy*i;	Z0 = 0 */
	double Zx, Zy;
	double Zx2, Zy2; /* Zx2 = Zx*Zx;  Zy2 = Zy*Zy  */
	/*  */
	int Iteration;
	const int IterationMax = 1000; // default

	/* bail-out value , radius of circle ;  */
	const double EscapeRadius = 400;
	double ER2 = EscapeRadius * EscapeRadius;


	// Custom Datatype Definition
	struct {
        int rowIndex;
        unsigned char row_color[iXmax * 3];
    } rowStruct;

	struct {
        int rowIndex;
        double Cy_Send;
    } CyStruct;

	int blocklengths_R[2];
	int blocklengths_C[2];
    MPI_Aint pos_R[2];
	MPI_Aint pos_S[2];

	MPI_Datatype RSTRUCT;
	MPI_Datatype CSTRUCT;
    MPI_Datatype old_types_R[2];
	MPI_Datatype old_types_S[2];

	blocklengths_R[0] = 1;
    blocklengths_R[1] = iXmax * 3;
	blocklengths_C[0] = 1;
    blocklengths_C[1] = 1;


    old_types_R[0] = MPI_INT;
    old_types_R[1] = MPI_UNSIGNED_CHAR;
	old_types_S[0] = MPI_INT;
    old_types_S[1] = MPI_DOUBLE;


    MPI_Get_address( &rowStruct.rowIndex, &pos_R[0] );
    MPI_Get_address( &rowStruct.row_color, &pos_R[1] );

    pos_R[1] = pos_R[1] - pos_R[0];
    pos_R[0] = 0;

	MPI_Get_address( &CyStruct.rowIndex, &pos_S[0] );
    MPI_Get_address( &CyStruct.Cy_Send, &pos_S[1] );

    pos_S[1] = pos_S[1] - pos_S[0];
    pos_S[0] = 0;

    MPI_Type_create_struct(2, blocklengths_R,pos_R,old_types_R,&RSTRUCT);
	MPI_Type_create_struct(2, blocklengths_C,pos_S,old_types_S,&CSTRUCT);

    MPI_Type_commit( &RSTRUCT );
	MPI_Type_commit( &CSTRUCT );

	// End Custom Datatype Definition

	
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
		const float DIVISION_FACTOR = 0.85;
		int lowerI = 0;
		int upperI = iYmax / 2 * DIVISION_FACTOR;
		int stopCount = 0;

		// printf("\n\n%lu\n",sizeof(Cy_Ar));

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
			MPI_Recv(&rowStruct, 1, RSTRUCT, i , 0, MPI_COMM_WORLD, &stat);
			
			if (i < numtasks / 2) {
				CyStruct.rowIndex = lowerI;
				CyStruct.Cy_Send = Cy_Ar[lowerI];
				// printf("LowerI ; %i\n",lowerI);
				lowerI += 1;

			}else{
				CyStruct.rowIndex = upperI;
				CyStruct.Cy_Send = Cy_Ar[upperI];
				// printf("UpperI ; %i\n",upperI);
				upperI += 1;
				
			}
			
			MPI_Send(&CyStruct, 1, CSTRUCT, i, 0, MPI_COMM_WORLD);
		}

		
		
		while (1){
			MPI_Recv(&rowStruct, 1, RSTRUCT, MPI_ANY_SOURCE , 0, MPI_COMM_WORLD, &stat);
			// printf("stopCount ; %i\n",stopCount);
			for (int t = 0;t < (iXmax) * 3;t++){ 

				// printf("Upper : %i \t Lower : %i \n", upperIndex, lowerIndex);
				ppm[ (rowStruct.rowIndex * iXmax * 3) + t ] = rowStruct.row_color[t];
				ppm[ ((iXmax - rowStruct.rowIndex - 1) * iXmax * 3) + t ] = rowStruct.row_color[t];
				
			}
			// Assigning task to the top half of the cores
			// printf("rank: %i\n", stat.MPI_SOURCE);
			if (stat.MPI_SOURCE < numtasks / 2) {
				CyStruct.rowIndex = lowerI;
				CyStruct.Cy_Send = Cy_Ar[lowerI];
				lowerI += 1;
				
				if (lowerI > (iYmax / 2) * DIVISION_FACTOR){
					// printf("LowerI ; %i\n",lowerI);
					stopCount += 1;
					CyStruct.rowIndex = iYmax + 1 ;
					MPI_Send(&CyStruct, 1, CSTRUCT, stat.MPI_SOURCE, 0, MPI_COMM_WORLD);
					if (stopCount == numtasks - 1){
						break;
					}
					continue;
				}
				MPI_Send(&CyStruct, 1, CSTRUCT, stat.MPI_SOURCE, 0, MPI_COMM_WORLD);

			}else{
				CyStruct.rowIndex = upperI;
				CyStruct.Cy_Send = Cy_Ar[upperI];
				upperI += 1;
				// printf("UpperI ; %i\n",upperI);
				// printf("rank: %i\n", stat.MPI_SOURCE);
				
				if (upperI > (iYmax / 2)){
					// printf("UpperI ; %i\n",upperI);
					stopCount += 1;
					CyStruct.rowIndex = iYmax + 1 ;
					MPI_Send(&CyStruct, 1, CSTRUCT, stat.MPI_SOURCE, 0, MPI_COMM_WORLD);
					if (stopCount == numtasks - 1){
						break;
					}
					continue;
				}
				MPI_Send(&CyStruct, 1, CSTRUCT, stat.MPI_SOURCE, 0, MPI_COMM_WORLD);

			}

			

		}
		
		/* write color to the file */
		fwrite(ppm, 1, iYmax * iXmax * 3 , fp);

		
		// Get the clock current time again
		// Subtract end from start to get the CPU time used.
		// end = clock();
		// cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

		end = MPI_Wtime(); 

		fclose(fp);

		MPI_Type_free( &CSTRUCT );
		MPI_Type_free( &RSTRUCT );


		printf("Completed Computing Mandelbrot Set.\n");
		printf("File: %s successfully closed.\n", filename);
		printf( "Mandelbrot computational process time %f\n", end - start ); 
		// printf("Mandelbrot computational process time: %lf\n", cpu_time_used);
		MPI_Finalize();
		return 0;


	}else{
		rowStruct.rowIndex = 0;
		
		while (1){
			
			// printf("\n\n%lu\n",sizeof(CyStruct));

			MPI_Send(&rowStruct,1, RSTRUCT, 0, 0, MPI_COMM_WORLD);
			
			MPI_Recv(&CyStruct,1, CSTRUCT, 0, 0, MPI_COMM_WORLD, &stat);
			
			rowStruct.rowIndex = CyStruct.rowIndex;

			if (CyStruct.rowIndex > iYmax / 2){
				
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
					Zy = (2 * Zx * Zy) + CyStruct.Cy_Send;
					Zx = Zx2 - Zy2 + Cx;
					Zx2 = Zx * Zx;
					Zy2 = Zy * Zy;
				};

				/* compute  pixel color (24 bit = 3 bytes) */
				if (Iteration == IterationMax)
				{
					// Point within the set. Mark it as black
					rowStruct.row_color[iX * 3] = 0;
					rowStruct.row_color[iX * 3 + 1] = 0;
					rowStruct.row_color[iX * 3 + 2] = 0;
				}
				else 
				{
					// Point outside the set. Mark it as white
					double c = 3*log((double)Iteration)/log((double)(IterationMax) - 1.0);
					// printf("%f\n",c);
					if (c < 1)
					{
						rowStruct.row_color[iX * 3] = 0;
						rowStruct.row_color[iX * 3 + 1] = 0;
						rowStruct.row_color[iX * 3 + 2] = 255*c;
						// printf("%c\n",rowStruct.row_color[iX * 3 + 2]);
					}
					else if (c < 2)
					{
						rowStruct.row_color[iX * 3] = 0;
						rowStruct.row_color[iX * 3 + 1] = 255*(c-1);
						rowStruct.row_color[iX * 3 + 2] = 255;
					}
					else
					{
						rowStruct.row_color[iX * 3] = 255*(c-2);
						rowStruct.row_color[iX * 3 + 1] = 255;
						rowStruct.row_color[iX * 3 + 2] = 255;
					}
				}

			}

		}

		MPI_Finalize();
		
		



	}

 }