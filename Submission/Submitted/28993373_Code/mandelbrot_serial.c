//////////////////////////////////////////////////////////////////////////////////////
// Bhanuka Manesha Samarasekera Vitharana Gamage
// 28993373
// bsam00002@student.monash.edu
// Parallel Computing - Assignment 1
// mandelbrot.c program: Mandelbort Set Fractal (Color Serial Code Implementation).
//
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

// Define the macro for the iXmax
#define iXmax 8000 // default
// Define the macro for the iYmax
#define iYmax 8000 // default
// Define the macro for the size of the final mandelbrot image
#define sizeOfppm iXmax * iYmax * 3

// Main program
int main()
 {
	// Defining the variables used for calculating the time
	clock_t startAll, endAll, startComp, endComp;
	// Defining the variables used for to store the calculated time
	double cpu_time_used_all, cpu_time_used_comp;

	// Get current clock time (to measure the overall program time).
	startAll = clock();

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
	FILE * fp;
	char *filename = "Mandelbrot.ppm";
	char *comment = "# ";	/* comment should start with # */

	// RGB color array
	static unsigned char color[3];

	// Define the variables for the Zx, Zy, Zx2 and Zy2 values
	double Zx, Zy;
	double Zx2, Zy2;
	
	// Define the variables to keep track of the Iteration and Iteration Max
	int Iteration;
	const int IterationMax = 2000; // default

	/* bail-out value , radius of circle ;  */
	const double EscapeRadius = 400;
	double ER2 = EscapeRadius * EscapeRadius;

	/*create new file,give it a name and open it in binary mode  */
	fp = fopen(filename, "wb"); /* b -  binary mode */

	/*write ASCII header to the file (PPM file format)*/
	fprintf(fp,"P6\n %s\n %d\n %d\n %d\n", comment, iXmax, iYmax, MaxColorComponentValue);

	printf("File: %s successfully opened for writing.\n", filename);
	printf("Computing Mandelbrot Set. Please wait...\n");

	// Define the ppm as an array
	static unsigned char ppm[sizeOfppm];

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

	// Get current clock time (to measure the multiplication time only).
	startComp = clock();

	/* compute and write image data bytes to the file */
	for(iY = 0; iY < iYmax; iY++)
	{
		Cy = CyMin + (iY * PixelHeight);
		if (fabs(Cy) < (PixelHeight / 2))
		{
			Cy = 0.0; /* Main antenna */
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
				Zy = (2 * Zx * Zy) + Cy_Ar[iY];
				Zx = Zx2 - Zy2 + Cx_Ar[iX];
				Zx2 = Zx * Zx;
				Zy2 = Zy * Zy;
			};

			/* compute  pixel color (24 bit = 3 bytes) */
			if (Iteration == IterationMax)
			{
				// Point within the set. Mark it as black
				color[0] = 0;
				color[1] = 0;
				color[2] = 0;
			}
			else 
			{
				// Point outside the set. Mark it as white
				double c = 3*log((double)Iteration)/log((double)(IterationMax) - 1.0);
				if (c < 1)
				{
					color[0] = 0;
					color[1] = 0;
					color[2] = 255*c;
				}
				else if (c < 2)
				{
					color[0] = 0;
					color[1] = 255*(c-1);
					color[2] = 255;
				}
				else
				{
					color[0] = 255*(c-2);
					color[1] = 255;
					color[2] = 255;
				}
			}

			// Copy the color into memory
			memcpy(&ppm[iY * iYmax * 3 + iX * 3],color, 3);
						
		}
	}
	
	// Get the clock current time again
	endComp = clock();
	
	// Subtract end from start to get the CPU time used.
	cpu_time_used_comp = ((double)(endComp - startComp)) / CLOCKS_PER_SEC;
	
	// print to console
	printf("Mandelbrot computational process time: %lf\n", cpu_time_used_comp);

	/* write color to the file */
	fwrite(ppm, 1, iXmax * iYmax * 3, fp);

	// Close the file
	fclose(fp);

	// Print to console
	printf("Completed Computing Mandelbrot Set.\n");
	printf("File: %s successfully closed.\n", filename);

	// Get the clock current time again
	endAll = clock();
	// Subtract end from start to get the CPU time used.
	cpu_time_used_all = ((double)(endAll - startAll)) / CLOCKS_PER_SEC;

	// Print to console
	printf("Mandelbrot Overall time: %lf\n", cpu_time_used_all);

	// Return 0
	return 0;
 }