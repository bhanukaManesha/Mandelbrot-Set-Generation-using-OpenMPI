ALL: mandelbrot_serial mandelbrot_parallel

mandelbrot_serial : mandelbrot_serial.c
	gcc mandelbrot_serial.c -o mandelbrot_serial -lm

mandelbrot_parallel : mandelbrot_parallel.c
	mpicc mandelbrot_parallel.c -o mandelbrot_parallel  -lm

run:
	mpirun -np 8 mandelbrot_parallel

clean :
	/bin/rm -f mandelbrot_serial *.o
	/bin/rm -f mandelbrot_parallel *.o
