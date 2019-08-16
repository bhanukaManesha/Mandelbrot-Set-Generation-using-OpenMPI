ALL: mandelbrot_serial mandelbrot_parallel

mandelbrot_serial : mandelbrot_serial.c
	mpicc -o mandelbrot_serial mandelbrot_serial.c

mandelbrot_parallel : mandelbrot_parallel.c
	mpicc -o mandelbrot_parallel mandelbrot_parallel.c

run:
	mpirun -np 4 mandelbrot_parallel

clean :
	/bin/rm -f mandelbrot_serial *.o
	/bin/rm -f mandelbrot_parallel *.o
