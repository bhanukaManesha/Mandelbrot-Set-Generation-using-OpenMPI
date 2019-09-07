ALL: mandelbrot_serial mandelbrot_parallel_naive mandelbrot_parallel_fcfs mandelbrot_parallel_alternating

mandelbrot_serial : mandelbrot_serial.c
	gcc mandelbrot_serial.c -o mandelbrot_serial -lm

mandelbrot_parallel_naive : mandelbrot_parallel_naive.c
	mpicc mandelbrot_parallel_naive.c -o mandelbrot_parallel_naive  -lm

mandelbrot_parallel_fcfs : mandelbrot_parallel_fcfs.c
	mpicc mandelbrot_parallel_fcfs.c -o mandelbrot_parallel_fcfs  -lm

mandelbrot_parallel_alternating : mandelbrot_parallel_alternating.c
	mpicc mandelbrot_parallel_alternating.c -o mandelbrot_parallel_alternating  -lm

run:
	./mandelbrot_serial
	mpirun -np 4 mandelbrot_parallel_naive
	mpirun -np 4 mandelbrot_parallel_fcfs
	mpirun -np 4 mandelbrot_parallel_alternating

clean :
	/bin/rm -f mandelbrot_serial *.o
	/bin/rm -f mandelbrot_parallel_naive *.o
	/bin/rm -f mandelbrot_parallel_fcfs *.o
	/bin/rm -f mandelbrot_parallel_alternating *.o
