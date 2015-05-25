all:
	gcc -Wall -I/usr/include -I./include  testfile.c -lgsl -lgslcblas -lm
clean:
	rm run
