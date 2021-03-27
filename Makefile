CC = gcc
LD = gcc
CFLAGS = -O2 -g -Wall -ffast-math -ftree-vectorize -fopt-info-vec -fopenmp
LDFLAGS = 
LDLIBS = -lm -lX11 -fopenmp
RM = /bin/rm -f
OBJS = barnes_hut.o graphics/graphics.o
EXECUTABLE = nbody

$(EXECUTABLE): $(OBJS)
	$(LD) $(LDFLAGS) $^ $(LDLIBS) -o $@

barnes_hut.o: barnes_hut.c graphics/graphics.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

graphics/graphics.o: graphics/graphics.c graphics/graphics.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	$(RM) $(EXECUTABLE) $(OBJS) result.gal
