coalseamsrc = $(filter-out src/pedsim.c, $(wildcard src/*.c))
pedsimsrc = $(filter-out src/coalseam.c, $(wildcard src/*.c))
headerfiles =  $(wildcard src/*.h)

all: coalseam

coalseam: $(coalseamsrc) $(headerfiles)
	gcc $(coalseamsrc) -DNDEBUG -lm -O3 -march=native -o coalseam

pedsim: $(pedsimsrc) $(headerfiles)
	gcc $(pedsimsrc) -DNDEBUG -lm -O3 -march=native -o pedsim 
