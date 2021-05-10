######## Please set this valuable BY YOURSELF ##############
VIENNA=/home/terai/local_tmp/
############################################################

all: compile link

compile:
	cd src; g++ -c -O3 -Wall CDSfold.cpp -I$(VIENNA)/include/ViennaRNA/ -I$(VIENNA)/include/ 
#	cd src; g++ -c -O3 -Wall CDSfold.cpp -I$(VIENNA)/H -I$(VIENNA)/include/ViennaRNA 

link:
	cd src; g++ -Wall -o CDSfold CDSfold.o -fmessage-length=0 -L$(VIENNA)/lib -fopenmp -lRNA

clean:
	-rm src/CDSfold.o
	-rm src/CDSfold
