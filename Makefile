# Makefile for the surface programs

hvb=hvb/
filam=filam/

FLAGS=-pedantic -Wall -g -I$(hvb) -I$(filam) -O3 
LIBS=-llapack -lblas -lm -L/usr/X11R6/lib -lX11

OBJECTS= $(hvb)matrix.o $(hvb)common.o $(hvb)easyx.o \
	$(filam)geom.o $(filam)filam.o $(filam)dynamics.o $(filam)graphics.o

objects: $(OBJECTS)

%.o: %.cc
	g++ -c $< $(FLAGS) -o $@

x%: x%.o objects
	g++ $(FLAGS) -o $@ $< $(OBJECTS) $(LIBS)

clean:
	ls x* | grep -v "\." | xargs rm -f
	rm -f *.o 

mrproper: clean
	rm -f $(hvb)*.o $(filam)*.o

