#Please change the following three lines to your own folders
GCC=/usr/bin/g++
INCLUDES= -I/usr/local/include/ -I$(BALL)/include/  -I$(BALL)/lcs
LIBS= -L /usr/local/lib/ -L $(LD_LIBRARY_PATH)

lcs: lcs.o mySurface.o
	$(GCC) $(LIBS) $(INCLUDES)  mySurface.o lcs.o -o lcs -lBALL -lnsl  -lm 

mySurface.o: mySurface.h mySurface.C
	$(GCC) $(INCLUDES) -c -g mySurface.h mySurface.C

lcs.o : lcs.C
	$(GCC) $(INCLUDES) -c -g lcs.C


.PHONY: clean
clean:
	-rm -f *.o
