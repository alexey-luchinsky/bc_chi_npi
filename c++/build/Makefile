CC=g++ -std=c++14 
INC=-I ~/Work/EvtGen/evtgen/ -I`root-config --incdir`
SRCS= ../src/EvtProbeModel.cpp ../src/segR.cc
LIBS=~/Work/EvtGen/evtgen/build/lib/libEvtGen.a /usr/local/Cellar/hepmc3/3.2.6/lib/libHepMC3.dylib `root-config --libs`

all: segR_test.exe bc_chiJ.exe


segR_test.exe:
	$(CC) ../src/segR_test.cc ../src/EvtProbeModel.cpp ../src/segR.cc $(INC) $(LIBS) -o segR_test.exe

bc_chiJ.exe: ../src/bc_chiJ.cc ../src/EvtBcVHad.cpp ../src/EvtWHad.cpp ../src/segR.cc
	$(CC) ../src/bc_chiJ.cc ../src/EvtBcVHad.cpp ../src/EvtWHad.cpp ../src/segR.cc $(INC) $(LIBS) -o bc_chiJ.exe

clean:
	rm *.exe

#  

# g++ ../src/segR.cc ../src/segR_test.cc ../src/EvtProbeModel.cpp -I ~/Work/EvtGen/evtgen/ -I`root-config --incdir` -std=c++14 ~/Work/EvtGen/evtgen/build/lib/libEvtGen.a /usr/local/Cellar/hepmc3/3.2.6/lib/libHepMC3.dylib `root-config --libs`
