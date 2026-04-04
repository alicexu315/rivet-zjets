g++ PythiaGen.cc -o PythiaGen `pythia8-config --cppflags --libs` -L${HEPMC2_DIR}/lib -Wl,-rpath,${HEPMC2_DIR}/lib -lHepMC
