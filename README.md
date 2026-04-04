# Rivet routine for the B-jet mass analysis

First, compile the Rivet routine into a `RivetAnalysis.so` file:
```
$ rivet-build LHCB_2019_I1730448.cc
```

Similarly, compile the provided `PythiaGen` script:
```
$ ./compile_PythiaGen_cc.cmd
```
and modify `PythiaGen.cmnd` according to the generation parameters you would like to use.

You can then use `PythiaGen` to create a .hepmc file or (preferred) generate events into a FIFO, which can then be read by Rivet:
```
$ mkfifo /tmp/alicexty/PythiaFIFO
$ ./PythiaGen PythiaGen.cmnd /tmp/alicexty/PythiaFIFO > PythiaOutput.log 2>&1 &
$ rivet --pwd --analysis=LHCB_2019_I1730448 /tmp/alicexty/PythiaFIFO

```
replacing `<username>` with your username. After it has finished running, a `Rivet.yoda` file will be created with the simulation histograms.

You can then compare the simulation to reference data by doing
```
$ rivet-mkhtml Rivet.yoda --errs --pwd
```
To plot a specific plot you can use the `-m` option, e.g. `-m /LHCB_2025_I2922449/d01-x01-y01`