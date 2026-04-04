// C++ script for generating Pythia events for testing Rivet
// Based on main42.cc by Mikhail Kirsanov <Mikhail.Kirsanov@cern.ch>
// Modified by Ezra D. Lesser <ezra.lesser@cern.ch>

// To compile: see compile_PythiaGen_cc.cmd

// To use:  (replace elesser with your username)
// $ mkfifo /tmp/elesser/PythiaFIFO
// $ ./PythiaGen PythiaGen.cmnd /tmp/elesser/PythiaFIFO > PythiaOutput.log 2>&1 &
// $ rivet --pwd --analysis=LHCB_2025_I2922449 /tmp/elesser/PythiaFIFO

#include <iostream>
#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"
#include "Pythia8Plugins/HepMC2.h"
//#include "Pythia8Plugins/Pythia8Rivet.h"  // For compiling PYTHIA with RIVET

// Search pythia event for a specific PID
bool search_pythia(const Pythia8::Pythia & pythia, const int require_pid, std::vector<int>& satisfier_ips) {
    bool pid_found = false;
    for (int ip = 0; ip < pythia.event.size(); ip++) {
        if (pythia.event[ip].idAbs() == require_pid) {
            pid_found = true;
            satisfier_ips.push_back(ip);
            break;
        }
    }
    return pid_found;
}

////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
    const int DEBUG_LEVEL = 1;

    if (DEBUG_LEVEL > 0) std::cout << "PYTHIA WORKING" << std::endl;
    // Check that correct number of command-line arguments
    if (argc != 3) {
        std::cerr << " Unexpected number of command-line arguments. \n You are"
                  << " expected to provide one input and one output file name. \n"
                  << " Program stopped! " << std::endl;
        return 1;
    }

    // Check that the provided input name corresponds to an existing file.
    std::ifstream is(argv[1]);
    if (!is) {
        std::cerr << " Command-line file " << argv[1] << " was not found. \n"
                  << " Program stopped! " << std::endl;
        return 1;
    }

    // Confirm that external files will be used for input and output.
    std::cout << "\n >>> PYTHIA settings will be read from file " << argv[1]
              << " <<< \n >>> HepMC events will be written to file "
              << argv[2] << " <<< \n" << std::endl;

    // Interface for conversion from Pythia8::Event to HepMC event.
    HepMC::Pythia8ToHepMC ToHepMC;

    // Specify file where HepMC events will be stored.
    HepMC::IO_GenEvent ascii_io(argv[2], std::ios::out);

    // Generator.
    Pythia8::Pythia pythia;

    // Read in commands from external file.
    pythia.readFile(argv[1]);

    // Extract settings to be used in the main program.
    int    nEvent    = pythia.mode("Main:numberOfEvents");
    int    nAbort    = pythia.mode("Main:timesAllowErrors");
    if (DEBUG_LEVEL > 0) std::cout << "P: settings extracted" << std::endl;

    // Initialization.
    pythia.init();
    if (DEBUG_LEVEL > 0) std::cout << "P: initialised" << std::endl;

    // Begin event loop.
    int iAbort = 0;
    int iEvent = 0;
    while (iEvent < nEvent) {

        // Generate event.
        if (!pythia.next()) {

            // If failure because reached end of file then exit event loop.
            if (pythia.info.atEndOfFile()) {
                std::cout << " Aborted since reached end of Les Houches Event File\n";
                break;
            }

            // First few failures write off as "acceptable" errors, then quit.
            if (++iAbort < nAbort) continue;
            std::cout << " Event generation aborted prematurely, owing to error!\n";
            break;
        }

        /* Generate events until you get one with a desired particle.
        const int DESIRED_PID = 5;  // b quark
        bool desired_parton_found = false;
        for (int iParticle = 0; iParticle < pythia.event.size(); iParticle++) {
          if (pythia.event[iParticle].idAbs() == DESIRED_PID) {
            desired_parton_found = true;
            break;
          }
        }
        if (!desired_parton_found) continue;
        */

        // Update hadronization.
        Pythia8::Event savedEvent = pythia.event;
        bool hstatus = pythia.forceHadronLevel();
        if (!hstatus) continue;
        if (DEBUG_LEVEL > 0) std::cout << "P: hadronization updated" << std::endl;

        // Check if there's a Z boson
        const int REQUIRE_PID = 23;  // Z boson
        std::vector<int> satisfier_ips = {};        // Index of particle with REQUIRE_PID
        bool pid_found = search_pythia(pythia, REQUIRE_PID, satisfier_ips);
        if (!pid_found) {
            continue;
        }
        if (DEBUG_LEVEL > 0) std::cout << "P: Z boson found" << std::endl;
        
        // If the particle is not in desired kinematic region, skip event & regenerate
        const double MIN_RAPIDITY_ACCEPTED = 1.5;
        const double MAX_RAPIDITY_ACCEPTED = 5.;
        bool rapidity_satisfied = false;
        for (int i = 0; i < satisfier_ips.size(); i++) {
          double satisfier_rapidity = pythia.event[satisfier_ips[i]].y();
          if ( (satisfier_rapidity >= MIN_RAPIDITY_ACCEPTED) && (satisfier_rapidity < MAX_RAPIDITY_ACCEPTED) ) {
              rapidity_satisfied = true;
              break;
          }
        }
        if (!rapidity_satisfied) continue;

        // Construct new empty HepMC event and fill it.
        // Units will be as chosen for HepMC build, but can be changed
        // by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
        HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
        ToHepMC.fill_next_event( pythia, hepmcevt );
        if (DEBUG_LEVEL>0) std::cout << "P: Event filled" << std::endl;

        // Write the HepMC event to file. Done with it.
        ascii_io << hepmcevt;
        delete hepmcevt;

        iEvent++;

    // End of event loop. Statistics.
    }
    pythia.stat();

    // Done.
    return 0;

}
