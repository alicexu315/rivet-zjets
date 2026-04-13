// -*- C++ -*-
#include <utility>   // std::pair
#include <assert.h>  // For checking collision system and COM energy
#include <string>
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"       // For creating Filter for jet Projector
//#include "fastjet/tools/Filter.hh"   // For applying reclustering/grooming in Projector
//#include "fastjet/contrib/LundGenerator.hh"   // For Soft Drop grooming

namespace Rivet {
  
  class JetInfo : public fastjet::PseudoJet::UserInfoBase {
  public:
    JetInfo(const int &_pid) : pid(_pid) {};
    int get_pid() const {
      return this->pid;
    }
  private:
    int pid = -999;
  };

  /// @brief Add a short analysis description here
  class LHCB_2019_I1730448 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2019_I1730448);


    /// @name Analysis methods
    /// @{

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /// Book histograms and initialise projections before the run
    void init() {

      // Check that we have the right collision system and sqrt(s)
      const ParticlePair &beam_pair = beams();
      assert(beam_pair.first.pid() == PID::PROTON && beam_pair.second.pid() == PID::PROTON);
      assert(this->isCompatibleWithSqrtS(8000 * GeV));

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within the given eta acceptance
      Cut particle_selector = (Cuts::eta > ETA_MIN_PARTICLES) & (Cuts::eta < ETA_MAX_PARTICLES);
      const FinalState fs(particle_selector);
      declare(fs, "fs_particles");

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.5
      // muons are included and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, JET_R, JetAlg::Muons::ALL, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");

      Cut Z_boson_selector = particle_selector & (Cuts::abspid == PID::ZBOSON) &
                             (Cuts::pT > PT_MIN_ZBOSONS) & (Cuts::pT < PT_MAX_ZBOSONS) &
                             (Cuts::mass > MASS_MIN_ZBOSONS) & (Cuts::mass < MASS_MAX_ZBOSONS);

      // Find *decayed* charged Z bosons from event
      const UnstableParticles decayed_Z_bosons(Z_boson_selector);
      declare(decayed_Z_bosons, "decayed_Z_bosons");

      // Missing momentum
      declare(MissingMomentum(fs), "MET");

      // Book histograms
      for (int i = 0; i < N_PT_BIN_EDGES - 1; i++) {
        std::string pT_name = std::to_string(PT_BIN_EDGES[i]) + '_' + std::to_string(PT_BIN_EDGES[i+1]);

        // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
        book(_h["unfold1d_z_" + pT_name],              1, 1, 1 + i);  // HEPData Table 0-2
        book(_h["unfold1d_jt_" + pT_name] ,            2, 1, 1 + i);  // HEPData Table 3-5
        book(_h["unfold1d_r_" + pT_name],              3, 1, 1 + i);  // HEPData Table 6-8
      }

      return;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles fs_particles = apply<FinalState>(event, "fs_particles").particles();

      //Retrieve *decayed* Z bosons
      Particles decayed_Z_bosons = apply<UnstableParticles>(event, "decayed_Z_bosons").particles();
      if (DEBUG_LEVEL > 0) std::cout << "Unstable Z count = " << decayed_Z_bosons.size() << std::endl;
      Particles Z_bosons;  // Z bosons which decayed to mu+ mu-
      for (const Particle& decayed_Z_boson : decayed_Z_bosons) {
        Z_bosons.push_back(decayed_Z_boson);
        /*bool has_mup = false, has_mum = false;
        Particles Z_boson_children = decayed_Z_boson.children(Cuts::OPEN);
        if (DEBUG_LEVEL > 0) std::cout << "children: " << decayed_Z_boson.children().size() << std::endl;
        for (const Particle& child: decayed_Z_boson.children()) {
          if (child.pid() == 13) has_mum = true;
          if (child.pid() == -13) has_mup = true;
        }
        if (has_mup && has_mum) {
          Z_bosons.push_back(decayed_Z_boson);
        }*/
      }

      if (Z_bosons.empty()) { vetoEvent; }

      // Retrieve clustered jets, sorted by pT, with applied rapidity and pT cuts
      Cut jet_selector = (Cuts::eta > ETA_MIN_JETS) & (Cuts::eta < ETA_MAX_JETS) &
                         (Cuts::pT > PT_MIN_JETS) & (Cuts::pT < PT_MAX_JETS);
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(jet_selector);
      if (jets.empty()) { vetoEvent; }

      // Create a vector of jets that are back to back with a Z boson
      Jets Z_jets;
      for (const Jet& jet : jets) {
        for (const Particle& Z_boson : Z_bosons) {
          if (std::abs(Z_boson.phi()-jet.phi()) > (7*M_PI)/8 && std::abs(Z_boson.phi()-jet.phi()) < (9*M_PI)/8) {
            Z_jets.push_back(jet);
          }
        }
      }
      if (Z_jets.empty()) { vetoEvent; }

      if (DEBUG_LEVEL > 1) std::cout << "Z jets found" << std::endl;

      ////////////////////////////////////////////////////////////////////////////////////////

      // Loop over jets and apply operations & fill histograms
      for (const Jet& Z_jet : Z_jets) {
        for (const PseudoJet& constituent : Z_jet.pseudojet().constituents()) {

        // Fill histograms
        if (DEBUG_LEVEL > 1) std::cout << "Filling histograms" << std::endl;
        if ((std::sqrt(constituent.modp2()) > P_MIN_HADRONS) && (std::sqrt(constituent.modp2()) < P_MAX_HADRONS) && (constituent.pt() > PT_MIN_HADRONS) &&
            (Z_jet.pseudojet().delta_R(constituent) < 0.5)) {
          double num_z = Z_jet.pseudojet().px()*constituent.px() + Z_jet.pseudojet().py()*constituent.py() + Z_jet.pseudojet().pz()*constituent.pz();
          double den_z = Z_jet.pseudojet().px()*Z_jet.pseudojet().px() + Z_jet.pseudojet().py()*Z_jet.pseudojet().py() + Z_jet.pseudojet().pz()*Z_jet.pseudojet().pz();
          double num_jt = std::sqrt(std::pow(Z_jet.pseudojet().py()*constituent.pz()-Z_jet.pseudojet().pz()*constituent.py(), 2.0) + std::pow(Z_jet.pseudojet().pz()*constituent.px()-Z_jet.pseudojet().px()*constituent.pz(), 2.0) + std::pow(Z_jet.pseudojet().px()*constituent.py()-Z_jet.pseudojet().py()*constituent.px(), 2.0));
          double den_jt = std::sqrt(std::pow(Z_jet.pseudojet().px(), 2.0) + std::pow(Z_jet.pseudojet().py(), 2.0) + std::pow(Z_jet.pseudojet().pz(), 2.0));
          _h[get_histo_name(Z_jet.pT(), true, false, false)]->fill(num_z / den_z); // Fill histogram with longitudinal momentum
          _h[get_histo_name(Z_jet.pT(), false, true, false)]->fill(num_jt / den_jt); // Fill histogram with transverse momentum
          _h[get_histo_name(Z_jet.pT(), false, false, true)]->fill(std::sqrt(std::pow(Z_jet.pseudojet().phi()-constituent.phi(), 2.0) + std::pow(Z_jet.pseudojet().rap()-constituent.rap(), 2.0))); // Fill histogram with radial profile distribution
          }
        }
      }

      return;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////
    /// Normalise histograms etc., after the run
    void finalize() {

      for (int i = 0; i < N_PT_BIN_EDGES - 1; i++) {
        std::string pT_name = std::to_string(PT_BIN_EDGES[i]) + '_' + std::to_string(PT_BIN_EDGES[i+1]);

        normalize(_h["unfold1d_z_" + pT_name]);  // HEPData Table 0-2
        normalize(_h["unfold1d_jt_" + pT_name]);  // HEPData Table 3-5
        normalize(_h["unfold1d_r_" + pT_name]);  // HEPData Table 6-8
      }

      /*
      normalize(_h["XXXX"]); // normalize to unity
      normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in pb (no cuts)
      scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)
      */

      return;
    }

    /// @}

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Helper functions

    // Function which returns internal histogram name for retrieving from map
    std::string get_histo_name(const double & jet_pt, const bool h_z, const bool h_jt, const bool h_r) {
      int i = get_max_pT_index(jet_pt);
      std::string name = "unfold1d_";
      if (h_z) { name += "z_"; }
      if (h_jt) { name += "jt_"; }
      if (h_r) { name += "r_"; }
      name += std::to_string(PT_BIN_EDGES[i-1]) + '_' + std::to_string(PT_BIN_EDGES[i]);
      return name;
    }

    int get_max_pT_index(const double & jet_pt) {
      assert(jet_pt >= PT_MIN_JETS && jet_pt < PT_MAX_JETS);
      int i = 0;
      while (jet_pt >= PT_BIN_EDGES[++i]) { }
      return i;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Analysis parameters

    const int DEBUG_LEVEL = 1;             // For some helpful debugging print functions. Between 0-2

    // Particle reconstruction
    const double ETA_MIN_PARTICLES = 2.;   // Min pseudorapidity of final-state particles
    const double ETA_MAX_PARTICLES = 4.5;   // Max "

    // Jet reconstruction
    const double JET_R = 0.5;              // Jet resolution parameter for anti-kT
    const double ETA_MIN_JETS = 2.5;  // Min rapidity of constructed B-jets
    const double ETA_MAX_JETS = 4.;  // Max "
    const double PT_MIN_JETS = 20.;        // Min jet transverse momentum
    const double PT_MAX_JETS = 100.;       // Max "

    // HF reconstruction
    const std::vector<PdgId> CH_B_MESON_PID = { PID::BPLUS, PID::BMINUS };

    // Z Boson
    const double MASS_MIN_ZBOSONS = 60.;
    const double MASS_MAX_ZBOSONS = 120.;
    const double PT_MIN_ZBOSONS = 0.;
    const double PT_MAX_ZBOSONS = 100.;

    // Charged Hadron
    const double P_MIN_HADRONS = 4.;
    const double P_MAX_HADRONS = 1000.;
    const double PT_MIN_HADRONS = 0.25;

    // Jet pT bin edges
    const int N_PT_BIN_EDGES = 4;
    const std::vector<int> PT_BIN_EDGES = {20, 30, 50, 100};

    /// @name Histograms
    /// @{
    map<std::string, Histo1DPtr> _h;
    map<std::string, Profile1DPtr> _p;
    map<std::string, CounterPtr> _c;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(LHCB_2019_I1730448);

}
