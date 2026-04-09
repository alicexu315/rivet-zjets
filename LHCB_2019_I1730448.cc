// -*- C++ -*-
#include <utility>   // std::pair
#include <assert.h>  // For checking collision system and COM energy
#include <string>
#include <cmath>
#include <math.h>
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
//#include "fastjet/Selector.hh"       // For creating Filter for jet Projector
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
    //FastJets jetfs(fs, FastJets::ANTIKT, JET_R, JetAlg::Muons::ALL, JetAlg::Invisibles::NONE);
    //declare(jetfs, "jets");

    // Create a list of jets reclustered using C/A with WTA axis
    /* Note: this method reclusters jets from full event, so it's not as efficient
    fastjet::Selector WTA_jet_selector =
      fastjet::SelectorPtMin(PT_MIN_JETS) && fastjet::SelectorRapRange(RAPIDITY_MIN_JETS, RAPIDITY_MAX_JETS);
    fastjet::Filter* WTA_reclusterer = new fastjet::Filter(WTA_JET_DEF, WTA_jet_selector);
    FastJets jetfs_WTA(fs, FastJets::ANTIKT, JET_R, JetAlg::Muons::ALL, JetAlg::Invisibles::NONE);
    jetfs_WTA.addTrf(WTA_reclusterer);
    declare(jetfs_WTA, "jets_WTA");
    */

    // FinalState of *undecayed* Z bosons
    // NOTE: To use this, Z bosons must be set to be undecayed in MC!
    Cut Z_boson_selector = particle_selector & (Cuts::abspid == PID::ZBOSON)
      & (Cuts::pT > PT_MIN_ZBOSONS) & (Cuts::pT < PT_MAX_ZBOSONS)
      & (Cuts::mass > MASS_MIN_ZBOSONS) & (Cuts::mass < MASS_MAX_ZBOSONS);
    //const FinalState fs_Z_bosons(Z_boson_selector);
    //declare(fs_Z_bosons, "Z_bosons");
    if (DEBUG_LEVEL > 1) std::cout << "Final state declared" << std::endl;

    // Find *decayed* Z bosons from event
    const UnstableParticles decayed_Z_bosons(Z_boson_selector);
    declare(decayed_Z_bosons, "decayed_Z_bosons");
    if (DEBUG_LEVEL > 1) std::cout << "Decayed Z bosons declared" << std::endl;

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
    //if(DEBUG_LEVEL > 1) std::cout << "Histograms booked" << std::endl;

    return;
  };


  /////////////////////////////////////////////////////////////////////////////////////////////////
  /// Perform the per-event analysis
  void analyze(const Event& event) {
    //if(DEBUG_LEVEL > 1) std::cout << "analyze() is called" << std::endl;
    /* Retrieve *undecayed* Z bosons
    Particles Z_bosons = apply<FinalState>(event, "Z_bosons").particles();
    if (Z_bosons.empty()) { vetoEvent; }

    if (DEBUG_LEVEL > 0) {
      std::cout << "Z bosons(s) found with ";
      bool add_comma = false;
      for (const Particle& Z_boson : Z_bosons) {
        if (add_comma) std::cout << ", ";
        std::cout << "rapidity=" << Z_boson.rap() << " and pT=" << Z_boson.pt();
        add_comma = true;
      }
      std::cout << std::endl;
    }
    */
    std::vector<PseudoJet> trimmed_particles;
    Particles fs_particles = apply<FinalState>(event, "fs_particles").particles();
    for(const Particle& fs_particle : fs_particles) {
      fastjet::PseudoJet AliceJet = fs_particle.pseudojet();
      //JetInfo* jet_pid = new JetInfo(fs_particle.pid());
      //AliceJet.set_user_info(jet_pid);
      AliceJet.set_user_info_shared_ptr(fastjet::SharedPtr<fastjet::PseudoJet::UserInfoBase>(new JetInfo(fs_particle.pid())));

      trimmed_particles.push_back(AliceJet);
    }

    JetDefinition jet_def(fastjet::antikt_algorithm, JET_R);
    ClusterSequence cs(trimmed_particles, jet_def);
    std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

    // Retrieve *decayed* Z bosons
    Particles decayed_Z_bosons = apply<UnstableParticles>(event, "decayed_Z_bosons").particles();
    Particles Z_bosons;  // Z bosons which decayed to muon antimuon
    for (const Particle& decayed_Z_boson : decayed_Z_bosons) {
      /*
      Particles Z_boson_children = decayed_Z_boson.children(Cuts::OPEN);
      std::cout << "Children: " << Z_boson_children.size() << endl;
      std::cout << decayed_Z_boson.abspid() << std::endl;
      if (Z_boson_children.size() == 2 &&
        Z_boson_children[0].eta() > ETA_MIN_PARTICLES && Z_boson_children[0].eta() < ETA_MAX_PARTICLES &&
        Z_boson_children[1].eta() > ETA_MIN_PARTICLES && Z_boson_children[1].eta() < ETA_MAX_PARTICLES) {
          //(Z_boson_children[0].abspid() == PID::MUON  && Z_boson_children[1].abspid() == PID::MUON) ) {
      */
          Z_bosons.push_back(decayed_Z_boson);
      //}
    }
    if (Z_bosons.empty()) { vetoEvent; }

    // Retrieve clustered jets, sorted by pT, with applied rapidity and pT cuts
    fastjet::Selector jet_selector = fastjet::SelectorPtRange(PT_MIN_JETS, PT_MAX_JETS) && fastjet::SelectorEtaRange(ETA_MIN_JETS, ETA_MAX_JETS);
    jets = jet_selector(jets);
    //Cut jet_selector = (Cuts::rap > RAPIDITY_MIN_JETS) & (Cuts::rap < RAPIDITY_MAX_JETS)
      //& (Cuts::pT > PT_MIN_JETS) & (Cuts::pT < PT_MAX_JETS);
    if (jets.empty()) { vetoEvent; }

    // Create a vector of jets recoiling against Z bosons
    Jets Z_jets;
    Particles Z_boson_matches;
    for (const PseudoJet& jet : jets) {
      //const Jet& jet = jets.front();
      for(const Particle& Z_boson : Z_bosons) {
        if (std::abs(Z_boson.phi()-jet.phi()) > (7*M_PI)/8 && std::abs(Z_boson.phi()-jet.phi()) < (9*M_PI)/8) {  // Jet is on the azimuthal away side of Z
          Z_jets.push_back(jet);
          Z_boson_matches.push_back(Z_boson);
        }
      }
    }
    if (Z_jets.empty()) { vetoEvent; }

    if (DEBUG_LEVEL > 1) std::cout << "Z jets found" << std::endl;

    ////////////////////////////////////////////////////////////////////////////////////////
      // Fill histograms
    if (DEBUG_LEVEL > 1) std::cout << "Filling histograms" << std::endl;

     for (long unsigned int i = 0; i < Z_jets.size(); i++) {
        const Jet& Z_jet = Z_jets[i];
        const Particle& Z_boson_match = Z_boson_matches[i];
        //std::vector<PseudoJet> trimmed_constituents;
        for (const PseudoJet& constituent : Z_jet.pseudojet().constituents()) {
          if (!constituent.has_user_info()) continue;
          const JetInfo& myinfo = constituent.user_info<JetInfo>();
          std::cout << "PID: " << myinfo.get_pid() << std::endl;
          if (((abs(myinfo.get_pid()) == 211) || (abs(myinfo.get_pid()) == 321) || (abs(myinfo.get_pid()) == 2212)) &&
              (std::sqrt(constituent.modp()) > P_MIN_HADRONS) && (std::sqrt(constituent.modp()) < P_MAX_HADRONS) && (constituent.pt() > PT_MIN_HADRONS) &&
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

    /*for (const PseudoJet& jet : jets) {
      std::cout << "jet pT=" << jet.pt() << " eta=" << jet.eta() << std::endl;
      for (const Particles& Z_boson : Z_bosons) {
        //if (std::abs(Z_boson.phi()-jet.phi()) > (7*M_PI)/8 && std::abs(Z_boson.phi()-jet.phi()) < (9*M_PI)/8) {  // Jet is on the azimuthal away side of Z
        if (deltaPhi(Z_boson, jet) > (7*M_PI)/8) {
          for (const PseudoJet& constituent : jet.constituents()) {
            //JetInfo myinfo = constituent.user_info<JetInfo>();
            if (!constituent.has_user_info()) continue;
            const JetInfo& myinfo = constituent.user_info<JetInfo>();
            std::cout << "PID: " << myinfo.get_pid() << std::endl;
            if (((abs(myinfo.get_pid()) == 211) || (abs(myinfo.get_pid()) == 321) || (abs(myinfo.get_pid()) == 2212)) &&
                (constituent.modp() > P_MIN_HADRONS) && (constituent.modp() < P_MAX_HADRONS) && (constituent.pt() > PT_MIN_HADRONS) &&
                (jet.delta_R(constituent) < 0.5)) {
              double num_z = jet.px()*constituent.px() + jet.py()*constituent.py() + jet.pz()*constituent.pz();
              double den_z = jet.px()*jet.px() + jet.py()*jet.py() + jet.pz()*jet.pz();
              double num_jt = std::sqrt(std::pow(jet.py()*constituent.pz()-jet.pz()*constituent.py(), 2.0) + std::pow(jet.pz()*constituent.px()-jet.px()*constituent.pz(), 2.0) + std::pow(jet.px()*constituent.py()-jet.py()*constituent.px(), 2.0));
              double den_jt = std::sqrt(std::pow(jet.px(), 2.0) + std::pow(jet.py(), 2.0) + std::pow(jet.pz(), 2.0));
              _h[get_histo_name(jet.perp(), true, false, false)]->fill(num_z / den_z); // Fill histogram with longitudinal momentum
              _h[get_histo_name(jet.perp(), false, true, false)]->fill(num_jt / den_jt); // Fill histogram with transverse momentum
              _h[get_histo_name(jet.perp(), false, false, true)]->fill(std::sqrt(std::pow(jet.phi()-constituent.phi(), 2.0) + std::pow(jet.rap()-constituent.rap(), 2.0))); // Fill histogram with radial profile distribution
            }
          }
        }
      }
    }*/
    
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

  const int DEBUG_LEVEL = 0;             // For some helpful debugging print functions. Between 0-2

  // Particle reconstruction
  const double ETA_MIN_PARTICLES = 2.;   // Min pseudorapidity of final-state particles
  const double ETA_MAX_PARTICLES = 4.5;   // Max "

  // Jet reconstruction
  const double JET_R = 0.5;              // Jet resolution parameter for anti-kT
  const double ETA_MIN_JETS = 2.5;  // Min rapidity of constructed Z-jets
  const double ETA_MAX_JETS = 4.0;  // Max "
  const double PT_MIN_ZBOSONS = 0;       // Min Z boson transverse momentum
  const double PT_MAX_ZBOSONS = 100;     // Max "
  const double MASS_MIN_ZBOSONS = 60;    // Min invariant mass of Z boson
  const double MASS_MAX_ZBOSONS = 120;   // Max "
  const double PT_MIN_JETS = 20.;        // Min jet transverse momentum
  const double PT_MAX_JETS = 100.;       // Max "
  const double PT_MIN_HADRONS = 0.250;   // Min charged hadron transverse momentum
  const double P_MIN_HADRONS = 4;
  const double P_MAX_HADRONS = 1000;

  // HF reconstruction
  const std::vector<PdgId> CH_B_MESON_PID = { PID::ZBOSON };

  // Jet pT bin edges
  const int N_PT_BIN_EDGES = 4;
  const std::vector<int> PT_BIN_EDGES = {20,30,50,100};

  /// @name Histograms
  /// @{
  map<std::string, Histo1DPtr> _h;
  map<std::string, Profile1DPtr> _p;
  map<std::string, CounterPtr> _c;
  /// @}

};


RIVET_DECLARE_PLUGIN(LHCB_2019_I1730448);

}
