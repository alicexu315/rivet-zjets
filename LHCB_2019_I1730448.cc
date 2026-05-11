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
#include "Rivet/Projections/ZFinder.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"       // For creating Filter for jet Projector

namespace Rivet {

  /// @brief pp collisions at 8 TeV for charged particles in LHCb acceptance
  class LHCB_2019_I1730448 : public Analysis {
  public:

    // Attaches constituent PID to a PseudoJet
    class JetInfo : public fastjet::PseudoJet::UserInfoBase {
    public:
      JetInfo(const int &_pid) : pid(_pid) {};
      int get_pid() const {
        return this->pid;
      }
    private:
      int pid = -999;
    };

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

      Cut Z_boson_selector = particle_selector & (Cuts::abspid == PID::ZBOSON) &
                            (Cuts::pT > PT_MIN_ZBOSONS) & (Cuts::pT < PT_MAX_ZBOSONS) &
                            (Cuts::mass > MASS_MIN_ZBOSONS) & (Cuts::mass < MASS_MAX_ZBOSONS);

      // Find *decayed* Z bosons from dimuon children
      ZFinder zfinder(fs, particle_selector, PID::MUON, MASS_MIN_ZBOSONS, MASS_MAX_ZBOSONS, 0.0);
      declare(zfinder, "ZFinder");

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
      
      if (DEBUG_LEVEL > 0) {
        int nmuons = 0;
          for (const Particle& p : fs_particles) {
            if (abs(p.pid()) == 13) nmuons++;
          }
          std::cout << "fs size: " << fs_particles.size() << " | muons: " << nmuons << std::endl;
      }

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.5
      std::vector<PseudoJet> trimmed_particles;
      for(const Particle& fs_particle : fs_particles) {
        fastjet::PseudoJet InfoJet = fs_particle.pseudojet();
        InfoJet.set_user_info_shared_ptr(fastjet::SharedPtr<fastjet::PseudoJet::UserInfoBase>(new JetInfo(fs_particle.pid())));
        trimmed_particles.push_back(InfoJet);
      }
      JetDefinition jet_def(fastjet::antikt_algorithm, JET_R);
      ClusterSequence cs(trimmed_particles, jet_def);
      std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

      // Retrieve *decayed* Z bosons and muon children
      const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");
      if (zfinder.bosons().empty()) vetoEvent;

      const Particles& Z_bosons = zfinder.bosons();
      const Particles& muons = zfinder.constituents();  // mu+ and mu-

      // Retrieve clustered jets, sorted by pT, with applied rapidity and pT cuts
      fastjet::Selector jet_selector = fastjet::SelectorPtRange(PT_MIN_JETS, PT_MAX_JETS) && fastjet::SelectorEtaRange(ETA_MIN_JETS, ETA_MAX_JETS);
      jets = jet_selector(jets);
      if (jets.empty()) { vetoEvent; }

      // Create a vector of jets that are back to back with a Z boson, with no muons in jet cone
      Jets Z_jets;
      for (const Jet& jet : jets) {
        for (const Particle& Z_boson : Z_bosons) {
          if ((Z_boson.pT() < PT_MIN_ZBOSONS) || (Z_boson.pT() > PT_MAX_ZBOSONS)) continue; // Check Z pT
          bool muon_in_jet = false;
          for (const Particle& muon : muons) {
            if (deltaR(muon, jet) < JET_R) {
              muon_in_jet = true;
              break;
            }
          }
          if (muon_in_jet) continue;
          if (deltaPhi(Z_boson, jet) > (7*M_PI)/8) {  // Azimuthal cut
            Z_jets.push_back(jet);
            // Recording the number of jets in each bin, for scaling purposes
            if ((jet.pT() >= 20.) && (jet.pT() < 30.)) {
              num_jets_20_30 += 1;
            }
            else if ((jet.pT() >= 30.) && (jet.pT() < 50.)) {
              num_jets_30_50 += 1;
            }
            else if ((jet.pT() >= 50.) && (jet.pT() < 100.)) {
              num_jets_50_100 += 1;
            }
            if (DEBUG_LEVEL > 0) std::cout << "deltaphi: " << deltaPhi(Z_boson, jet) << std::endl;
            if (DEBUG_LEVEL > 0) std::cout << "jet eta: " << jet.eta() << std::endl;
            if (DEBUG_LEVEL > 0) std::cout << "jet pt: " << jet.pT() << std::endl;
          }
        }
      }
      if (Z_jets.empty()) { vetoEvent; }

      if (DEBUG_LEVEL > 1) std::cout << "Z jets found" << std::endl;

      ////////////////////////////////////////////////////////////////////////////////////////

      // Loop over jets and apply operations & fill histograms
      for (const Jet& Z_jet : Z_jets) {
        for (const PseudoJet& constituent : Z_jet.pseudojet().constituents()) {
          const JetInfo& myinfo = constituent.user_info<JetInfo>();
        // Fill histograms
        if (DEBUG_LEVEL > 1) std::cout << "Filling histograms" << std::endl;
        if ((std::sqrt(constituent.modp2()) > P_MIN_HADRONS) && (std::sqrt(constituent.modp2()) < P_MAX_HADRONS) &&
            ((abs(myinfo.get_pid()) == PID::PIPLUS) || (abs(myinfo.get_pid()) == PID::KPLUS) || (abs(myinfo.get_pid()) == PID::PROTON)) &&
            (constituent.pt() > PT_MIN_HADRONS) &&
            (std::sqrt(std::pow(Z_jet.pseudojet().delta_phi_to(constituent), 2) + std::pow(Z_jet.pseudojet().eta()-constituent.eta(), 2)) < JET_R)) {
          double num_z = Z_jet.pseudojet().px()*constituent.px() + Z_jet.pseudojet().py()*constituent.py() + Z_jet.pseudojet().pz()*constituent.pz();
          double den_z = Z_jet.pseudojet().px()*Z_jet.pseudojet().px() + Z_jet.pseudojet().py()*Z_jet.pseudojet().py() + Z_jet.pseudojet().pz()*Z_jet.pseudojet().pz();
          double num_jt = std::sqrt(std::pow(Z_jet.pseudojet().py()*constituent.pz()-Z_jet.pseudojet().pz()*constituent.py(), 2)
           + std::pow(Z_jet.pseudojet().pz()*constituent.px()-Z_jet.pseudojet().px()*constituent.pz(), 2)
           + std::pow(Z_jet.pseudojet().px()*constituent.py()-Z_jet.pseudojet().py()*constituent.px(), 2));
          double den_jt = std::sqrt(std::pow(Z_jet.pseudojet().px(), 2) + std::pow(Z_jet.pseudojet().py(), 2) + std::pow(Z_jet.pseudojet().pz(), 2));
          double r = std::sqrt(std::pow(Z_jet.pseudojet().delta_phi_to(constituent), 2) + std::pow(Z_jet.pseudojet().rap()-constituent.rap(), 2));
          _h[get_histo_name(Z_jet.pT(), true, false, false)]->fill(num_z / den_z); // Fill histogram with longitudinal momentum
          _h[get_histo_name(Z_jet.pT(), false, true, false)]->fill(num_jt / den_jt); // Fill histogram with transverse momentum
          _h[get_histo_name(Z_jet.pT(), false, false, true)]->fill(r); // Fill histogram with radial profile distribution
          }
        }
      }

      return;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////
    /// Scale histograms etc., after the run
    void finalize() {
      for (int i = 0; i < N_PT_BIN_EDGES - 1; i++) {
        std::string pT_name = std::to_string(PT_BIN_EDGES[i]) + '_' + std::to_string(PT_BIN_EDGES[i+1]);
        double njet = 0.;
        if (PT_BIN_EDGES[i] == 20) njet = num_jets_20_30;
        if (PT_BIN_EDGES[i] == 30) njet = num_jets_30_50;
        if (PT_BIN_EDGES[i] == 50) njet = num_jets_50_100;
        if (njet == 0.) continue;
        scale(_h["unfold1d_z_" + pT_name], 1.0/njet);  // HEPData Table 0-2
        scale(_h["unfold1d_jt_" + pT_name], 1.0/njet);  // HEPData Table 3-5
        scale(_h["unfold1d_r_" + pT_name], 1.0/njet);  // HEPData Table 6-8
      }

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
    const double ETA_MAX_PARTICLES = 4.5;  // Max "

    // Jet reconstruction
    const double JET_R = 0.5;              // Jet resolution parameter for anti-kT
    const double ETA_MIN_JETS = 2.5;       // Min rapidity of constructed Z-jets
    const double ETA_MAX_JETS = 4.;        // Max "
    const double PT_MIN_JETS = 20.;        // Min jet transverse momentum
    const double PT_MAX_JETS = 100.;       // Max "
    double num_jets_20_30 = 0.;            // Total number of Z-tagged jets from the 20<pT<30 GeV bin
    double num_jets_30_50 = 0.;            // " 30<pT<50 GeV bin
    double num_jets_50_100 = 0.;           // " 50<pT<100 GeV bin

    // Z Boson
    const double MASS_MIN_ZBOSONS = 60.;   // Min invariant mass of Z bosons
    const double MASS_MAX_ZBOSONS = 120.;  // Max "
    const double PT_MIN_ZBOSONS = 0.;      // Min Z boson transverse momentum
    const double PT_MAX_ZBOSONS = 100.;    // Max "

    // Charged Hadron
    const double P_MIN_HADRONS = 4.;       // Min momentum of charged hadrons
    const double P_MAX_HADRONS = 1000.;    // Max "
    const double PT_MIN_HADRONS = 0.25;    // Min transverse momentum of charged hadrons

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
