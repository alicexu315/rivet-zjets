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
#include "fastjet/contrib/LundGenerator.hh"   // For Soft Drop grooming

namespace Rivet {


  /// @brief Add a short analysis description here
  class LHCB_2025_I2922449 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2025_I2922449);


    /// @name Analysis methods
    /// @{

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /// Book histograms and initialise projections before the run
    void init() {

      // Check that we have the right collision system and sqrt(s)
      const ParticlePair &beam_pair = beams();
      assert(beam_pair.first.pid() == PID::PROTON && beam_pair.second.pid() == PID::PROTON);
      assert(this->isCompatibleWithSqrtS(13000 * GeV));

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within the given eta acceptance
      Cut particle_selector = (Cuts::eta > ETA_MIN_PARTICLES) & (Cuts::eta < ETA_MAX_PARTICLES);
      const FinalState fs(particle_selector);

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.5
      // muons are included and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, JET_R, JetAlg::Muons::ALL, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");

      // Create a list of jets reclustered using C/A with WTA axis
      /* Note: this method reclusters jets from full event, so it's not as efficient
      fastjet::Selector WTA_jet_selector =
        fastjet::SelectorPtMin(PT_MIN_JETS) && fastjet::SelectorRapRange(RAPIDITY_MIN_JETS, RAPIDITY_MAX_JETS);
      fastjet::Filter* WTA_reclusterer = new fastjet::Filter(WTA_JET_DEF, WTA_jet_selector);
      FastJets jetfs_WTA(fs, FastJets::ANTIKT, JET_R, JetAlg::Muons::ALL, JetAlg::Invisibles::NONE);
      jetfs_WTA.addTrf(WTA_reclusterer);
      declare(jetfs_WTA, "jets_WTA");
      */

      // FinalState of *undecayed* charged Z bosons
      // NOTE: To use this, charged Z bosons must be set to be undecayed in MC!
      Cut Z_boson_selector = particle_selector & (Cuts::abspid == PID::BPLUS);
      const FinalState fs_Z_bosons(Z_boson_selector);
      declare(fs_Z_bosons, "Z_bosons");

      /* Find *decayed* charged Z bosons from event
      const UnstableParticles decayed_Z_bosons(Z_boson_selector);
      */

      // FinalState of direct photons and bare muons and electrons in the event
      DirectFinalState photons(Cuts::abspid == PID::PHOTON);
      DirectFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);

      // Dress the bare direct leptons with direct photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons
      Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
      DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
      declare(dressed_leps, "leptons");

      // Missing momentum
      declare(MissingMomentum(fs), "MET");

      // Book histograms
      for (int i = 0; i < N_PT_BIN_EDGES - 1; i++) {
        std::string pT_name = "pT" + std::to_string(PT_BIN_EDGES[i]) + '-' + std::to_string(PT_BIN_EDGES[i+1]);

        // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
        book(_h[pT_name],              13 + i, 1, 1);  // HEPData Table 12-17
        book(_h[pT_name + "_gr"] ,      1 + i, 1, 1);  // HEPData Table 0-5
        book(_h[pT_name + "_WTA"],     37 + i, 1, 1);  // HEPData Table 36-41
        book(_h[pT_name + "_gr_WTA"],  19 + i, 1, 1);  // HEPData Table 18-23
      }

      return;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////
    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve *undecayed* charged Z bosons
      Particles Z_bosons = apply<FinalState>(event, "Z_bosons").particles();
      if (Z_bosons.empty()) { vetoEvent; }

      if (DEBUG_LEVEL > 0) {
        std::cout << "Z boson(s) found with ";
        bool add_comma = false;
        for (const Particle& Z_boson : Z_bosons) {
          if (add_comma) std::cout << ", ";
          std::cout << "rapidity=" << Z_boson.rap() << " and pT=" << Z_boson.pt();
          add_comma = true;
        }
        std::cout << std::endl;
      }

      Retrieve *decayed* charged Z bosons
      Particles decayed_Z_bosons = apply<UnstableParticles>(event, "decayed_Z_bosons").particles();
      Particles Z_bosons;  // Charged Z bosons which decayed to Jpsi K
      for (const Particle& decayed_Z_boson : decayed_Z_bosons) {
        Particles Z_boson_children = decayed_Z_boson.children(Cuts::OPEN);
        if (Z_boson_children.size() == 2) &&
            ( (Z_boson_children[0].abspid() == PID::JPSI  && Z_boson_children[1].abspid() == PID::KPLUS) ||
            ( (Z_boson_children[0].abspid() == PID::KPLUS && Z_boson_children[1].abspid() == PID::JPSI ) )) {
          Z_bosons.push_back(decayed_Z_boson);
          // TODO: replace its children in event before jet reconstruction
        }
      }

      // Retrieve clustered jets, sorted by pT, with applied rapidity and pT cuts
      Cut jet_selector = (Cuts::rap > RAPIDITY_MIN_JETS) & (Cuts::rap < RAPIDITY_MAX_JETS)
        & (Cuts::pT > PT_MIN_JETS) & (Cuts::pT < PT_MAX_JETS);
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(jet_selector);
      if (jets.empty()) { vetoEvent; }

      // Create a vector of jets containing a charged Z boson
      Jets B_jets;
      for (const Jet& jet : jets) {
        // for (const Particle& Z_boson : Z_bosons) {
        // if (jet.containsParticle(Besmon)) {
        if (jet.containsParticleId(CH_Z_boson_PID)) {  // Jet contains B+ or B-
          B_jets.push_back(jet);
        }
      }
      if (B_jets.empty()) { vetoEvent; }

      if (DEBUG_LEVEL > 0) std::cout << "B jets found" << std::endl;

      ////////////////////////////////////////////////////////////////////////////////////////

      // Loop over jets and apply operations & fill histograms
      Jets B_jets_WTA;
      for (const Jet& B_jet : B_jets) {
        // Find the WTA axis of the B jet using C/A algorithm
        fastjet::ClusterSequence B_jet_WTA_cs(B_jet.pseudojet().constituents(), WTA_JET_DEF);
        Jet B_jet_WTA = fastjet::sorted_by_pt(B_jet_WTA_cs.inclusive_jets())[0];
        assert(B_jet.pseudojet().constituents().size() == B_jet_WTA.pseudojet().constituents().size());  // Make sure that reclustering preserves all particles

        // Check if charged Z boson is on WTA axis (arXiv:2205.01117)
        bool WTA_tagged = false;
        for (const Particle& Z_boson : Z_bosons) {
          double WTA_distance = B_jet_WTA.pseudojet().delta_R(Z_boson);
          if (WTA_distance < MAX_WTA_DISTANCE) {
            WTA_tagged = true;
            break;
          }
        }

        // Groom the jets with Soft Drop
        Jet B_jet_gr;
        for (fastjet::contrib::LundDeclustering & ld : B_JET_LG.result(B_jet)) {
          if (ld.z() > SD_ZCUT * std::pow(ld.Delta() / JET_R, SD_BETA)) {  // SD condition
            B_jet_gr = ld.pair();
            break;
          }
        }
        // Check that charged Z boson survives
        //bool gr_tagged = B_jet_gr.containsParticleId(CH_Z_boson_PID);  // Doesn't work; particle information is lost
        bool gr_tagged = false;
        for (const Particle& Z_boson : Z_bosons) {
          if (gr_tagged || !B_jet_gr.pseudojet().has_constituents()) break;
          for (const fastjet::PseudoJet& particle : B_jet_gr.pseudojet().constituents()) {
            if ((particle.delta_R(Z_boson) < SMALL_NUMBER) && (abs(particle.pt() - Z_boson.pt()) < SMALL_NUMBER)) {
              // Should be the same particle
              gr_tagged = true;
              break;
            }
          }
        }
        if (DEBUG_LEVEL > 1) {
          if (gr_tagged) { std::cout << "** B jet survived grooming" << std::endl; }
          else { std::cout << "** B jet did NOT survive grooming" << std::endl; }
        }

        // Fill histograms
        if (DEBUG_LEVEL > 0) std::cout << "Filling histograms" << std::endl;
        _h[get_histo_name(B_jet.pT(), false, false)]->fill(B_jet.pseudojet().m() / B_jet.pT());

        if (WTA_tagged) {
          _h[get_histo_name(B_jet.pT(), true, false)]->fill(B_jet.pseudojet().m() / B_jet.pT());
        }

        if (gr_tagged) {
          // Note: Groomed observable uses ungroomed jet pT
          _h[get_histo_name(B_jet.pT(), false, true)]->fill(B_jet_gr.pseudojet().m() / B_jet.pT());
          if (WTA_tagged) {
            _h[get_histo_name(B_jet.pT(), true, true)]->fill(B_jet_gr.pseudojet().m() / B_jet.pT());
          }
        }
      }

      return;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////
    /// Normalise histograms etc., after the run
    void finalize() {

      for (int i = 0; i < N_PT_BIN_EDGES - 1; i++) {
        std::string pT_name = "pT" + std::to_string(PT_BIN_EDGES[i]) + '-' + std::to_string(PT_BIN_EDGES[i+1]);

        double groomed_WTA_tagging_fraction = _h[pT_name + "_gr_WTA"]->integral() / _h[pT_name + "_gr"]->integral();

        normalize(_h[pT_name]);           // HEPData Table 12-17
        normalize(_h[pT_name + "_gr"]);   // HEPData Table 0-5
        normalize(_h[pT_name + "_WTA"]);  // HEPData Table 36-41
        normalize(_h[pT_name + "_gr_WTA"], groomed_WTA_tagging_fraction);  // HEPData Table 18-23
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
    std::string get_histo_name(const double & jet_pt, const bool h_wta, const bool h_gr) {
      int i = get_max_pT_index(jet_pt);
      std::string name = "pT" + std::to_string(PT_BIN_EDGES[i-1]) + '-' + std::to_string(PT_BIN_EDGES[i]);
      if (h_gr) { name += "_gr"; }
      if (h_wta) { name += "_WTA"; }
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
    const double ETA_MIN_PARTICLES = 0.;   // Min pseudorapidity of final-state particles
    const double ETA_MAX_PARTICLES = 7.;   // Max "

    // Jet reconstruction
    const double JET_R = 0.5;              // Jet resolution parameter for anti-kT
    const double RAPIDITY_MIN_JETS = 2.5;  // Min rapidity of constructed B-jets
    const double RAPIDITY_MAX_JETS = 4.0;  // Max "
    const double PT_MIN_JETS = 10.;        // Min jet transverse momentum
    const double PT_MAX_JETS = 100.;       // Max "

    // HF reconstruction
    const std::vector<PdgId> CH_Z_boson_PID = { PID::BPLUS, PID::BMINUS };

    // WTA reclustering
    const double SMALL_NUMBER = 1e-4;
    const double & MAX_WTA_DISTANCE = SMALL_NUMBER;
    const fastjet::JetDefinition WTA_JET_DEF = fastjet::JetDefinition(
      fastjet::cambridge_algorithm, fastjet::JetDefinition::max_allowable_R, fastjet::WTA_pt_scheme);

    // Soft Drop grooming
    const double SD_ZCUT = 0.1;
    const double SD_BETA = 0.;
    const fastjet::contrib::LundGenerator B_JET_LG =
      fastjet::contrib::LundGenerator(fastjet::cambridge_algorithm);

    // Jet pT bin edges
    const int N_PT_BIN_EDGES = 7;
    const std::vector<int> PT_BIN_EDGES = {10, 12, 15, 20, 30, 50, 100};

    /// @name Histograms
    /// @{
    map<std::string, Histo1DPtr> _h;
    map<std::string, Profile1DPtr> _p;
    map<std::string, CounterPtr> _c;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(LHCB_2025_I2922449);

}
