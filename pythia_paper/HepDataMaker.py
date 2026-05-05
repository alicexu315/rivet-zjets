#! /usr/bin/env python

#----------------------------------------------------------------------
# Write HEPData submission
# You can test it here: https://www.hepdata.net/record/sandbox
#
# You will want to edit the tables slightly to add observable-specific info:
#   units, description, etc.
#----------------------------------------------------------------------

import os           # Access to file system
import numpy as np  # Numerical operations
#import ROOT         # Ability to read ROOT files
import hepdata_lib  # Helpful library for HepData output


#----------------------------------------------------------------------
def add_hepdata_table(i_count, observable, min_pt, max_pt):

    table = hepdata_lib.Table(f'Table {i_count}')
    table.keywords["reactions"] = ['P P --> Z jet X']
    table.keywords["cmenergies"] = ['8000.']

    # Set observable-specific info in table
    x_label, y_label = set_hepdata_table_descriptors(
        table, observable, min_pt, max_pt)

    # Create readers to read histograms
    f_final_name = "allhist.root"
    f_final_path = os.path.join("/home", "alicexty", "rivet", "Rivet", "pythia_paper", f_final_name)
    #f_final = ROOT.TFile(f_final_path, "READ")
    hepdata_reader = hepdata_lib.RootFileReader(f_final_path)

    # Define variables
    h_name = "%s_%i_%i" % (observable, min_pt, max_pt)
    h = hepdata_reader.read_hist_1d(h_name)

    n_significant_digits = 3  # Note: use 2 significant digits for uncertainties

    x = hepdata_lib.Variable(x_label, is_independent=True, is_binned=True, units='')
    x.digits = n_significant_digits
    x.values = h['x_edges']

    y = hepdata_lib.Variable(y_label, is_independent=False, is_binned=False, units='')
    y.digits = n_significant_digits
    y.values = h['y']
    y.add_qualifier('RE', 'P P --> Z jet X')
    y.add_qualifier('SQRT(S)', 8., 'TeV')
    y.add_qualifier('ETARAP', r'2.5 < $y_{\textrm{jet}}$ < 4.0')
    y.add_qualifier('jet radius', jetR)
    y.add_qualifier('jet method', r'Anti-$k_{\textrm{T}}$')
    y.add_qualifier('jet pT', r'%.0f < $p_{\textrm{T,jet}}$ < %.0f GeV' % (min_pt, max_pt))

    # Define uncertainties
    stat = hepdata_lib.Uncertainty('stat', is_symmetric=True)
    stat.values = h['dy'] #float('{:.2g}'.format(dy))

    # Add tables to submission
    table.add_variable(x)
    table.add_variable(y)
    y.add_uncertainty(stat)

    ''' TODO: systematic breakdown
    # Add unfolding systematic
    name = 'hSystematic_Unfolding_R{}_{}_{}-{}'.format(self.utils.remove_periods(jetR), obs_label,
                                                      int(min_pt), int(max_pt))
    h_sys_unfolding = hepdata_reader_systematics.read_hist_1d(getattr(self, name).GetName())
    sys_unfolding = hepdata_lib.Uncertainty('sys,unfolding', is_symmetric=True)
    #sys_unfolding_percent = h_sys_unfolding['y']  # Percent uncertainty
    sys_unfolding.values = [h_sys_unfolding['y'][i] * y.values[i] / 100. for i in range(len(y.values))]
    y.add_uncertainty(sys_unfolding)

    # Add generator systematic
    name = 'hSystematic_generator_R{}_{}_{}-{}'.format(self.utils.remove_periods(jetR), obs_label,
                                                      int(min_pt), int(max_pt))
    h_sys_generator = hepdata_reader_systematics.read_hist_1d(getattr(self, name).GetName())
    sys_generator = hepdata_lib.Uncertainty('sys,generator', is_symmetric=True)
    #sys_generator_percent = h_sys_generator['y']  # Percent uncertainty
    sys_generator.values = [h_sys_generator['y'][i] * y.values[i] / 100. for i in range(len(y.values))]
    y.add_uncertainty(sys_generator)

    # Add systematic uncertainty breakdown
    for systematic in self.systematics_list:

        if systematic in ['main', 'prior1', 'truncation', 'binning'] or 'generator' in systematic:
            continue

        h_sys = self.retrieve_systematic(
          systematic, jetR, obs_label, None, min_pt, max_pt)
        if not h_sys:
            continue

        h_sys = hepdata_reader_systematics.read_hist_1d(h_sys.GetName())
        sys = hepdata_lib.Uncertainty('sys,{}'.format(systematic), is_symmetric=True)
        #sys_percent = h_sys['y']  # Percent uncertainty
        sys.values = [h_sys['y'][i] * y.values[i] / 100. for i in range(len(y.values))]
        y.add_uncertainty(sys)
    

    # Add total systematic
    g_name_totalsys = "sys_%s_%i_%i" % (observable, min_pt, max_pt)
    g_totalsys = hepdata_reader.read_graph(g_name_totalsys)
    totalsys = hepdata_lib.Uncertainty('total sys.', is_symmetric=True)
    totalsys.values = g_totalsys['dy'] #float('{:.2g}'.format(dy))
    y.add_uncertainty(totalsys)
    '''

    # Add table to the submission
    hepdata_submission.add_table(table)


#----------------------------------------------------------------------
def set_hepdata_table_descriptors(table, observable, min_pt, max_pt):

    if observable == 'z':
      if np.isclose(min_pt, 20.):
        table.location = "Figure 5a"
      if np.isclose(min_pt, 30.):
        table.location = "Figure 5b"
      if np.isclose(min_pt, 50.):
        table.location = "Figure 5c"

      table.description = "Distributions of the longitudinal momentum fraction of the hadron with respect to the jet.\n"
      table.description += r'${} < p_{{\textrm{{T,jet}}}} < {}$ GeV'.format(min_pt, max_pt)

      x_label = r'${}$'.format(observable)
      y_label = r'$\frac{1}{N_{\textrm{Z+jet}}} \frac{\textrm{d} N}{\textrm{d}z}$'

    elif observable == "jt":
      if np.isclose(min_pt, 20.):
        table.location = "Figure 6a"
      if np.isclose(min_pt, 30.):
        table.location = "Figure 6b"
      if np.isclose(min_pt, 50.):
        table.location = "Figure 6c"

      table.description = "Distributions of the transverse momentum of charged hadrons with respect to the jet axis.\n"
      table.description += r'${} < p_{{\textrm{{T,jet}}}} < {}$ GeV'.format(min_pt, max_pt)

      x_label = r'${}$'.format(observable)
      y_label = r'$\frac{1}{N_{\textrm{Z+jet}}} \frac{\textrm{d} N}{\textrm{d}j_{\textrm{T}}}[\textrm{$GeV^{-1}$}]$'

    elif observable == "r":
      if np.isclose(min_pt, 20.):
        table.location = "Figure 7a"
      if np.isclose(min_pt, 30.):
        table.location = "Figure 7b"
      if np.isclose(min_pt, 50.):
        table.location = "Figure 7c"

      table.description = "Radial profile distributions of hadrons with respect to the jet axis.\n"
      table.description += r'${} < p_{{\textrm{{T,jet}}}} < {}$ GeV.'.format(min_pt, max_pt)

      x_label = r'${}$'.format(observable)
      y_label = r'$\frac{1}{N_{\textrm{Z+jet}}} \frac{\textrm{d} N}{\textrm{d}r}$'

    else:
      raise ValueError("Observable %s not implemented!" % observable)

    return x_label, y_label


#----------------------------------------------------------------------
# Create submission
hepdata_dir = os.path.join("/home", "alicexty", "rivet", "Rivet", "pythia_paper")
if not os.path.exists(hepdata_dir):
  os.makedirs(hepdata_dir)
hepdata_submission = hepdata_lib.Submission()

# Settings from header file
observables = ["z","jt","r"]                      # Observable name (for internal ref.)
jetR = 0.5                                        # Jet resolution parameter
n_iter = 3                                        # Number of iterations through unfolding
pt_binedgesfinal = [20, 30, 50, 100]  # Reported pT bin edges

# Keep track of count for the HepData table number
i_count = 0

# Loop through WTA subconfigurations
for observable in observables:
  
  # Loop through pt slices, and add table for each result
  for i_pt in range(len(pt_binedgesfinal) - 1):
    min_pt = pt_binedgesfinal[i_pt]
    max_pt = pt_binedgesfinal[i_pt+1]

    add_hepdata_table(i_count, observable, min_pt, max_pt)

    # Update count for next table
    i_count += 1

# Write submission files
hepdata_submission.create_files(hepdata_dir, remove_old=False)
