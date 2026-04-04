# For documentation see:
# https://gitlab.com/hepcedar/rivet/blob/release-4-1-x/doc/tutorials/makeplots.md

####################################################################################
# Longitudinal momentum fraction

BEGIN PLOT /LHCB_2019_I1730448/d01-x01-y01
Title=Distributions of the longitudinal momentum fraction of the hadron with respect to the jet.
XLabel=${z}$
YLabel=$\frac{1}{N_{\textrm{Z+jet}}} \frac{\textrm{d} N}{\textrm{d}z}$
YLabelSep=5.5
ErrorBars=1
NormalizeToIntegral=1
LogY=1
#YMin=0
#YMax=200
LogX=1
#XMin=0
#XMax=0.7
Legend=1
CustomLegend=${20} < p_{{\textrm{{T,jet}}}} < {30}$ GeV
#LegendXPos=0.65
END PLOT

BEGIN PLOT /LHCB_2019_I1730448/d01-x01-y02
Title=Distributions of the longitudinal momentum fraction of the hadron with respect to the jet.
XLabel=${z}$
YLabel=$\frac{1}{N_{\textrm{Z+jet}}} \frac{\textrm{d} N}{\textrm{d}z}$
YLabelSep=5.5
ErrorBars=1
NormalizeToIntegral=1
LogY=1
#YMin=0
#YMax=5
LogX=1
#XMin=0
#XMax=1
Legend=1
CustomLegend=${30} < p_{{\textrm{{T,jet}}}} < {50}$ GeV
#LegendXPos=0.65
END PLOT

BEGIN PLOT /LHCB_2019_I1730448/d01-x01-y03
Title=Distributions of the longitudinal momentum fraction of the hadron with respect to the jet.
XLabel=${z}$
YLabel=$\frac{1}{N_{\textrm{Z+jet}}} \frac{\textrm{d} N}{\textrm{d}z}$
YLabelSep=5.5
ErrorBars=1
NormalizeToIntegral=1
LogY=1
#YMin=0
#YMax=5
LogX=1
#XMin=0
#XMax=1
Legend=1
CustomLegend=${50} < p_{{\textrm{{T,jet}}}} < {100}$ GeV
#LegendXPos=0.65
END PLOT

####################################################################################
# Transverse momentum

BEGIN PLOT /LHCB_2019_I1730448/d02-x01-y01
Title=Distributions of the transverse momentum of charged hadrons with respect to the jet axis.
XLabel=${jt}$
YLabel=$\frac{1}{N_{\textrm{Z+jet}}} \frac{\textrm{d} N}{\textrm{d}j_{\textrm{T}}}[\textrm{$GeV^{-1}$}]$
YLabelSep=5.5
ErrorBars=1
NormalizeToIntegral=1
LogY=1
#YMin=0
#YMax=5
LogX=0
XMin=0
XMax=3
Legend=1
CustomLegend=${20} < p_{{\textrm{{T,jet}}}} < {30}$ GeV
#LegendXPos=0.65
END PLOT

BEGIN PLOT /LHCB_2019_I1730448/d02-x01-y02
Title=Distributions of the transverse momentum of charged hadrons with respect to the jet axis.
XLabel=${jt}$
YLabel=$\frac{1}{N_{\textrm{Z+jet}}} \frac{\textrm{d} N}{\textrm{d}j_{\textrm{T}}}[\textrm{$GeV^{-1}$}]$
YLabelSep=5.5
ErrorBars=1
NormalizeToIntegral=1
LogY=1
#YMin=0
#YMax=5
LogX=0
XMin=0
XMax=3
Legend=1
CustomLegend=${30} < p_{{\textrm{{T,jet}}}} < {50}$ GeV
#LegendXPos=0.65
END PLOT

BEGIN PLOT /LHCB_2019_I1730448/d02-x01-y03
Title=Distributions of the transverse momentum of charged hadrons with respect to the jet axis.
XLabel=${jt}$
YLabel=$\frac{1}{N_{\textrm{Z+jet}}} \frac{\textrm{d} N}{\textrm{d}j_{\textrm{T}}}[\textrm{$GeV^{-1}$}]$
YLabelSep=5.5
ErrorBars=1
NormalizeToIntegral=1
LogY=1
#YMin=0
#YMax=5
LogX=0
XMin=0
XMax=3
Legend=1
CustomLegend=${50} < p_{{\textrm{{T,jet}}}} < {100}$ GeV
#LegendXPos=0.65
END PLOT

####################################################################################
# Radial profile distributions

BEGIN PLOT /LHCB_2019_I1730448/d03-x01-y01
Title=Radial profile distributions of hadrons with respect to the jet axis.
XLabel=${r}$
YLabel=$\frac{1}{N_{\textrm{Z+jet}}} \frac{\textrm{d} N}{\textrm{d}r}$
YLabelSep=5.5
ErrorBars=1
NormalizeToIntegral=1
LogY=1
#YMin=0
#YMax=5
LogX=0
XMin=0
XMax=0.5
Legend=1
CustomLegend=${20} < p_{{\textrm{{T,jet}}}} < {30}$ GeV
#LegendXPos=0.65
END PLOT

BEGIN PLOT /LHCB_2019_I1730448/d03-x01-y02
Title=Radial profile distributions of hadrons with respect to the jet axis.
XLabel=${r}$
YLabel=$\frac{1}{N_{\textrm{Z+jet}}} \frac{\textrm{d} N}{\textrm{d}r}$
YLabelSep=5.5
ErrorBars=1
NormalizeToIntegral=1
LogY=1
#YMin=0
#YMax=5
LogX=0
XMin=0
XMax=0.5
Legend=1
CustomLegend=${30} < p_{{\textrm{{T,jet}}}} < {50}$ GeV
#LegendXPos=0.65
END PLOT

BEGIN PLOT /LHCB_2019_I1730448/d03-x01-y03
Title=Radial profile distributions of hadrons with respect to the jet axis.
XLabel=${r}$
YLabel=$\frac{1}{N_{\textrm{Z+jet}}} \frac{\textrm{d} N}{\textrm{d}r}$
YLabelSep=5.5
ErrorBars=1
NormalizeToIntegral=1
LogY=1
#YMin=0
#YMax=5
LogX=0
XMin=0
XMax=0.5
Legend=1
CustomLegend=${50} < p_{{\textrm{{T,jet}}}} < {100}$ GeV
#LegendXPos=0.65
END PLOT