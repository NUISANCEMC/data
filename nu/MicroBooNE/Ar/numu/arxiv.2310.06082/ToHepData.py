#!/usr/bin/env python3

import yaml
import csv
import os
import ROOT
import re
from math import sqrt

from hepdata_lib import Submission, Table, Variable, Uncertainty, RootFileReader

NUISANCE_DATA_ROOT = os.environ.get("NUISANCE_DATA_ROOT")

if not NUISANCE_DATA_ROOT:
  print("[ERROR]: NUISANCE_DATA_ROOT is not set.")
  exit(1)

probe = "nu"
expt = "MicroBooNE"
target = "Ar"
species = "numu"
ref = "arxiv.2310.06082"

release_dir = "/".join([NUISANCE_DATA_ROOT,"released", probe, expt, target, species, ref])
nuisance_dir = "/".join([NUISANCE_DATA_ROOT, probe, expt, target, species, ref])


def build_flux_table(hname, tname):
  inFileName = "MicroBooNE_FHC_numu_flux.root"

  reader = RootFileReader("/".join([nuisance_dir, inFileName]))
  
  fh = reader.read_hist_1d(hname)

  #### Build Submission
  EnuVar = Variable("e_nu", is_independent=True, is_binned=True, units="GeV")
  EnuVar.values = fh["x_edges"]

  FluxVar = Variable("flux_nu", is_independent=False, is_binned=False, units="$/cm^{2}$")
  FluxVar.values = fh["y"]

  FluxVar.add_qualifier("probe_species", "numu")
  FluxVar.add_qualifier("bin_content_type", "count")

  FluxTable = Table(tname)

  FluxTable.add_variable(EnuVar)
  FluxTable.add_variable(FluxVar)

  return FluxTable

def single_diff_measurement(varname, projname, xsunits, covunits, *args, **kwargs):
  binFileName = "DataRelease.root"
  binFile ="/".join([release_dir, binFileName])
  reader = RootFileReader(binFile)

  inHist = reader.read_hist_1d("TotalUnc_%s" % varname)
  inCov = reader.read_hist_2d("Cov_%s" % varname)
  inAc = reader.read_hist_2d("Ac_%s" % varname)

#### Build Submission
  XVar = Variable(projname, is_independent=True, is_binned=True, units="")
  XVar.values = inHist["x_edges"]

  if "x_pretty_name" in kwargs:
    XVar.add_qualifier("pretty_name", kwargs["x_pretty_name"])

  CrossSection = Variable("cross_section", is_independent=False, is_binned=False, units=xsunits)
  CrossSection.values = inHist["y"]

  CrossSection.add_qualifier("select", "MicroBooNE_CC0Pi_GKI_nu_SelectSignal")
  CrossSection.add_qualifier("project:%s" % projname, "MicroBooNE_CC0Pi_GKI_nu_%s" % projname)
  CrossSection.add_qualifier("target", "Ar")
  CrossSection.add_qualifier("probe_species", "numu")
  CrossSection.add_qualifier("probe_spectrum", "microboone_flux_numu")
  CrossSection.add_qualifier("variable_type", "cross-section-measurement")
  CrossSection.add_qualifier("covariance", "covariance-%s" % projname)
  CrossSection.add_qualifier("smearing", "smearing-%s" % projname)
  if "pretty_name" in kwargs:
    CrossSection.add_qualifier("pretty_name", kwargs["pretty_name"])

  TotalUncertainty = Uncertainty("total", is_symmetric=True)
  TotalUncertainty.values = inHist["dy"]

  CrossSection.add_uncertainty(TotalUncertainty)

  xsTable = Table("cross_section-%s" % projname)
  xsTable.description = ""
  xsTable.location = ""

  xsTable.add_variable(XVar)
  xsTable.add_variable(CrossSection)

  if "observables" in kwargs:
    xsTable.keywords["observables"] = kwargs["observables"]
  if "reactions" in kwargs:
    xsTable.keywords["reactions"] = kwargs["reactions"]
  if "phrases" in kwargs:
    xsTable.keywords["phrases"] = kwargs["phrases"]

##### matrices

  bin_i = Variable("bin_i", is_independent=True, is_binned=False, units="")
  bin_i.values = []

  bin_j = Variable("bin_j", is_independent=True, is_binned=False, units="")
  bin_j.values = []

  for j in range(len(inHist["x"])):
    for i in range(len(inHist["x"])):
      bin_i.values.append(i)
      bin_j.values.append(j)

  Covariance = Variable("covariance", is_independent=False, is_binned=False, units=covunits)
  Covariance.values = inCov["z"]

  WSmear = Variable("wiener_svd-smearing-matrix", is_independent=False, is_binned=False)
  WSmear.values = inAc["z"]

  Covariance.add_qualifier("variable_type", "covariance")
  WSmear.add_qualifier("variable_type", "wiener_svd-smearing-matrix")

  covTable = Table("covariance-%s" % projname)
  covTable.description = ""
  covTable.location = ""

  covTable.add_variable(bin_i)
  covTable.add_variable(bin_j)
  covTable.add_variable(Covariance)

  smearTable = Table("smearing-%s" % projname)
  smearTable.description = ""
  smearTable.location = ""

  smearTable.add_variable(bin_i)
  smearTable.add_variable(bin_j)
  smearTable.add_variable(WSmear)

  return (xsTable,covTable,smearTable)

submission = Submission()
submission.read_abstract("/".join([nuisance_dir,"abstract.txt"]))

pd_tables = single_diff_measurement(varname="DeltaPn", projname="pn", xsunits=r"$\text{cm}^{2}\ c/\text{GeV}\ /\text{Nucleon}$", covunits=r"$(\text{cm}^{2}\ c/\text{GeV}\ /\text{Nucleon})^{2}$")
submission.add_table(pd_tables[0])
submission.add_table(pd_tables[1])
submission.add_table(pd_tables[2])

pd_tables = single_diff_measurement(varname="DeltaAlpha3Dq", projname="alpha3d", xsunits=r"degrees", covunits=r"$\text{degrees}^{2}$")
submission.add_table(pd_tables[0])
submission.add_table(pd_tables[1])
submission.add_table(pd_tables[2])

pd_tables = single_diff_measurement(varname="DeltaPhi3D", projname="phi3d", xsunits=r"degrees", covunits=r"$\text{degrees}^{2}$")
submission.add_table(pd_tables[0])
submission.add_table(pd_tables[1])
submission.add_table(pd_tables[2])

pd_tables = single_diff_measurement(varname="DeltaPnPar", projname="pn_para", xsunits=r"$\text{cm}^{2}\ c/\text{GeV}\ /\text{Nucleon}$", covunits=r"$(\text{cm}^{2}\ c/\text{GeV}\ /\text{Nucleon})^{2}$")
submission.add_table(pd_tables[0])
submission.add_table(pd_tables[1])
submission.add_table(pd_tables[2])

pd_tables = single_diff_measurement(varname="DeltaPnPerp", projname="pn_perp", xsunits=r"$\text{cm}^{2}\ c/\text{GeV}\ /\text{Nucleon}$", covunits=r"$(\text{cm}^{2}\ c/\text{GeV}\ /\text{Nucleon})^{2}$")
submission.add_table(pd_tables[0])
submission.add_table(pd_tables[1])
submission.add_table(pd_tables[2])

pd_tables = single_diff_measurement(varname="DeltaPnPerpx", projname="pn_perp_x", xsunits=r"$\text{cm}^{2}\ c/\text{GeV}\ /\text{Nucleon}$", covunits=r"$(\text{cm}^{2}\ c/\text{GeV}\ /\text{Nucleon})^{2}$")
submission.add_table(pd_tables[0])
submission.add_table(pd_tables[1])
submission.add_table(pd_tables[2])

pd_tables = single_diff_measurement(varname="DeltaPnPerpy", projname="pn_perp_y", xsunits=r"$\text{cm}^{2}\ c/\text{GeV}\ /\text{Nucleon}$", covunits=r"$(\text{cm}^{2}\ c/\text{GeV}\ /\text{Nucleon})^{2}$")
submission.add_table(pd_tables[0])
submission.add_table(pd_tables[1])
submission.add_table(pd_tables[2])

submission.add_table(build_flux_table("numu_hist", "microboone_flux_numu"))

submission.add_additional_resource(description="Python conversion script used to build this submisson. Part of NUISANCE.",
    location="/".join([nuisance_dir,"ToHepData.py"]),
    copy_file=True, file_type="NUISANCE")

submission.add_additional_resource(description="ROOT version of the flux provided by the MicroBooNE collaboration.",
    location="/".join([nuisance_dir,"MicroBooNE_FHC_numu_flux.root"]),
    copy_file=True)

submission.add_additional_resource(description="Selection and projection function examples. Can be executued in the ProSelecta environment v1.0.",
    location="/".join([nuisance_dir,"analysis.cxx"]),
    copy_file=True)

submission.add_additional_resource(description="2D Binning scheme",
    location="/".join([release_dir,"BinScheme.txt"]),
    copy_file=True)
submission.add_additional_resource(description="Official data release documentation",
    location="/".join([release_dir,"README.txt"]),
    copy_file=True)

submission.add_link(description="pre-print", location="https://doi.org/10.48550/arXiv.2310.06082")
submission.add_link(description="use with NUISANCE3", location="https://github.com/NUISANCEMC/nuisance3")

submission.create_files("testout", remove_old=True)
