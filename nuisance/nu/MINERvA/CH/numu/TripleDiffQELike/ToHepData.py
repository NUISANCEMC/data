#!/usr/bin/env python3

import yaml
import csv
import os
import ROOT

from hepdata_lib import Submission, Table, Variable, Uncertainty

NUISANCE_DATA_ROOT = os.environ.get("NUISANCE_DATA_ROOT")

if not NUISANCE_DATA_ROOT:
  print("[ERROR]: NUISANCE_DATA_ROOT is not set.")
  exit(1)

probe = "nu"
expt = "MINERvA"
target = "CH"
species = "numu"
ref = "arXiv-2203-08022v1"

release_dir = "/".join([NUISANCE_DATA_ROOT,"released", probe, expt, target, species, ref])

def sumtpptpz_xsec():
  inFileName = "MINERvA_TripleDiffQELike_ptpzsumtp.root"
  inFile = ROOT.TFile.Open("/".join([release_dir,inFileName]) ,"READ")

  ptpzsumtp_data_cross_section_with_total_unc = inFile.Get("ptpzsumtp_data_cross_section_with_total_unc")
  ptpzsumtp_data_cross_section_with_stat_unc = inFile.Get("ptpzsumtp_data_cross_section_with_stat_unc")

  data = {
    "values": [],
    "staterror": [],
    "toterror": []
  }

  SumTpAxis = []
  MuonPtAxis = []
  MuonPpAxis = []

  for x in range(ptpzsumtp_data_cross_section_with_total_unc.GetXaxis().GetNbins()):
    xl = ptpzsumtp_data_cross_section_with_total_unc.GetXaxis().GetBinLowEdge(x+1)
    xu = ptpzsumtp_data_cross_section_with_total_unc.GetXaxis().GetBinUpEdge(x+1)

    for y in range(ptpzsumtp_data_cross_section_with_total_unc.GetYaxis().GetNbins()):
      yl = ptpzsumtp_data_cross_section_with_total_unc.GetYaxis().GetBinLowEdge(y+1)
      yu = ptpzsumtp_data_cross_section_with_total_unc.GetYaxis().GetBinUpEdge(y+1)

      for z in range(ptpzsumtp_data_cross_section_with_total_unc.GetZaxis().GetNbins()):
        zl = ptpzsumtp_data_cross_section_with_total_unc.GetZaxis().GetBinLowEdge(z+1)
        zu = ptpzsumtp_data_cross_section_with_total_unc.GetZaxis().GetBinUpEdge(z+1)

        bv = (xu-xl) * (yu - yl) * (zu - zl)
        data["values"].append(ptpzsumtp_data_cross_section_with_total_unc.GetBinContent(x+1,y+1,z+1)/bv)
        data["staterror"].append(ptpzsumtp_data_cross_section_with_total_unc.GetBinError(x+1,y+1,z+1)/bv)
        data["toterror"].append(ptpzsumtp_data_cross_section_with_stat_unc.GetBinError(x+1,y+1,z+1)/bv)
        SumTpAxis.append((xl, xu))
        MuonPtAxis.append((yl, yu))
        MuonPpAxis.append((zl, zu))

  #### Build Submission
  MuonPt = Variable("MuonPt", is_independent=True, is_binned=True, units="GeV/c")
  MuonPt.values = MuonPtAxis
  MuonPp = Variable("MuonPp", is_independent=True, is_binned=True, units="GeV/c")
  MuonPp.values = MuonPpAxis
  SumTp = Variable("SumTp", is_independent=True, is_binned=True, units="GeV")
  SumTp.values = SumTpAxis

  StatisticalUncertainty = Uncertainty("Statistical Uncertainty", is_symmetric=True)
  StatisticalUncertainty.values = data["staterror"]
  TotalUncertainty = Uncertainty("Total Uncertainty", is_symmetric=True)
  TotalUncertainty.values = data["toterror"]

  CrossSection = Variable("CrossSection", is_independent=False, is_binned=False, units="cm${}^2$/nucleon/(GeV${}^3$/c${}^2$)")
  CrossSection.values = data["values"]
  CrossSection.add_uncertainty(StatisticalUncertainty)
  CrossSection.add_uncertainty(TotalUncertainty)

  table = Table("MINERvA_TripleDiffQELike_sumtpptpz")
  table.description = "MINERvA_TripleDiffQELike_ptpzsumtp"
  table.location = "MINERvA_TripleDiffQELike_ptpzsumtp"

  table.add_variable(SumTp)
  table.add_variable(MuonPt)
  table.add_variable(MuonPp)
  table.add_variable(CrossSection)

  return table

def sumtpq0qeemu_xsec():
  inFileName = "MINERvA_TripleDiffQELike_q0qeemusumtp.root"
  inFile = ROOT.TFile.Open("/".join([release_dir,inFileName]) ,"READ")

  q0qeemusumtp_data_cross_section_with_total_unc = inFile.Get("q0qeemusumtp_data_cross_section_with_total_unc")
  q0qeemusumtp_data_cross_section_with_stat_unc = inFile.Get("q0qeemusumtp_data_cross_section_with_stat_unc")

  data = {
    "values": [],
    "staterror": [],
    "toterror": []
  }

  SumTpAxis = []
  q0QEAxis = []
  MuonEAxis = []

  for x in range(q0qeemusumtp_data_cross_section_with_total_unc.GetXaxis().GetNbins()):
    xl = q0qeemusumtp_data_cross_section_with_total_unc.GetXaxis().GetBinLowEdge(x+1)
    xu = q0qeemusumtp_data_cross_section_with_total_unc.GetXaxis().GetBinUpEdge(x+1)

    for y in range(q0qeemusumtp_data_cross_section_with_total_unc.GetYaxis().GetNbins()):
      yl = q0qeemusumtp_data_cross_section_with_total_unc.GetYaxis().GetBinLowEdge(y+1)
      yu = q0qeemusumtp_data_cross_section_with_total_unc.GetYaxis().GetBinUpEdge(y+1)

      for z in range(q0qeemusumtp_data_cross_section_with_total_unc.GetZaxis().GetNbins()):
        zl = q0qeemusumtp_data_cross_section_with_total_unc.GetZaxis().GetBinLowEdge(z+1)
        zu = q0qeemusumtp_data_cross_section_with_total_unc.GetZaxis().GetBinUpEdge(z+1)

        bv = (xu-xl) * (yu - yl) * (zu - zl)
        data["values"].append(q0qeemusumtp_data_cross_section_with_total_unc.GetBinContent(x+1,y+1,z+1)/bv)
        data["staterror"].append(q0qeemusumtp_data_cross_section_with_total_unc.GetBinError(x+1,y+1,z+1)/bv)
        data["toterror"].append(q0qeemusumtp_data_cross_section_with_stat_unc.GetBinError(x+1,y+1,z+1)/bv)
        SumTpAxis.append((xl, xu))
        q0QEAxis.append((yl, yu))
        MuonEAxis.append((zl, zu))

  #### Build Submission
  q0QE = Variable("q0QE", is_independent=True, is_binned=True, units="GeV")
  q0QE.values = q0QEAxis
  MuonE = Variable("MuonE", is_independent=True, is_binned=True, units="GeV")
  MuonE.values = MuonEAxis
  SumTp = Variable("SumTp", is_independent=True, is_binned=True, units="GeV")
  SumTp.values = SumTpAxis

  StatisticalUncertainty = Uncertainty("Statistical Uncertainty", is_symmetric=True)
  StatisticalUncertainty.values = data["staterror"]
  TotalUncertainty = Uncertainty("Total Uncertainty", is_symmetric=True)
  TotalUncertainty.values = data["toterror"]

  CrossSection = Variable("CrossSection", is_independent=False, is_binned=False, units="cm${}^2$/nucleon/(GeV${}^3$)")
  CrossSection.values = data["values"]
  CrossSection.add_uncertainty(StatisticalUncertainty)
  CrossSection.add_uncertainty(TotalUncertainty)

  table = Table("MINERvA_TripleDiffQELike_sumtpq0qeemu")
  table.description = "MINERvA_TripleDiffQELike_ptpzsumtp"
  table.location = "MINERvA_TripleDiffQELike_ptpzsumtp"

  table.add_variable(SumTp)
  table.add_variable(q0QE)
  table.add_variable(MuonE)
  table.add_variable(CrossSection)

  return table

submission = Submission()
submission.add_table(sumtpq0qeemu_xsec())
submission.add_table(sumtpptpz_xsec())
submission.create_files("testout", remove_old=True)
