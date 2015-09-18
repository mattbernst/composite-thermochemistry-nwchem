#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
from BeautifulSoup import BeautifulSoup as BS
import csv
import os
import urllib
import sys

prefix = "test_data/G2-97/"
gm = {"LiH" : ("lithiumhydride.xyz", 0, "singlet"),
      "BeH" : ("berylliumhydride.xyz", 0, "doublet"),
      "CH" : ("ch_rad.xyz", 0, "doublet"),
      "CH2(3B1)" : ("methylene_triplet.xyz", 0, "triplet"),
      "CH2(1A1)" : ("methylene_singlet.xyz", 0, "singlet"),
      "CH3" : ("methyl_rad.xyz", 0, "doublet"),
      "CH4" : ("methane.xyz", 0, "singlet"),
      "NH" : ("nh_rad.xyz", 0, "triplet"),
      "NH2" : ("nh2_rad.xyz", 0, "doublet"),
      "NH3" : ("ammonia.xyz", 0, "singlet"),
      "OH" : ("oh_rad.xyz", 0, "doublet"),
      "OH2" : ("water.xyz", 0, "singlet"),
      "FH" : ("hydrogenfluoride.xyz", 0, "singlet"),
      "SiH2(1A1)" : ("silylene_singlet.xyz", 0, "singlet"),
      "SiH2(3B1)" : ("silylene_triplet.xyz", 0, "triplet"),
      "SiH3" : ("silyl_rad.xyz", 0, "doublet"),
      "SiH4" : ("silane.xyz", 0, "singlet"),
      "PH2" : ("ph2_rad.xyz", 0, "doublet"),
      "PH3" : ("phosphine.xyz", 0, "singlet"),
      "SH2" : ("hydrogensulfide.xyz", 0, "singlet"),
      "ClH" : ("hydrogenchloride.xyz", 0, "singlet"),
      "Li2" : ("li2.xyz", 0, "singlet"),
      "LiF" : ("lithiumfluoride.xyz", 0, "singlet"),
      "C2H2" : ("acetylene.xyz", 0, "singlet"),
      "C2H4" : ("ethylene.xyz", 0, "singlet"),
      "C2H6" : ("ethane.xyz", 0, "singlet"),
      "CN" : ("cyano_rad.xyz", 0, "doublet"),
      "HCN" : ("hydrogencyanide.xyz", 0, "singlet"),
      "CO" : ("carbonmonoxide.xyz", 0, "singlet"),
      "HCO" : ("hco_rad.xyz", 0, "doublet"),
      "H2CO" : ("formaldehyde.xyz", 0, "singlet"),
      "H3COH" : ("methylalcohol.xyz", 0, "singlet"),
      "N2" : ("n2.xyz", "singlet"),
      "H2NNH2" : ("hydrazine.xyz", 0, "singlet"),
      "NO" : ("nitricoxide.xyz", 0, "doublet"),
      "O2" : ("o2.xyz", 0, "triplet"),
      "HOOH" : ("hydrogenperoxide.xyz", 0, "singlet"),
      "F2" : ("f2.xyz", "singlet"),
      "CO2" : ("carbondioxide.xyz", 0, "singlet"),
      "Na2" : ("na2.xyz", 0, "singlet"),
      "Si2" : ("si2.xyz", 0, "triplet"),
      "P2" : ("p2.xyz", 0, "singlet"),
      "S2" : ("s2.xyz", 0, "triplet"),
      "Cl2" : ("cl2.xyz", 0, "singlet"),
      "NaCl" : ("sodiumchloride.xyz", 0, "singlet"),
      "SiO" : ("siliconmonoxide.xyz", 0, "singlet"),
      "SC" : (),
      "SO" : ("sulfurmonoxide.xyz", 0, "triplet"),
      "ClO" : ("chlorinemonoxide.xyz", 0, "doublet"),
      "FCl" : ("chlorinemonofluoride.xyz", 0, "singlet"),
      "Si2H6" : ("disilane.xyz", 0, "singlet"),
      "CH3Cl" : ("methylchloride.xyz", 0, "singlet"),
      "H3CSH" : ("methanethiol.xyz", 0, "singlet"),
      "HOCl" : ("hypochlorousacid.xyz", 0, "singlet"),
      "SO2" : ("sulfurdioxide.xyz", 0, "singlet"),
      "BF3" : ("trifluoro-borane.xyz", 0, "singlet"),
      "BCl3" : ("trichloro-borane.xyz", 0, "singlet"),
      "AlF3" : ("aluminumtrifluoride.xyz", 0, "singlet"),
      "AlCl3" : ("aluminumtrichloride.xyz", 0, "singlet"),
      "CF4" : ("carbontetrafluoride.xyz", 0, "singlet"),
      "CCl4" : ("carbontetrachloride.xyz", 0, "singlet"),
      "COS" : ("carbonylsulfide.xyz", 0, "singlet"),
      "CS2" : ("carbondisulfide.xyz", 0, "singlet"),
      "COF2" : ("carbonicdifluoride.xyz", 0, "singlet"),
      "SiF4" : ("silicontetrafluoride.xyz", 0, "singlet"),
      "SiCl4" : ("tetrachloro-silane.xyz", 0, "singlet"),
      "N2O" : ("nitrousoxide.xyz", 0, "singlet"),
      "ClNO" : ("nitrosylchloride.xyz", 0, "singlet"),
      "NF3" : ("nitrogentrifluoride.xyz", 0, "singlet"),
      "PF3" : ("phosphorustrifluoride.xyz", 0, "singlet"),
      "O3" : ("ozone.xyz", 0, "singlet"),
      "F2O" : ("difluorinemonoxide.xyz", 0, "singlet"),
      "ClF3" : ("chlorinetrifluoride.xyz", 0, "singlet"),
      "C2F4" : ("tetrafluoroethylene.xyz", 0, "singlet"),
      "C2Cl4" : ("tetrachloroethylene.xyz", 0, "singlet"),
      "CF3CN" : ("trifluoro-acetonitrile.xyz", 0, "singlet"),
      "CH3CCH (propyne)" : ("propyne.xyz", 0, "singlet"),
      "CH2=C=CH2(allene)" : ("allene.xyz", 0, "singlet"),
      "C3H4(cyclopropene)" : ("cyclopropene.xyz", 0, "singlet"),
      "CH3CH=CH2(propylene)" : ("propene.xyz", 0, "singlet"),
      "C3H6(cyclopropane)" : ("cyclopropane.xyz", 0, "singlet"),
      "C3H8(propane)" : ("propane.xyz", 0, "singlet"),
      "CH2CHCHCH2(butadiene)" : ("1,3-butadiene.xyz", 0, "singlet"),
      "C4H6(2-butyne)" : ("2-butyne.xyz", 0, "singlet"),
      "C4H6(methylene cyclopropane)" : ("methylenecyclopropane.xyz", 0, "singlet"),
      "C4H6(bicyclobutane)" : ("bicyclobutane.xyz", 0, "singlet"),
      "C4H6(cyclobutene)" : ("cyclobutene.xyz", 0, "singlet"),
      "C4H8(cyclobutane)" : ("cyclobutane.xyz", 0, "singlet"),
      "C4H8(isobutene)" : (),
      "C4H10(trans butane)" : (),
      "C4H10(isobutane)" : ("isobutane.xyz", 0, "singlet"),
      "C5H8(spiropentane)" : ("spiropentane.xyz", 0, "singlet"),
      "C6H6(benzene)" : ("benzene.xyz", 0, "singlet"),
      "CH2F2" : ("difluoromethane.xyz", 0, "singlet"),
      "CHF3" : ("trifluoro-methane.xyz", 0,"singlet"),
      "CH2Cl2" : ("methylenechloride.xyz", 0, "singlet"),
      "CHCl3" : ("chloroform.xyz", 0, "singlet"),
      "CH3NH2(methylamine)" : ("methylamine.xyz", 0, "singlet"),
      "CH3CN (methyl cyanide)" : ("acetonitrile.xyz", 0, "singlet"),
      "CH3NO2(nitromethane)" : ("nitro-methane.xyz", 0, "singlet"),
      "CH3ONO (methyl nitrite)" : ("methylnitrite.xyz", 0, "singlet"),
      "CH3SiH3(methyl silane)" : ("methylsilane.xyz", 0, "singlet"),
      "HCOOH (formic acid)" : ("formicacid.xyz", 0, "singlet"),
      "HCOOCH3(methyl formate)" : ("methylformate.xyz", 0, "singlet"),
      "CH3CONH2(acetamide)" : ("acetamide.xyz", 0, "singlet"),
      "C2H4NH (aziridine)" : ("aziridine.xyz", 0, "singlet"),
      "NCCN (cyanogen)" : ("cyanogen.xyz", 0, "singlet"),
      "(CH3)2NH (dimethylamine)" : ("dimethylamine.xyz", 0, "singlet"),
      "CH3CH2NH2(trans ethylamine)" : ("ethylamine.xyz", 0, "singlet"),
      "CH2CO (ketene)" : ("ketene.xyz", 0, "singlet"),
      "C2H4O (oxirane)" : ("ethyleneoxide.xyz", 0, "singlet"),
      "CH3CHO (acetaldehyde)" : ("acetaldehyde.xyz", 0, "singlet"),
      "HCOCOH (glyoxal)" : ("ethanedial.xyz", 0, "singlet"),
      "CH3CH2OH (ethanol)" : ("ethanol.xyz", 0, "singlet"),
      "CH3OCH3(dimethylether)" : ("dimethylether.xyz", 0, "singlet"),
      "C2H4S (thiirane)" : ("thiirane.xyz", 0, "singlet"),
      "(CH3)2SO (dimethyl sulfoxide)" : ("dimethylsulfoxide.xyz", 0, "singlet"),
      "C2H5SH (ethanethiol)" : ("ethanethiol.xyz", 0, "singlet"),
      "CH3SCH3(dimethyl sulfide)" : ("dimethylsulfide.xyz", 0, "singlet"),
      "CH2=CHF (vinyl fluoride)" : ("fluoro-ethene.xyz", 0, "singlet"),
      "C2H5Cl (ethyl chloride)" : ("ethylchloride.xyz", 0, "singlet"),
      "CH2=CHCl (vinyl chloride)" : ("chloro-ethene.xyz", 0, "singlet"),
      "CH2=CHCN (acrylonitrile)" : ("acrylonitrile.xyz", 0, "singlet"),
      "CH3COCH3(acetone)" : ("acetone.xyz", 0, "singlet"),
      "CH3COOH (acetic acid)" : ("aceticacid.xyz", 0, "singlet"),
      "CH3COF (acetyl fluoride)" : ("acetylfluoride.xyz", 0, "singlet"),
      "CH3COCl (acetyl chloride)" : ("acetylchloride.xyz", 0, "singlet"),
      "CH3CH2CH2Cl (propyl chloride)" : ("1-chloro-propane.xyz", 0, "singlet"),
      "(CH3)2CHOH (isopropanol)" : ("isopropanol.xyz", 0, "singlet"),
      "C2H5OCH3(methyl ethyl ether)" : ("methoxy-ethane.xyz", 0, "singlet"),
      "(CH3)3N (trimethylamine)" : ("trimethylamine.xyz", 0, "singlet"),
      "C4H4O (furan)" : ("furan.xyz", 0, "singlet"),
      "C4H4S (thiophene)" : ("thiophene.xyz", 0, "singlet"),
      "C4H5N (pyrrole)" : ("pyrrole.xyz", 0, "singlet"),
      "C5H5N (pyridine)" : ("pyridine.xyz", 0, "singlet"),
      "H2" : ("h2.xyz", 0, "singlet"),
      "HS" : ("sh_rad.xyz", 0, "doublet"),
      "CCH" : ("cch_rad.xyz", 0, "doublet"),
      "C2H3(2A')" : ("c2h3_rad.xyz", 0, "doublet"),
      "CH3CO (2A')" : ("ch3co_rad.xyz", 0, "doublet"),
      "H2COH (2A)" : ("h2coh_rad.xyz", 0, "doublet"),
      "CH3O CS (2A')" : ("ch3o_rad.xyz", 0, "doublet"),
      "CH3CH2O (2A'')" : ("ch3ch2o_rad.xyz", 0, "doublet"),
      "CH3S (2A')" : ("ch3s_rad.xyz", 0, "doublet"),
      "C2H5(2A')" : ("c2h5_rad.xyz", 0, "doublet"),
      "(CH3)2CH (2A')" : ("isopropyl_rad.xyz", 0, "doublet"),
      "(CH3)3C (t-butyl radical)" : ("t-butyl_rad.xyz", 0, "doublet"),
      "NO2" : ("no2.xyz", 0, "doublet")
      }

def main():
    urls = ["http://www.cse.anl.gov/OldCHMwebsiteContent/compmat/g3energies/g3neut.htm",
            "http://www.cse.anl.gov/OldCHMwebsiteContent/compmat/g3energies/G3_MP3_DHF.htm",
            "http://www.cse.anl.gov/OldCHMwebsiteContent/compmat/g3energies/g3mp2neut.htm",
            "http://www.cse.anl.gov/OldCHMwebsiteContent/compmat/g3energies/G3_B3lyp_DHF.htm",
            "http://www.cse.anl.gov/OldCHMwebsiteContent/compmat/g3energies/G3_mp2_b3lyp_DHF.htm"]
    for url in urls:
        url2csv(url)

def url2csv(url):
    rows = []
    header = None
    htmlfile = url.split("/")[-1]
    if not os.path.exists(htmlfile):
        data = urllib.urlopen(url).read()
        with open(htmlfile, "w") as outfile:
            outfile.write(data)
    else:
        with open(htmlfile) as infile:
            data = infile.read()

    soup = BS(data, convertEntities=BS.HTML_ENTITIES)
    for tr in soup.findAll("tr"):
        row = [c.text for c in tr.findAll("td")]
        if not header:
            header = row[:1] + ["System", "Charge", "Multiplicity"] + row[1:]
        else:
            try:
                s, c, m = gm[row[0]]
                s = prefix + s
            except ValueError:
                s, c, m = "none", 0, "singlet"
            r = dict(zip(header, row[:1] + [s, c, m] + row[1:]))
            rows.append(r)

    csvfile = htmlfile.replace(".htm", ".csv")
    csvwrite(csvfile, rows, header)

def csvwrite(outname, rows, header):
    with open(outname, "w") as outfile:
        d = csv.DictWriter(outfile, fieldnames=header)
        d.writeheader()
        for row in rows:
            d.writerow(row)

if __name__ == "__main__":
    main()
