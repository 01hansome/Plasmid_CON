# -*- coding: utf-8 -*-
# File : plasmid.PY
# Author: HANSOME
# Date : 2021/11/22

from Bio.Seq import Seq

def plasmidMCal(an, cn, gn, tn, concentration):
    mo_weight = an * 313.2 + tn * 304.2 + cn * 289.2 + gn * 329.2
    print("Plasmid_M", mo_weight)
    copies = (concentration * 6.022 * (10 ** 14)) / mo_weight
    print("Copies(copies/μL): %.3e" % copies)

def get_dNmpNum(positive_dna_strand):
    n_num = {}
    negative_dna_strand = positive_dna_strand.reverse_complement()
    n_num["A"] = positive_dna_strand.count("A") + negative_dna_strand.count("A")
    n_num["C"] = positive_dna_strand.count("C") + negative_dna_strand.count("C")
    n_num["G"] = positive_dna_strand.count("G") + negative_dna_strand.count("G")
    n_num["T"] = positive_dna_strand.count("T") + negative_dna_strand.count("T")
    return n_num
def getDNAParameter():
    return {"one_plasmid_strand": input("Input one plasmid strand:").upper(),
            "rna_concentration": float(input("plasmid concentration (ng/μL)："))}

if __name__ == '__main__':
    concentration = float(input("Concentration(ng/μL)："))
    single_dna = Seq(input("Input P/N Plasmid strand ：").upper())
    dntp = get_dNmpNum(single_dna)
    plasmidMCal(dntp["A"], dntp["C"], dntp["G"], dntp["T"], concentration )

