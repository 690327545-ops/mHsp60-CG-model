#!/bin/bash

function get_resname {
echo $1 |
  sed "s/GLY/GA/g" |
  sed "s/ALA/AA/g" |
  sed "s/VAL/VA/g" |
  sed "s/LEU/LA/g" |
  sed "s/NLE/LB/g" |
  sed "s/ILE/IA/g" |
  sed "s/PRO/PA/g" |
  sed "s/HYP/PB/g" |
  sed "s/PHE/FA/g" |
  sed "s/TYR/YA/g" |
  sed "s/LYN/YB/g" |
  sed "s/TRP/WA/g" |
  sed "s/SER/SA/g" |
  sed "s/THR/TA/g" |
  sed "s/CYS/CA/g" |
  sed "s/CYM/CB/g" |
  sed "s/CYX/CC/g" |
  sed "s/MET/MA/g" |
  sed "s/ASN/NA/g" |
  sed "s/GLN/OA/g" |
  sed "s/ASP/DA/g" |
  sed "s/ASH/DB/g" |
  sed "s/GLU/EA/g" |
  sed "s/GLH/EB/g" |
  sed "s/LYS/KA/g" |
  sed "s/ARG/RA/g" |
  sed "s/HIS/HA/g" |
  sed "s/HIE/HB/g" |
  sed "s/HID/HC/g" |
  sed "s/HIP/HD/g"
}

prot_res_aa="ALA ARG ASH ASN ASP CALA CARG CASN CASP CCYS CCYX CGLN CGLU CGLY CHID CHIE CHIP CHYP CILE CLEU CLYS CMET CPHE CPRO CSER CTHR CTYR CVAL CYM CYS CYX GLH GLN GLU GLY HID HIE HIP HYP ILE LEU LYN LYS MET NALA NARG NASN NASP NCYS NCYX NGLN NGLU NGLY NHE NHID NHIE NHIP NILE NLEU NLYS NMET NPHE NPRO NSER NTHR NTRP NTYR NVAL PHE PRO SER THR TRP TYR VAL"
prot_res_ua=""
prot_res_ma=""

for f in $prot_res_aa; do
  prot_res_ua="$prot_res_ua""U`get_resname $f` "
  prot_res_ma="$prot_res_ma""M`get_resname $f` "
done

echo resname $prot_res_aa
echo resname $prot_res_ua
echo resname $prot_res_ma

