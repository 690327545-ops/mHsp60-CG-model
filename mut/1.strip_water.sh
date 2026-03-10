#!/bin/bash

MA_PATH=$XCG_PATH/templates/martini-harmonic
UA_PATH=$XCG_PATH/templates/united-14sb-lfy
AA_PATH=$XCG_PATH/templates/ff14sb
LI_PATH=`pwd`/top_ligands

#active_resname="ATA ATB ATC ATD ATE ATF ATG KAA KBB KCC KDD KEE KFF KGG MGA MGB MGC MGD MGE MGF MGG"
active_resname="ATC KAC MGC"

function get_resid {
  resname=$1
  pdbfile=$2
  grep ATOM $pdbfile | grep $resname | awk '{print $5}' | uniq
}

function strip_solvents_outside {

  mkdir -p top; cd top
  
  ln -sf ../top_orig/*.pdb    orig.pdb
  ln -sf ../top_orig/*.prmtop orig.prmtop
  ln -sf ../top_orig/*.rst7   orig.rst7

  rad="$1"
  cat << EOF > cpptraj.in
parm orig.prmtop
trajin orig.rst7
reference orig.rst7
strip "(:WAT,Na+) & (:${active_resname// /,} >: $rad)" parmout prot_watsph.prmtop
trajout prot_watsph.rst7
trajout prot_watsph.pdb
run
EOF
  cpptraj -i cpptraj.in &> /dev/null
  rm -f cpptraj.in cpptraj.log

  cd ..
}

function call_xcg {
  # ENM限制XCG-martini的Calpha碳，以及28.0-20.0A的UA的Calpha碳
  enm_cutoff=$1
  enm_k=10.0 #kcal/A2
  mkdir -p prot_watsph_relax_openmm; cd prot_watsph_relax_openmm

  ln -sf ../top/prot_watsph.* ./

  cat << EOF > xcg.1.in
assembling
  template_dirs {
    $MA_PATH
    $UA_PATH
    $AA_PATH
    $LI_PATH
  }
  mixdir mix.mix
  auto_mix 1
  # 默认1号划分（XCG-Martini-harmonic）
  id_default 1
  id_except_region {
     # 例：这些:1825周围20.0A以内采用2号划分
     # "1825" 20.0 2
EOF

  # XCG-UA-LFY
  for r in $active_resname; do
    for id in `get_resid $r prot_watsph.pdb`; do
      echo "     \"$id\" 28.0 2" >> xcg.1.in
    done
  done

  # AA
  for r in $active_resname; do
    for id in `get_resid $r prot_watsph.pdb`; do
      echo "     \"$id\" 20.0 3" >> xcg.1.in
    done
  done

  cat << EOF >> xcg.1.in
  }
  # 这些resid只采用4号划分（AA）
  id_except_molname {
    ATA 4 
    ATB 4
    ATC 4
    ATD 4
    ATE 4
    ATF 4
    ATG 4
    KAA 4
    KBB 4
    KCC 4
    KDD 4
    KEE 4
    KFF 4
    KGG 4
    MGA 4
    MGB 4
    MGC 4
    MGD 4
    MGE 4
    MGF 4
    MGG 4
    Na+ 4
    WAT 4
  }
/
EOF

  xcg assembling -i xcg.1.in -p prot_watsph.pdb -g 1

  cat << EOF > xcg.2.in
amber
  prmtop prot_watsph.prmtop
  inpcrds {
    prot_watsph.rst7
  }
/

restart
  # 额外增加弹性网络模型
  read_indexes 0
  file mix.gz
  mixdir mix.mix
/

sbbonded
/

enm
  style 1
  k $enm_k
  cutoff $enm_cutoff
  mask {
EOF

grep "\-XA[1-9]\-Xt[1-9]" mix.data | awk '{print $1}' >> xcg.2.in
#grep "\-CA-CX" mix.data | grep "# U" | awk '{print $1}' >> xcg.2.in

cat << EOF >> xcg.2.in
  }
/
EOF

  xcg mapping -i xcg.2.in -o multiscale

  rm -rf mix.* topol_*
  cd ..
}

# main
for f in atp atp_mg_k mg_k; do
  cd $f

  strip_solvents_outside 28.0
  call_xcg 10.0 10.0

  cd ..
done
