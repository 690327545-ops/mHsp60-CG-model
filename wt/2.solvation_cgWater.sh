#!/bin/bash

CGWAT_PATH=$XCG_PATH/templates/solvents/p1
SOD_PATH=`pwd`/top_ligands
#active_resname="ATA ATB ATC ATD ATE ATF ATG KAA KBB KCC KDD KEE KFF KGG MGA MGB MGC MGD MGE MGF MGG"
active_resname="ATC KAC MGC"

function packmol_box {

  N_sod_watSph=`cat top/prot_watsph.pdb | grep NA | wc -l`
  N_sod_orig=`cat top/orig.pdb | grep NA | wc -l`
  N_wat_watSph=`cat top/prot_watsph.pdb | grep "O   WAT" | wc -l`
  N_wat_orig=`cat top/orig.pdb | grep "O   WAT" | wc -l`

  N_sod_add=`echo "${N_sod_orig}-${N_sod_watSph}" | bc`
  N_wat_add=`echo "${N_wat_orig}-${N_wat_watSph}" | bc`
  N_wat_add=`echo "$N_wat_add/2" | bc`
  
  cryst_x_half=`head -n 1 top/orig.pdb | awk '{print $2}'`
  cryst_y_half=`head -n 1 top/orig.pdb | awk '{print $3}'`
  cryst_z_half=`head -n 1 top/orig.pdb | awk '{print $4}'`

  cryst_x_half=`echo "scale=3;$cryst_x_half/2+10.0" | bc`
  cryst_y_half=`echo "scale=3;$cryst_y_half/2+10.0" | bc`
  cryst_z_half=`echo "scale=3;$cryst_z_half/2+10.0" | bc`

  cat << EOF > packmol.in
tolerance 3.0
filetype pdb
output mixture.pdb

structure prot_watsph_relax_openmm/multiscale.pdb
  number 1
  fixed 0. 0. 0. 0. 0. 0.
  center
end structure
EOF
  packmol < packmol.in &> /dev/null

  exc_idx=()
  for res in $active_resname; do
    grep $res top/orig.pdb | grep ATOM | awk '{print substr($0,31,24)}' > tmp
    while read mask; do
      idx=`grep "$mask" prot_watsph_relax_openmm/multiscale.pdb | awk '{print $2}'`
      exc_idx+=($idx)
    done < tmp
    rm -f tmp
  done

  exc_mask=""
  for i in ${exc_idx[*]}; do
    mask=`printf "ATOM  %-6d" $i`
    exc_mask="$exc_mask  outside sphere `grep "$mask" mixture.pdb | awk '{print substr($0,31,24)}'` 28.0\n"
  done

  cat << EOF > packmol.in
tolerance 3.0
filetype pdb
output mixture.pdb

structure $CGWAT_PATH/p1p.pdb
  number $N_wat_add
  inside box -$cryst_x_half -$cryst_y_half -$cryst_z_half $cryst_x_half $cryst_y_half $cryst_z_half
`echo -e "$exc_mask"`
end structure

structure $CGWAT_PATH/p1n.pdb
  number $N_wat_add
  inside box -$cryst_x_half -$cryst_y_half -$cryst_z_half $cryst_x_half $cryst_y_half $cryst_z_half
`echo -e "$exc_mask"`
end structure

structure $SOD_PATH/Na+.pdb
  number $N_sod_add
  inside box -$cryst_x_half -$cryst_y_half -$cryst_z_half $cryst_x_half $cryst_y_half $cryst_z_half
`echo -e "$exc_mask"`
end structure

structure prot_watsph_relax_openmm/multiscale.pdb
  number 1
  fixed 0. 0. 0. 0. 0. 0.
  center
end structure

EOF
  packmol < packmol.in
}

function call_xcg {
  cat << EOF > xcg.in
assembling
  mollist {
    CM1 prot_watsph_relax_openmm/multiscale.gz
    Na+ $SOD_PATH/Na+.gz
    P1N $CGWAT_PATH/p1n.gz
    P1P $CGWAT_PATH/p1p.gz
  }
  auto_mix 1
  mixdir mix.mix
/
EOF
  xcg assembling -i xcg.in -p mixture.pdb
}

# main
for f in atp atp_mg_k mg_k; do
  cd $f
  packmol_box
  call_xcg
  cd ..
done
