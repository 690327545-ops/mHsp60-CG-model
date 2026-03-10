#!/bin/bash

for f in ../atp/top_orig/atp.pdb ../atp_mg_k/top_orig/atp+mg+k.pdb ../mg_k/top_orig/mg+k.pdb; do
  for res in `grep ATOM $f | awk '{print substr($0,18,3)}' | uniq`; do
    if ! [ -f "$XCG_PATH/templates/united-14sb-lfy/$res.gz" ]; then
      resid=`grep "$res" ${f} | head -n 1 | awk '{print $5}'`
      cat << EOF > cpptraj.in
parm ${f%.*}.prmtop
trajin ${f%.*}.rst7
strip "!(:$resid)" parmout ${res}.prmtop
trajout ${res}.rst7
run
EOF
      cpptraj -i cpptraj.in &> /dev/null
      rm -rf cpptraj.in cpptraj.log
    fi
  done
done
