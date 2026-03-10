#!/bin/bash

for f in `ls -1 *.prmtop`; do
  cat << EOF > xcg.in
amber
  molname ${f%.*}
  prmtop $f
  inpcrds {
    ${f%.*}.rst7
  }
  tolammps ${f%.*}
/
EOF
  xcg mapping -i xcg.in
  rm -rf *.in *.data
done
