#!/bin/bash
mol="Glucose"

gmx editconf -f $mol.gro -o centered.gro -box 2.3 2.3 2.3
gmx solvate -cs tip4p -cp centered.gro -p $mol.top -o solved.gro

gmx grompp -f min.mdp -po min-run.mdp -p $mol.top -c solved.gro -pp min.top -o min.tpr -r solved.gro
gmx mdrun -v -deffnm min-run.mdp

gmx grompp -f eqi.mdp -po eqi-run.mdp -p $mol.top -c min.gro -r min.gro -pp eqi.top -o eqi.tpr
gmx mdrun -v -deffnm eqi-run.mdp

gmx grompp -f md.mdp -po md-run.mdp -p $mol.top -c eqi.gro -r eqi.gro -t eqi.trr -pp md.top -o md.tpr
gmx mdrun -v -deffnm md-run.mdp
gmx trjconv -f md.trr -s md.tpr -pbc mol -o md-nojump.trr << EOF
0
EOF

python3 T1Calculator.py $mol
