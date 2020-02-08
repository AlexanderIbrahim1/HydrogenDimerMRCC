# HydrogenDimerMRCC

The package 'HydrogenDimerMRCC' provides a method of constructing input files for the hydrogen dimer system, used with the MRCC electronic structure software [M. Kallay *et al.*, *J. Chem. Phys.* **139**, 094105 (2013)]. It creates the hydrogen dimer system using seven parameters:

*R* - the centre-of-mass intermolecular distance between the two hydrogen molecules
*r~1~*, *r~2~* - the bond lengths of the two molecules
*&theta;~1~*, *&theta;~2~* - the space-fixed polar angles of the two molecules
*&phi;~1~*, *&phi;~2~* - the space-fixed azimuthal angles of the two molecules

In the context of this project, a large number (several thousands) of electronic structure energy calculations for the hydrogen dimer system must be performed on a large cluster. The technical details (basis set, memory usage, etc.) of the input files are, for the most part, identical. These details are hard-coded based on calculations done in a paper by R. J. Hinde [R. J. Hinde, *J. Chem. Phys.* **128**, 154308 (2008)]. The input files differ only by the positions of the hydrogen atoms.

## Examples

### Example 1

Create a single MRCC input file in the directory "MRCCinputfiles/example1".
```Julia
using HydrogenDimerMRCC

r0 = 0.776777 # ground state vibrational bond length
R  = 6.0      # centre-of-mass intermolecular distance
dimer1 = PolarCoords(r0, pi/2, pi/4)
dimer2 = PolarCoords(r0, pi/4, 4*pi/3)
dparam = DimerParameters(R, dimer1, dimer2)

filename = "MRCCinputfiles/MINP"
make_mrcc_input_file(dparam, filename)
```

The resulting file will look like:

```txt
# hydrogen dimer
calc=ccsd(t)
mem=8192MB
cctol=10
basis=aug-cc-pVQZ

unit=angs
geom=xyz
4
H     0.274632142    0.274632142    0.000000000
H    -0.274632142   -0.274632142   -0.000000000
H    -0.137316071   -0.237838412    6.274632142
H     0.137316071    0.237838412    5.725367858
```

### Example 2

Create a set of MRCC input files, where one molecule lies along the z-axis, while the other molecule goes through several polar angles.
```Julia
using HydrogenDimerMRCC

r0 = 0.776777
R  = 6.0
theta_values = range(0.0, stop = pi, length = 7)

dimer1 = PolarCoords(r0, 0.0, 0.0)
for (i, theta) in enumerate(theta_values)
	dimer2 = PolarCoords(r0, theta, 0.0)
	dparam = DimerParameters(R, dimer1, dimer2)
	make_mrcc_input_file(dparam, "MRCCinputfiles/example2/MINP$i")
end
```

## Future Work

- allow finer control of the calculation details (basis set, etc.)
- implement ghost atoms (atoms whose charges have been removed from the calculation, but whose basis set remains)
