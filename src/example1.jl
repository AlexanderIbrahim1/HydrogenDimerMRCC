using HydrogenDimerMRCC

r0 = 0.776777 # ground state vibrational bond length
R  = 6.0      # centre-of-mass intermolecular distance
dimer1 = PolarCoords(r0, pi/2, pi/4)
dimer2 = PolarCoords(r0, pi/4, 4*pi/3)
dparam = DimerParameters(R, dimer1, dimer2)

filename = "MRCCinputfiles/example1/MINP"
make_mrcc_input_file(dparam, filename)
