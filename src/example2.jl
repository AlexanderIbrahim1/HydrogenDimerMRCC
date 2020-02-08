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
