using HydrogenDimerMRCC
using Test

@testset "HydrogenDimerMRCC" begin
	@testset "Invalid Inputs" begin
		# negative bond length
		@test_throws ErrorException p = PolarCoords(-1.0, pi/2, pi/2)
		# zero bond length
		@test_throws ErrorException p = PolarCoords(0.0, pi/2, pi/4)
		# polar angle below bounds
		@test_throws ErrorException p = PolarCoords(1.0, -0.5*pi, pi/2)
		# polar angle above bounds
		@test_throws ErrorException p = PolarCoords(1.0, 1.5*pi, pi/2)
		# azimuthal angle below bounds
		@test_throws ErrorException p = PolarCoords(1.0, pi/2, -0.5*pi)
		# azimuthal angle above bounds
		@test_throws ErrorException p = PolarCoords(1.0, pi/2, 2.5*pi)
		# invalid distance directly into DimerParameters
		@test_throws ErrorException d = DimerParameters(-5.0, 1.0, pi/2, pi/2, 1.0, pi/2, pi/2)
		@test_throws ErrorException d = DimerParameters(0.0, 1.0, pi/2, pi/2, 1.0, pi/2, pi/2)
	end
	
	@testset "Ideal Case 1" begin
		import HydrogenDimerMRCC: atom_positions, Cartesian3D, distance
		
		# rotate dimer (at origin) of length 2 onto x-y plane
		# one atom should be at  ( 1, 0, 0)
		# the other should be at (-1, 0, 0)
		polar = PolarCoords(2.0, pi/2, 0.0)
		com   = Cartesian3D(0.0, 0.0, 0.0)
		atom1, atom2 = atom_positions(polar, com)
		
		ideal1 = Cartesian3D( 1.0, 0.0, 0.0)
		ideal2 = Cartesian3D(-1.0, 0.0, 0.0)
		
		tolerance = 1.0e-15
		@test distance(atom1, ideal1) < tolerance
		@test distance(atom2, ideal2) < tolerance
	end
	
	@testset "Ideal Case 2" begin
		import HydrogenDimerMRCC: atom_positions, Cartesian3D, distance
		
		# rotate dimer (at origin) of length 2 onto x-y plane
		# then rotate 45 degrees in x-y place
		# one atom should be at  ( 1/sqrt(2),  1/sqrt(2), 0)
		# the other should be at (-1/sqrt(2), -1/sqrt(2), 0)
		polar = PolarCoords(2.0, pi/2, pi/4)
		com   = Cartesian3D(0.0, 0.0, 0.0)
		atom1, atom2 = atom_positions(polar, com)
		
		ideal1 = Cartesian3D( 1.0/sqrt(2.0),  1.0/sqrt(2.0), 0.0)
		ideal2 = Cartesian3D(-1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0)
		
		tolerance = 1.0e-15
		@test distance(atom1, ideal1) < tolerance
		@test distance(atom2, ideal2) < tolerance
	end
end
