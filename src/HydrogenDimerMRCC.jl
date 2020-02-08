module HydrogenDimerMRCC

using Printf
import Base.@kwdef

@kwdef struct SystemParameters
	R::Float64
	r1::Float64
	theta1::Float64
	phi1::Float64
	r2::Float64
	theta2::Float64
	phi2::Float64
end

@kwdef struct MRCCInputFileInfo
	titlecomment::String
	calctype::String
	memoryMB::Int64
	threshold::Int64
	basis::String
	geometry::String
end

struct PolarInfo
	bondlength::Float64  # bond length of the hydrogen molecule
	theta::Float64       # space-fixed polar angle of the hydrogen molecule
	phi::Float64         # space-fixed azimuthal angle of the hydrogen molecule
	
	function PolarInfo(bondlength::Float64, theta::Float64, phi::Float64)
		if bondlength <= 0
			error("The molecule's bond length must be positive. ENTERED: $bondlength")
		end
		
		if !(0.0 <= theta <= pi)
			error("The polar angle must be between 0 and pi. ENTERED: $theta")
		end
		
		if !(0.0 <= phi <= 2*pi)
			error("The azimuthal angle must be between 0 and 2*pi. ENTERED: $phi")
		end
		
		new(bondlength, theta, phi)
	end
end

struct Cartesian3D
	x::Float64
	y::Float64
	z::Float64
end

function get_coordinate_line(pos::Cartesian3D)
	# returns a string for the cartesian coordinates of a hydrogen atom
	# for an MRCC input file
	
	xpos_str = "$(@sprintf("% .9f", pos.x))"
	ypos_str = "$(@sprintf("% .9f", pos.y))"
	zpos_str = "$(@sprintf("% .9f", pos.z))"
	
	hatom_str = "H    $xpos_str   $ypos_str   $zpos_str\n"
	return hatom_str
end

function get_atom_positions(hmol::PolarInfo, com::Cartesian3D)
	# find the positions of the atoms relative to the centre of mass (com)
	shift = 0.5 * hmol.bondlength
	
	dx = shift * sin(hmol.theta) * cos(hmol.phi)
	dy = shift * sin(hmol.theta) * sin(hmol.phi)
	dz = shift * cos(hmol.theta)
	
	# return the positions of the atoms in this hydrogen molecule
	atom_pos1 = Cartesian3D(com.x + dx, com.y + dy, com.z + dz)
	atom_pos2 = Cartesian3D(com.x - dx, com.y - dy, com.z - dz)
	
	return atom_pos1, atom_pos2
end

function write_input_file(filename::AbstractString, mrcc::MRCCInputFileInfo)
	# Create the file
	fout = open(filename, "w")
	
	# Create the contents of the file
	write(fout, "# $(mrcc.titlecomment)\n")
	write(fout, "calc=$(mrcc.calctype)\n")
	write(fout, "mem=$(mrcc.memoryMB)MB\n")
	write(fout, "cctol=$(mrcc.threshold)\n")
	write(fout, "basis=$(mrcc.basis)\n")
	write(fout, "\n")
	write(fout, "unit=angs\n")
	write(fout, "geom=xyz\n")
	write(fout, "4\n")
	write(fout, "$(mrcc.geometry)\n")
	
	# Close it
	close(fout)
end

function make_mrcc_input_file(sparam::SystemParameters, pathtofilename::AbstractString)
	# Create the polar coordinates and centre of masses of the two molecules
	# The first is at the origin, the second is at a distance R on the z-axis
	mol1 = PolarInfo(sparam.r1, sparam.theta1, sparam.phi1)
	com1 = Cartesian3D(0.0, 0.0, 0.0)
	mol2 = PolarInfo(sparam.r2, sparam.theta2, sparam.phi2)
	com2 = Cartesian3D(0.0, 0.0, sparam.R)
	
	# Get the 3D cartesian coordinates for the four atoms
	atomA, atomB = get_atom_positions(mol1, com1)
	atomC, atomD = get_atom_positions(mol2, com2)
	
	# Get the strings for the cartesian coordinates for the MRCC file
	# They will be used to create the 'geometry block' of the MRCC file
	lines = [get_coordinate_line(atom) for atom in [atomA, atomB, atomC, atomD]]
	geometryblock = join(lines)
	
	# Collect the information needed for the MRCC Input File
	# Most of this information is fixed for this set of calculations
	mrcc = MRCCInputFileInfo(titlecomment = "hydrogen dimer",
	                         calctype = "ccsd(t)",
	                         memoryMB = 8192,
	                         threshold = 10,
	                         basis = "aug-cc-pVQZ",
	                         geometry = geometryblock)
	
	# Create the MRCC input file
	write_input_file(pathtofilename, mrcc)
end

sparam = SystemParameters(6.0, 0.776777, pi/4.0, pi/4.0, 0.776777, pi/4.0, pi/2.0)
filename = "MRCCinputfiles/MINP"
make_mrcc_input_file(sparam, filename)


end # module
