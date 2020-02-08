# HydrogenDimerMRCC

The package 'HydrogenDimerMRCC' provides a method of constructing input files for the hydrogen dimer system, used with the MRCC electronic structure software. It accepts seven parameters, ...

```Julia
using HydrogenDimerMRCC
blah blah
```

In the context of this project, a large number (several thousands) of electronic structure energy calculations for the hydrogen dimer system must be performed on a large cluster. The technical details (basis set, memory usage, etc.) of the input files are, for the most part, identical. These details are hard-coded based on the parameters given by [this paper]. The input files differ only by the positions of the hydrogen atoms.

## Examples

Create a set of MRCC input files for a single rotational orientation, but different centre-of-mass separations.
```Julia
using HydrogenDimerMRCC
blah blah
```

Changes in the hydrogen molecule bond lengths can be used to calculate vibrational excitations.
```Julia
using HydrogenDimerMRCC
blah blah
```

## Future Work

- allow finer control of the calculation details (basis set, etc.)
- implement ghost atoms (atoms whose charges have been removed from the calculation, but whose basis set remains)
