LatticeReconstructionTools
2025.6.10
===========================

This project contains a collection of auxiliary scripts and programs designed to assist the main program "LatticeReconstruct". These tools help with generating and modify input files, performing Wigner-Seitz (WS) analysis, and filtering true point defects. Each module operates independently and can be used as needed.



1. GetInputsByOvito_py: Generate Input Files for LatticeReconstruct
-------------------------------------------------------------------
This module uses Ovito's Python API to extract information from a LAMMPS dump file and generate the three required input files for LatticeReconstruct.

Directory structure:
1_GetInputsByOvito_py/
├── datafile_input/                # Input folder containing the displaced lattice (LAMMPS dump format)
│   └── input.dump
├── datafile_output/              # Output folder for the generated input files
└── GetInitialInputs.py           # Main script: run with 'python GetInitialInputs.py'

Note: Requires the Ovito Python module. Install with `pip install ovito`.
Note: A sample input is provided in the input folder.



2. GetNewCrystalAnalysisFile: Modify disline.input File
--------------------------------------------------------
When the target defect is a dislocation loop, it is necessary to modify the disline.input file in the input files to remove the target dislocation loop, thereby replacing the dislocation loop region with a perfect crystal lattice.
This module reads an initial displaced alttice configuration (using its periodic boundary only) and a disline.input file (i.e. a Crystal Analysis file), remove target dislocations, and then generating a new disline.input file and a dump file for visualization.

Directory structure:
2_GetNewCrystalAnalysisFile/
├── datafiles_inputs/             # Input folder
│   ├── disline.input             # Original disline input file
│   └── dump.0                    # Displaced lattice configuration
├── datafiles_outputs/            # Output folder
│   ├── disline_new.input         # Modified disline.input
│   └── nodes.dump                # Dump file for visualizing dislocation structures (containing all dislocation line vertex)
├── GetNewCrystalAnalysis.sh      # Run with 'bash GetNewCrystalAnalysis.sh'
├── loopremoving.f90              # Source code
├── a.out                         # Executable
├── log.dat                       # Execution log
├── word_count.sh, words.dat      # Auxiliary scripts and intermediate files

Note: In this section, the user needs to manually specify the location of the dislocation loops or dislocation lines to be removed in "GetNewCrystalAnalysis.sh" file. Note that selecting just a single point is sufficient — all dislocations within a 15 Å radius around that point, along with any dislocations connected to them, will be removed.
Note: Sample inputs/outputs are provided in the input/output folder.



3. WignerSeitz_py: Perform Wigner-Seitz (WS) Analysis
-----------------------------------------------------
This module compares the reference and displaced lattices using the WS method, identifying atoms with occupancy ≠ 1.

Directory structure:
3_WignerSeitz_py/
├── datafile_input/
│   ├── DisplacedLattice
│   └── ReferenceLattice
├── datafile_output/
│   └── ws.dump                   # WS analysis result, contains only atoms with Occupancy ≠ 1
└── ws.py                         # Main script: run with 'python ws.py'

Note: Sample inputs/outputs are provided in the input/output folder.



4. DefectFilter: Filter True Point Defects
------------------------------------------
This module filters out WS identified defect atoms near dislocation lines (often misidentified as defects when using LatticeReconstruct program) and retains only true point defects.

Directory structure:
4_DefectFilter/
├── AtomsNearDisline_f/
│   ├── SelectingAtomsAwayDislines.f90   # Fortran source code
│   ├── a.out                            # Compiled executable
│   ├── Select.sh                        # Main script: run with 'bash Select.sh'
│   ├── disline.input, ws.dump           # Input files
│   ├── output*.dat, *.dump              # Intermediate and output files
│   ├── word_count.sh, words.dat         # Auxiliary scripts and intermediate files
├── datafiles-input/
│   ├── disline.input                    # Input file 1, dislocation line positions
│   └── ws.dump                          # Input file 2, WS identified defects, in dump format
├── datafiles-output/
│   └── ws-new.dump                      # Final output: only true defects
└── GetFilteredDefects.sh                # Run with 'bash GetFilteredDefects.sh'

Note: Sample inputs/outputs are provided in the input/output folder.



Environment Requirements
------------------------
- Python 3.x
- Ovito Python Module 
- gfortran (to compile Fortran code)
- Bash Shell



Usage Recommendation
--------------------
The modules can be used in the following order:

1. Use "1_GetInputsByOvito_py" to generate the input files.
2. Use "2_GetNewCrystalAnalysisFile" to modify the disline.input file.
3. Use "3_WignerSeitz_py" to perform WS analysis.
4. Use "4_DefectFilter" to filter out false defect atoms and retain true point defects.

