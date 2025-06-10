from ovito.io import import_file, export_file
from ovito.modifiers import *
from ovito.data import DislocationNetwork

import ovito
ovito.enable_logging()

# Import
pipeline = import_file("./datafile_input/input.dump")

# Export dump.0:
export_file(pipeline, "./datafile_output/dump.0", "lammps/dump",
    columns = ["Particle Identifier", "Particle Type",'Position.X', 'Position.Y', 'Position.Z'],
    multiple_frames = True)

# DXA
dxa= DislocationAnalysisModifier()
dxa.input_crystal_structure = DislocationAnalysisModifier.Lattice.BCC
pipeline.modifiers.append(dxa)
data = pipeline.compute()

total_line_length = data.attributes['DislocationAnalysis.total_line_length']
cell_volume = data.attributes['DislocationAnalysis.cell_volume']
print("Dislocation density: %f" % (total_line_length / cell_volume))

# Expression Selection 1
ep_1 = ExpressionSelectionModifier(
    expression = 'StructureType != 0')
pipeline.modifiers.append(ep_1)
data = pipeline.compute()
print("Atoms in ep_1:",data.attributes['ExpressionSelection.count'])

# Delete Selected
pipeline.modifiers.append(DeleteSelectedModifier())

# Export type0.data:
export_file(pipeline, "./datafile_output/type0.data", "lammps/data")

# Export CA file:
export_file(pipeline, "./datafile_output/disline.input", "ca")


