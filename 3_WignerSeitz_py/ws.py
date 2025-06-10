from ovito.io import *
from ovito.data import *
from ovito.modifiers import *
from ovito.pipeline import *
import numpy as np

import ovito
ovito.enable_logging()

# Import Data File
pipeline = import_file("./datafile_input/DisplacedLattice")

# Perform Wigner-Seitz analysis:
ws = WignerSeitzAnalysisModifier(
    per_type_occupancies = False,
    affine_mapping = ReferenceConfigurationModifier.AffineMapping.ToReference)
ws.reference = FileSource()
ws.reference.load("./datafile_input/ReferenceLattice")
pipeline.modifiers.append(ws)

# Expression Selection 1
ep_1 = ExpressionSelectionModifier(
    expression = 'Occupancy != 1')
pipeline.modifiers.append(ep_1)
data = pipeline.compute()
print("Atoms in ep_1:",data.attributes['ExpressionSelection.count'])

# Invert Selection
pipeline.modifiers.append(InvertSelectionModifier())

# Delete Selected
pipeline.modifiers.append(DeleteSelectedModifier())

# Export the Lammps/dump coordinates of just the antisites by removing all other atoms.
export_file(pipeline, "./datafile_output/ws.dump", "lammps/dump",
    columns = ["Particle Identifier", "Particle Type",'Position.X', 'Position.Y', 'Position.Z',"Occupancy"],
    multiple_frames = True)
