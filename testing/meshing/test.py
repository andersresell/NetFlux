import gmsh

gmsh.initialize()

# Create a new model
model = gmsh.model.create()

# Create a new geometry
gmsh.model.geo.addPoint(0, 0, 0, 0.1, 1)
gmsh.model.geo.addPoint(1, 0, 0, 0.1, 2)
gmsh.model.geo.addPoint(1, 1, 0, 0.1, 3)
gmsh.model.geo.addPoint(0, 1, 0, 0.1, 4)
gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(2, 3, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 1, 4)
gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
gmsh.model.geo.addPlaneSurface([1], 1)

# Set the model to be the active model
gmsh.model.setCurrent(model)

# Set the meshing algorithm
gmsh.option.setNumber("Mesh.Algorithm", 6)

# Generate the mesh
gmsh.model.mesh.generate(2)

# Save the mesh to a file
gmsh.write("mesh.msh")

# Finalize Gmsh
gmsh.finalize()
