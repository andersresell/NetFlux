from meshpy.tet import MeshInfo, build

import meshpy.triangle as mt







import os
print(os.getcwd())

mesh_info = MeshInfo()
mesh_info.set_points([
    (0,0,0), (12,0,0), (12,12,0), (0,12,0),
    (0,0,12), (12,0,12), (12,12,12), (0,12,12),
    ])
mesh_info.set_facets([
    [0,1,2,3],
    [4,5,6,7],
    [0,4,5,1],
    [1,5,6,2],
    [2,6,7,3],
    [3,7,4,0],
    ],markers=[1,2,3,4,5,6])


mesh = build(mesh_info,)
mesh = build(mesh_info=mesh_info)

print("Mesh Points:")
for i, p in enumerate(mesh.points):
    print(i, p)
print("Point numbers in tetrahedra:")
for i, t in enumerate(mesh.elements):
    print(i, t)
mesh.write_vtk("test.vtk")
