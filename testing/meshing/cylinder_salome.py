#!/usr/bin/python3

###
### This file is generated automatically by SALOME v9.9.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/anders/dev/Compress3D/meshing')

###
### SHAPER component
###

from salome.shaper import model

model.begin()
partSet = model.moduleDocument()
model.end()

###
### SHAPERSTUDY component
###

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Cylinder_1 = geompy.MakeCylinderRH(100, 300)
side = geompy.CreateGroup(Cylinder_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(side, [3])
top = geompy.CreateGroup(Cylinder_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(top, [10])
bottom = geompy.CreateGroup(Cylinder_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(bottom, [12])
[side, top, bottom] = geompy.GetExistingSubObjects(Cylinder_1, False)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Cylinder_1, 'Cylinder_1' )
geompy.addToStudyInFather( Cylinder_1, side, 'side' )
geompy.addToStudyInFather( Cylinder_1, top, 'top' )
geompy.addToStudyInFather( Cylinder_1, bottom, 'bottom' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Cylinder_1,'Mesh_1')
NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
side_1 = Mesh_1.GroupOnGeom(side,'side',SMESH.FACE)
top_1 = Mesh_1.GroupOnGeom(top,'top',SMESH.FACE)
bottom_1 = Mesh_1.GroupOnGeom(bottom,'bottom',SMESH.FACE)
isDone = Mesh_1.Compute()
[ side_1, top_1, bottom_1 ] = Mesh_1.GetGroups()


## Set names of Mesh objects
smesh.SetName(top_1, 'top')
smesh.SetName(bottom_1, 'bottom')
smesh.SetName(side_1, 'side')
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()

from cfdmsh import *
ExportSU2File(mesh=Mesh_1,file="cylinder")