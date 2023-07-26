#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.9.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/anders/dev/NetFlux/testing/meshing')

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
Box_1 = geompy.MakeBoxDXDYDZ(100, 10, 10)
Partition_1 = geompy.MakePartition([Box_1], [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
walls = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(walls, [3, 13, 23, 27, 31, 33])
[walls] = geompy.GetExistingSubObjects(Partition_1, False)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Box_1, 'Box_1' )
geompy.addToStudy( Partition_1, 'Partition_1' )
geompy.addToStudyInFather( Partition_1, walls, 'walls' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Partition_1,'Mesh_1')
Regular_1D = Mesh_1.Segment()
Local_Length_1 = Regular_1D.LocalLength(1,None,1e-07)
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Hexa_3D = Mesh_1.Hexahedron(algo=smeshBuilder.Hexa)
walls_1 = Mesh_1.GroupOnGeom(walls,'walls',SMESH.FACE)
isDone = Mesh_1.Compute()
[ walls_1 ] = Mesh_1.GetGroups()
isDone = Mesh_1.Compute()
[ walls_1 ] = Mesh_1.GetGroups()
status = Mesh_1.RemoveHypothesis(Local_Length_1)
Local_Length_2 = Regular_1D.LocalLength(0.2,None,1e-07)
isDone = Mesh_1.Compute()
[ walls_1 ] = Mesh_1.GetGroups()
isDone = Mesh_1.Compute()
[ walls_1 ] = Mesh_1.GetGroups()
Local_Length_2.SetLength( 0.5 )
Local_Length_2.SetPrecision( 1e-07 )
isDone = Mesh_1.Compute()
[ walls_1 ] = Mesh_1.GetGroups()
Local_Length_2.SetLength( 0.3 )
Local_Length_2.SetPrecision( 1e-07 )
isDone = Mesh_1.Compute()
[ walls_1 ] = Mesh_1.GetGroups()


## Set names of Mesh objects
smesh.SetName(walls_1, 'walls')
smesh.SetName(Hexa_3D.GetAlgorithm(), 'Hexa_3D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Local_Length_1, 'Local Length_1')
smesh.SetName(Local_Length_2, 'Local Length_2')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
