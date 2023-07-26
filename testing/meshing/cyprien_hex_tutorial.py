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
Box_1 = geompy.MakeBoxDXDYDZ(200, 200, 200)
[Edge_1,Edge_2,Edge_3,Edge_4,Edge_5,Edge_6,Edge_7,Edge_8,Edge_9,Edge_10,Edge_11,Edge_12] = geompy.ExtractShapes(Box_1, geompy.ShapeType["EDGE"], True)
Cylinder_1 = geompy.MakeCylinderRH(50, 200)
geompy.TranslateDXDYDZ(Cylinder_1, 100, 100, 0)
[Edge_13,Edge_14,Edge_15] = geompy.ExtractShapes(Cylinder_1, geompy.ShapeType["EDGE"], True)
Box_2 = geompy.MakeBoxDXDYDZ(40, 40, 200)
geompy.TranslateDXDYDZ(Box_2, 80, 80, 0)
[Edge_13,Edge_14,Edge_15,Edge_16,Edge_17,Edge_18,Edge_19,Edge_20,Edge_21,Edge_22,Edge_23,Edge_24] = geompy.ExtractShapes(Box_2, geompy.ShapeType["EDGE"], True)
Edge_17_vertex_2 = geompy.GetSubShape(Edge_17, [2])
Edge_5_vertex_2 = geompy.GetSubShape(Edge_5, [2])
Extrusion_1 = geompy.MakePrism(Edge_21, Edge_17_vertex_2, Edge_5_vertex_2)
Edge_14_vertex_3 = geompy.GetSubShape(Edge_14, [3])
Edge_1_vertex_3 = geompy.GetSubShape(Edge_1, [3])
Extrusion_2 = geompy.MakePrism(Edge_13, Edge_14_vertex_3, Edge_1_vertex_3)
Edge_19_vertex_2 = geompy.GetSubShape(Edge_19, [2])
Edge_10_vertex_2 = geompy.GetSubShape(Edge_10, [2])
Extrusion_3 = geompy.MakePrism(Edge_24, Edge_19_vertex_2, Edge_10_vertex_2)
Edge_14_vertex_2 = geompy.GetSubShape(Edge_14, [2])
Edge_2_vertex_2 = geompy.GetSubShape(Edge_2, [2])
Extrusion_4 = geompy.MakePrism(Edge_16, Edge_14_vertex_2, Edge_2_vertex_2)
Cylinder_2 = geompy.MakeCylinderRH(50, 200)
geompy.Rotate(Cylinder_2, OZ, 45*math.pi/180.0)
geompy.TranslateDXDYDZ(Cylinder_2, 100, 100, 0)
Partition_1 = geompy.MakePartition([Box_1, Box_2, Extrusion_1, Extrusion_2, Extrusion_3, Extrusion_4, Cylinder_2], [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
walls = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(walls, [168, 86, 124, 21, 45, 55, 79, 141, 170, 129, 146, 161, 93, 112, 38, 69, 90, 62, 158, 4, 107])
geompy.DifferenceIDs(walls, [168, 86, 124, 21, 45, 55, 79, 141, 170, 129, 146, 161, 93, 112, 38, 69, 90, 62, 158, 4, 107])
geompy.UnionIDs(walls, [168, 86, 124, 21, 45, 55, 79, 141, 170, 129, 146, 161, 93, 112, 38, 69, 90, 62, 158, 4, 107, 31])
[walls] = geompy.GetExistingSubObjects(Partition_1, False)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Box_1, 'Box_1' )
geompy.addToStudy( Cylinder_1, 'Cylinder_1' )
geompy.addToStudyInFather( Box_1, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Box_1, Edge_2, 'Edge_2' )
geompy.addToStudy( Box_2, 'Box_2' )
geompy.addToStudyInFather( Box_1, Edge_3, 'Edge_3' )
geompy.addToStudyInFather( Box_1, Edge_4, 'Edge_4' )
geompy.addToStudyInFather( Box_1, Edge_5, 'Edge_5' )
geompy.addToStudyInFather( Box_1, Edge_6, 'Edge_6' )
geompy.addToStudyInFather( Box_1, Edge_7, 'Edge_7' )
geompy.addToStudyInFather( Box_1, Edge_8, 'Edge_8' )
geompy.addToStudyInFather( Box_1, Edge_9, 'Edge_9' )
geompy.addToStudyInFather( Box_1, Edge_10, 'Edge_10' )
geompy.addToStudyInFather( Box_1, Edge_11, 'Edge_11' )
geompy.addToStudyInFather( Box_1, Edge_12, 'Edge_12' )
geompy.addToStudyInFather( Box_2, Edge_13, 'Edge_13' )
geompy.addToStudyInFather( Box_2, Edge_14, 'Edge_14' )
geompy.addToStudyInFather( Box_2, Edge_15, 'Edge_15' )
geompy.addToStudyInFather( Box_2, Edge_16, 'Edge_16' )
geompy.addToStudyInFather( Box_2, Edge_17, 'Edge_17' )
geompy.addToStudyInFather( Box_2, Edge_18, 'Edge_18' )
geompy.addToStudyInFather( Box_2, Edge_19, 'Edge_19' )
geompy.addToStudyInFather( Box_2, Edge_20, 'Edge_20' )
geompy.addToStudyInFather( Box_2, Edge_21, 'Edge_21' )
geompy.addToStudyInFather( Box_2, Edge_22, 'Edge_22' )
geompy.addToStudyInFather( Box_2, Edge_23, 'Edge_23' )
geompy.addToStudyInFather( Box_2, Edge_24, 'Edge_24' )
geompy.addToStudyInFather( Edge_17, Edge_17_vertex_2, 'Edge_17:vertex_2' )
geompy.addToStudyInFather( Edge_5, Edge_5_vertex_2, 'Edge_5:vertex_2' )
geompy.addToStudy( Extrusion_1, 'Extrusion_1' )
geompy.addToStudyInFather( Edge_14, Edge_14_vertex_3, 'Edge_14:vertex_3' )
geompy.addToStudyInFather( Edge_1, Edge_1_vertex_3, 'Edge_1:vertex_3' )
geompy.addToStudy( Extrusion_2, 'Extrusion_2' )
geompy.addToStudyInFather( Edge_19, Edge_19_vertex_2, 'Edge_19:vertex_2' )
geompy.addToStudyInFather( Edge_10, Edge_10_vertex_2, 'Edge_10:vertex_2' )
geompy.addToStudy( Extrusion_3, 'Extrusion_3' )
geompy.addToStudyInFather( Edge_14, Edge_14_vertex_2, 'Edge_14:vertex_2' )
geompy.addToStudyInFather( Edge_2, Edge_2_vertex_2, 'Edge_2:vertex_2' )
geompy.addToStudy( Extrusion_4, 'Extrusion_4' )
geompy.addToStudy( Cylinder_2, 'Cylinder_2' )
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

Local_Length_1 = smesh.CreateHypothesis('LocalLength')
Local_Length_1.SetLength( 5 )
Local_Length_1.SetPrecision( 1e-07 )
Regular_1D = smesh.CreateHypothesis('Regular_1D')
Quadrangle_2D = smesh.CreateHypothesis('Quadrangle_2D')
Hexa_3D = smesh.CreateHypothesis('Hexa_3D')
Local_Length_2 = smesh.CreateHypothesis('LocalLength')
Local_Length_2.SetLength( 5 )
Local_Length_2.SetPrecision( 1e-07 )
Local_Length_3 = smesh.CreateHypothesis('LocalLength')
Local_Length_3.SetLength( 5 )
Local_Length_3.SetPrecision( 1e-07 )
Mesh_1 = smesh.Mesh(Partition_1,'Mesh_1')
status = Mesh_1.AddHypothesis(Local_Length_3)
status = Mesh_1.AddHypothesis(Regular_1D)
status = Mesh_1.AddHypothesis(Quadrangle_2D)
status = Mesh_1.AddHypothesis(Hexa_3D)
walls_1 = Mesh_1.GroupOnGeom(walls,'walls',SMESH.FACE)
isDone = Mesh_1.Compute()
[ walls_1 ] = Mesh_1.GetGroups()


## Set names of Mesh objects
smesh.SetName(Hexa_3D, 'Hexa_3D')
smesh.SetName(Quadrangle_2D, 'Quadrangle_2D')
smesh.SetName(Regular_1D, 'Regular_1D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(walls_1, 'walls')
smesh.SetName(Local_Length_1, 'Local Length_1')
smesh.SetName(Local_Length_2, 'Local Length_2')
smesh.SetName(Local_Length_3, 'Local Length_3')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
