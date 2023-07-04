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
box = geompy.MakeBoxDXDYDZ(100, 10, 10)
sidewalls = geompy.CreateGroup(box, geompy.ShapeType["FACE"])
geompy.UnionIDs(sidewalls, [33, 23, 27, 31])
outlet = geompy.CreateGroup(box, geompy.ShapeType["FACE"])
geompy.UnionIDs(outlet, [13])
inlet = geompy.CreateGroup(box, geompy.ShapeType["FACE"])
geompy.UnionIDs(inlet, [3])
[sidewalls, outlet, inlet] = geompy.GetExistingSubObjects(box, False)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( box, 'box' )
geompy.addToStudyInFather( box, sidewalls, 'sidewalls' )
geompy.addToStudyInFather( box, outlet, 'outlet' )
geompy.addToStudyInFather( box, inlet, 'inlet' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(box,'Mesh_1')
NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
sidewalls_1 = Mesh_1.GroupOnGeom(sidewalls,'sidewalls',SMESH.FACE)
outlet_1 = Mesh_1.GroupOnGeom(outlet,'outlet',SMESH.FACE)
inlet_1 = Mesh_1.GroupOnGeom(inlet,'inlet',SMESH.FACE)
isDone = Mesh_1.Compute()
[ sidewalls_1, outlet_1, inlet_1 ] = Mesh_1.GetGroups()
NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()
NETGEN_3D_Parameters_1.SetMaxSize( 5 )
NETGEN_3D_Parameters_1.SetMinSize( 3.3665 )
NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 2 )
NETGEN_3D_Parameters_1.SetChordalError( -1 )
NETGEN_3D_Parameters_1.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
NETGEN_3D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_1.SetCheckChartBoundary( 3 )
isDone = Mesh_1.Compute()
[ sidewalls_1, outlet_1, inlet_1 ] = Mesh_1.GetGroups()
status = Mesh_1.RemoveHypothesis(NETGEN_3D_Parameters_1)
NETGEN_3D_Parameters_2 = NETGEN_1D_2D_3D.Parameters()
NETGEN_3D_Parameters_2.SetMaxSize( 10.0995 )
NETGEN_3D_Parameters_2.SetMinSize( 1 )
NETGEN_3D_Parameters_2.SetSecondOrder( 0 )
NETGEN_3D_Parameters_2.SetOptimize( 1 )
NETGEN_3D_Parameters_2.SetFineness( 2 )
NETGEN_3D_Parameters_2.SetChordalError( -1 )
NETGEN_3D_Parameters_2.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_2.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_2.SetFuseEdges( 1 )
NETGEN_3D_Parameters_2.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_2.SetCheckChartBoundary( 3 )
isDone = Mesh_1.Compute()
[ sidewalls_1, outlet_1, inlet_1 ] = Mesh_1.GetGroups()
status = Mesh_1.RemoveHypothesis(NETGEN_3D_Parameters_2)
NETGEN_3D_Parameters_3 = NETGEN_1D_2D_3D.Parameters()
NETGEN_3D_Parameters_3.SetMaxSize( 2 )
NETGEN_3D_Parameters_3.SetMinSize( 0.1 )
NETGEN_3D_Parameters_3.SetSecondOrder( 0 )
NETGEN_3D_Parameters_3.SetOptimize( 1 )
NETGEN_3D_Parameters_3.SetFineness( 2 )
NETGEN_3D_Parameters_3.SetChordalError( -1 )
NETGEN_3D_Parameters_3.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_3.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_3.SetFuseEdges( 1 )
NETGEN_3D_Parameters_3.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_3.SetCheckChartBoundary( 3 )
isDone = Mesh_1.Compute()
[ sidewalls_1, outlet_1, inlet_1 ] = Mesh_1.GetGroups()


## Set names of Mesh objects
smesh.SetName(outlet_1, 'outlet')
smesh.SetName(inlet_1, 'inlet')
smesh.SetName(sidewalls_1, 'sidewalls')
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
smesh.SetName(NETGEN_3D_Parameters_2, 'NETGEN 3D Parameters_2')
smesh.SetName(NETGEN_3D_Parameters_3, 'NETGEN 3D Parameters_3')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
