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
all_walls = geompy.CreateGroup(Box_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(all_walls, [27, 33, 23, 3, 13, 31])
[all_walls] = geompy.GetExistingSubObjects(Box_1, False)
[all_walls] = geompy.GetExistingSubObjects(Box_1, False)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Box_1, 'Box_1' )
geompy.addToStudyInFather( Box_1, all_walls, 'all_walls' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

NETGEN_3D_Parameters_1 = smesh.CreateHypothesis('NETGEN_Parameters', 'NETGENEngine')
NETGEN_3D_Parameters_1.SetMaxSize( 0.5 )
NETGEN_3D_Parameters_1.SetMinSize( 0.01 )
NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 4 )
NETGEN_3D_Parameters_1.SetChordalError( -1 )
NETGEN_3D_Parameters_1.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
NETGEN_3D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_1.SetCheckChartBoundary( 3 )
NETGEN_1D_2D_3D = smesh.CreateHypothesis('NETGEN_2D3D', 'NETGENEngine')
NETGEN_3D_Parameters_2 = smesh.CreateHypothesis('NETGEN_Parameters', 'NETGENEngine')
NETGEN_3D_Parameters_2.SetMaxSize( 10.0995 )
NETGEN_3D_Parameters_2.SetMinSize( 3.3665 )
NETGEN_3D_Parameters_2.SetSecondOrder( 0 )
NETGEN_3D_Parameters_2.SetOptimize( 1 )
NETGEN_3D_Parameters_2.SetFineness( 2 )
NETGEN_3D_Parameters_2.SetChordalError( -1 )
NETGEN_3D_Parameters_2.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_2.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_2.SetFuseEdges( 1 )
NETGEN_3D_Parameters_2.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_2.SetCheckChartBoundary( 3 )
NETGEN_3D_Parameters_3 = smesh.CreateHypothesis('NETGEN_Parameters', 'NETGENEngine')
NETGEN_3D_Parameters_3.SetMaxSize( 1 )
NETGEN_3D_Parameters_3.SetMinSize( 0.5 )
NETGEN_3D_Parameters_3.SetSecondOrder( 0 )
NETGEN_3D_Parameters_3.SetOptimize( 1 )
NETGEN_3D_Parameters_3.SetFineness( 2 )
NETGEN_3D_Parameters_3.SetChordalError( -1 )
NETGEN_3D_Parameters_3.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_3.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_3.SetFuseEdges( 1 )
NETGEN_3D_Parameters_3.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_3.SetCheckChartBoundary( 3 )
NETGEN_3D_Parameters_4 = smesh.CreateHypothesis('NETGEN_Parameters', 'NETGENEngine')
NETGEN_3D_Parameters_4.SetMaxSize( 1 )
NETGEN_3D_Parameters_4.SetMinSize( 0.5 )
NETGEN_3D_Parameters_4.SetSecondOrder( 0 )
NETGEN_3D_Parameters_4.SetOptimize( 1 )
NETGEN_3D_Parameters_4.SetFineness( 2 )
NETGEN_3D_Parameters_4.SetChordalError( -1 )
NETGEN_3D_Parameters_4.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_4.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_4.SetFuseEdges( 1 )
NETGEN_3D_Parameters_4.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_4.SetCheckChartBoundary( 3 )
NETGEN_3D_Parameters_5 = smesh.CreateHypothesis('NETGEN_Parameters', 'NETGENEngine')
NETGEN_3D_Parameters_5.SetMaxSize( 2 )
NETGEN_3D_Parameters_5.SetMinSize( 1 )
NETGEN_3D_Parameters_5.SetSecondOrder( 0 )
NETGEN_3D_Parameters_5.SetOptimize( 1 )
NETGEN_3D_Parameters_5.SetFineness( 4 )
NETGEN_3D_Parameters_5.SetChordalError( -1 )
NETGEN_3D_Parameters_5.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_5.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_5.SetFuseEdges( 1 )
NETGEN_3D_Parameters_5.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_5.SetCheckChartBoundary( 3 )
NETGEN_3D_Parameters_6 = smesh.CreateHypothesis('NETGEN_Parameters', 'NETGENEngine')
NETGEN_3D_Parameters_6.SetMaxSize( 1 )
NETGEN_3D_Parameters_6.SetMinSize( 0.4 )
NETGEN_3D_Parameters_6.SetSecondOrder( 0 )
NETGEN_3D_Parameters_6.SetOptimize( 1 )
NETGEN_3D_Parameters_6.SetFineness( 4 )
NETGEN_3D_Parameters_6.SetChordalError( -1 )
NETGEN_3D_Parameters_6.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_6.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_6.SetFuseEdges( 1 )
NETGEN_3D_Parameters_6.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_6.SetCheckChartBoundary( 3 )
Mesh_1 = smesh.Mesh(Box_1,'Mesh_1')
status = Mesh_1.AddHypothesis(NETGEN_3D_Parameters_6)
status = Mesh_1.AddHypothesis(NETGEN_1D_2D_3D)
all_walls_1 = Mesh_1.GroupOnGeom(all_walls,'all_walls',SMESH.FACE)
isDone = Mesh_1.Compute()
[ all_walls_1 ] = Mesh_1.GetGroups()


## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D_3D, 'NETGEN 1D-2D-3D')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(all_walls_1, 'all_walls')
smesh.SetName(NETGEN_3D_Parameters_4, 'NETGEN 3D Parameters_4')
smesh.SetName(NETGEN_3D_Parameters_5, 'NETGEN 3D Parameters_5')
smesh.SetName(NETGEN_3D_Parameters_6, 'NETGEN 3D Parameters_6')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
smesh.SetName(NETGEN_3D_Parameters_2, 'NETGEN 3D Parameters_2')
smesh.SetName(NETGEN_3D_Parameters_3, 'NETGEN 3D Parameters_3')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
