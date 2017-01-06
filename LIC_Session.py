### FAnToM Session
### API Version: 20140515
### Used core version:    GITDIR-NOTFOUND (GITDIR-NOTFOUND)
### Used toolbox version: GITDIR-NOTFOUND (GITDIR-NOTFOUND)

################################################################
###                  Reset GUI                               ###
################################################################
fantom.ui.setCamera( 0, fantom.ui.Camera( fantom.math.Vector3(17.8712, -2.70261, 52.2401), fantom.math.Vector3(17.6914, -2.57705, 51.2644), fantom.math.Vector3(-0.443396, 0.87501, 0.194314), 1, -3.46058e-310 ) )
fantom.ui.setCamera( 1, fantom.ui.Camera( fantom.math.Vector3(0, 3.8637, 0), fantom.math.Vector3(0, 2.8637, 0), fantom.math.Vector3(0, 0, 1), 0, 1 ) )
fantom.ui.setCamera( 2, fantom.ui.Camera( fantom.math.Vector3(0, 0, 3.8637), fantom.math.Vector3(0, 0, 2.8637), fantom.math.Vector3(0, 1, 0), 0, 1 ) )
fantom.ui.setCamera( 3, fantom.ui.Camera( fantom.math.Vector3(3.8637, 0, 0), fantom.math.Vector3(2.8637, 0, 0), fantom.math.Vector3(0, 0, 1), 0, 1 ) )

fantom.ui.setClippingPlane( fantom.ui.ClippingPlane( 0, fantom.math.Matrix4( (1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1) ), False ) )
fantom.ui.setClippingPlane( fantom.ui.ClippingPlane( 1, fantom.math.Matrix4( (1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1) ), False ) )
fantom.ui.setClippingPlane( fantom.ui.ClippingPlane( 2, fantom.math.Matrix4( (1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1) ), False ) )
fantom.ui.setClippingPlane( fantom.ui.ClippingPlane( 3, fantom.math.Matrix4( (1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1) ), False ) )
fantom.ui.setClippingPlane( fantom.ui.ClippingPlane( 4, fantom.math.Matrix4( (1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1) ), False ) )
fantom.ui.setClippingPlane( fantom.ui.ClippingPlane( 5, fantom.math.Matrix4( (1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1) ), False ) )

fantom.ui.setBackgroundColor( fantom.math.Color(1, 1, 1, 1) )

fantom.ui.setRotationCenter( fantom.ui.RotationCenter( fantom.math.Vector3(0, 0, 0), True, False, True ) )


################################################################
###                  Create algorithms                       ###
################################################################
VisPraktikum_LIC = fantom.makeAlgorithm("VisPraktikum/LIC")
VisPraktikum_LIC.setName( "VisPraktikum/LIC" )
fantom.ui.setAlgorithmPosition(VisPraktikum_LIC, fantom.math.Vector2(-3, 190))
VisPraktikum_LIC.setVisualOutputVisible('licTexture', True)
VisPraktikum_LIC.setVisualOutputVisible('noiseTexture', False)

Load_VTK = fantom.makeAlgorithm("Load/VTK")
Load_VTK.setName( "Load/VTK" )
Load_VTK.setOption("Input File", "/home/yves/Schreibtisch/Vis/TestData/markCellSmall.vtk")
Load_VTK.setOption("Big Endian", True)
Load_VTK.setOption("Dimension", "3D")
Load_VTK.setOption("Time List", "")
fantom.ui.setAlgorithmPosition(Load_VTK, fantom.math.Vector2(105, 50))

Grid_ShowGrid = fantom.makeAlgorithm("Grid/Show Grid")
Grid_ShowGrid.setName( "Grid/Show Grid" )
Grid_ShowGrid.setOption("Line color", fantom.math.Color(0, 0, 1, 1))
Grid_ShowGrid.setOption("Random jittering of color", True)
Grid_ShowGrid.setOption("Random seed", 0)
Grid_ShowGrid.setOption("Line width", 1.5)
fantom.ui.setAlgorithmPosition(Grid_ShowGrid, fantom.math.Vector2(204, 180))
Grid_ShowGrid.setVisualOutputVisible('grid', False)



################################################################
###                     Make Connections                     ###
################################################################
Load_VTK.connect("Fields", VisPraktikum_LIC, "Field")
Load_VTK.connect("Grid", Grid_ShowGrid, "Grid")


################################################################
###                      Run algorithms                      ###
################################################################
fantom.scheduleAllNecessary()