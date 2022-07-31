# OFW17---Training-Profile-Extrusion-Die-Design
Training session from OFW17 - Using OpenFOAM to Design Extrusion Dies for Thermoplastic Profiles

# Tools used
* ### FreeCAD v0.20
* ### Windows Subsystem Linux - Ubuntu 20.04
* ### OpenFOAM v21.12


# Main steps:
1. in the Module "Part Desing" edit the geometry file in FreeCAD and use the export command to create an stl file or each patch (use ast format to export the stl file in ascii, then change the extension to stl), using the shape boundary binder
2. use the uniquSTL.X script (stl folder) to create the total.fms file (the myList file should contain the list of the individual stl files to join)
3. copy the total.fms file to the root of the case folder
4. compile the fobj funtion object (fobj folder) with "wmake libso"
5. Check functions subdictionary in controlDict file provided, and define there the number of elemental sections of your case
6. Run the case with simpleFoam
7. The fobj file created in the case root folder will contain the evolution of the Objective Function values
