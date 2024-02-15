import gmsh
gmsh.initialize()

include("FEM.jl")
using .FEM

gmsh.model.add("post")
gmsh.merge("/home/perebal/Dokumentumok/GitHub/Oktatas/src/postproc-test.geo")
gmsh.merge("/home/perebal/Dokumentumok/GitHub/Oktatas/src/postproc-test.msh")

elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(2, -1)
elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
S = gmsh.view.add("field")
s = []
push!(s, [1.,1.,1.])
push!(s, [2.,2.,2.,2.])

e = [13,12]

gmsh.view.addModelData(S, 0, "post", "ElementNodeData", e, s, 0., 1)

display(e)
display(s)


gmsh.fltk.run()
gmsh.finalize()
