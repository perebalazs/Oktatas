```julia
using LowLevelFEM
gmsh.initialize()
gmsh.open("cube.geo")
mat = material("body")
problem = Problem([mat])
γ = 0.5
x(X, Y, Z) = X + γ * Y
y(X, Y, Z) = Y
z(X, Y, Z) = Z
deformation = field("body", fx=x, fy=y, fz=z)
r = vectorField(problem, [deformation])
showDeformationResults(r, :vector, visible=true)
openPostProcessor()
F = ∇(r)
atP = "P"
FatP = probe(F, atP)
R = nodePositionVector(problem)

R = nodesToElements(R)

showDeformationResults(F * R, :vector)
r = nodesToElements(r)
showDeformationResults(inv(F) * r, :vector)
openPostProcessor()
FatP * [1, 0, 0]
FatP * [0, 1, 0]
FatP * [1, 1, 0]
gmsh.finalize()
```
