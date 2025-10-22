```julia
using LowLevelFEM
gmsh.initialize()

gmsh.open("truss1.geo")
mat = material("rod", E=2e5)
prob = Problem([mat], type=:Truss)
supp = displacementConstraint("left", ux=0)
load1 = load("right", fx=100)
load2 = load("rod", fx=1)

supp0 = displacementConstraint("rod", uy=0, uz=0)
K = stiffnessMatrix(prob)
F = loadVector(prob, [load1])
F.a
fx = loadVector(prob, [load2])
fx.a
K1, f1 = applyBoundaryConditions(K, F + fx, [supp, supp0]);
K1
f1.a
q = K1 \ f1
q.a
N = solveAxialForce(q)
N.A
ux = showDoFResults(q, :ux)
N0 = showElementResults(N, :scalar)
plotOnPath("rod", ux)
plotOnPath("rod", N0)
openPostProcessor()
gmsh.finalize()
```
