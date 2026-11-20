using LowLevelFEM, LinearAlgebra, Statistics, Printf

"""Make a structured Q4 mesh and optionally move its nodes."""
function q4_grid(nx, ny; lx=1.0, ly=1.0, map=nothing)
    structured_rect_mesh(; lx, ly, nx, ny, order=1)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    if map !== nothing
        for (tag, xyz) in zip(tags, eachcol(reshape(coords, 3, :)))
            x, y = map(xyz[1], xyz[2])
            gmsh.model.mesh.setNode(tag, [x, y, 0.0], Float64[])
        end
    end
    U = Field([Material("body")]; type=:VectorField, dim=2, fieldName=:u)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    xy = reshape(coords, 3, :)[1:2, sortperm(tags)]
    return U, xy
end

"""Assemble standard and selective reduced-integration Q4 stiffness matrices."""
function q4_stiffness(U, E, ν)
    μ = E / (2(1 + ν))
    κ = E / (3(1 - 2ν))
    Ddev = μ * [4/3 -2/3 0; -2/3 4/3 0; 0 0 1]
    Kdev = ∫(SymGrad(U) ⋅ Ddev ⋅ SymGrad(U))
    Kfi  = Kdev + ∫(Div(U) ⋅ κ ⋅ Div(U))
    Ksri = Kdev + ∫(Div(U) ⋅ κ ⋅ Div(U); gauss=:reduced)
    return Kfi, Ksri
end

"""Return the smallest signed Q4 Jacobian at 2×2 Gauss points or at corners."""
function min_jacobian(xy; corners=false)
    a = corners ? 1.0 : inv(sqrt(3.0))
    _, nodes = gmsh.model.mesh.getElementsByType(3)
    points = [-a,-a,0, a,-a,0, a,a,0, -a,a,0]
    _, dN, _ = gmsh.model.mesh.getBasisFunctions(3, points, "GradLagrange")
    dN = reshape(dN, 3, 4, 4)[1:2, :, :]
    return minimum(det(xy[:, tags] * dN[:, :, q]')
                   for tags in eachcol(reshape(nodes, 4, :)), q in 1:4)
end
