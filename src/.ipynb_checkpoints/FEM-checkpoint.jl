module FEM

using LinearAlgebra, SparseArrays
using Arpack
import gmsh

struct Problem
    name::String
    type::String
    dim::Int64
    E::Float64
    ν::Float64
    ρ::Float64
    thickness::Float64
    Problem(name; E=2e5,ν=0.3,ρ=7.85e-9,thickness=1,dim=3,type="Solid") = new(name,type,dim,E,ν,ρ,thickness)
end

function displacementConstraint(name; ux=1im, uy=1im, uz=1im)
    bc0 = name, ux, uy, uz
    return bc0
end

function traction(name; fx=0, fy=0, fz=0)
    ld0 = name, fx, fy, fz
    return ld0
end

#=
function displacementConstraint(groupOfLines; ux=1im, uy=1im)
    phg = gmsh.model.addPhysicalGroup(1, groupOfLines)
    gmsh.model.setPhysicalName(1, phg, "support")
    bc0 = phg, ux, uy
    return bc0
end
=#

#=
function traction(groupOfLines; fx=0, fy=0, thickness=1)
    phg = gmsh.model.addPhysicalGroup(1, groupOfLines)
    gmsh.model.setPhysicalName(1, phg, "load")
    ld0 = phg, fx, fy, thickness
    return ld0
end
=#

"""
Végeselemes felosztás elvégzése
"""
function generateMesh(problem, surf, elemSize; approxOrder=1, algorithm=6, quadrangle=0, internalNodes=0)
    gmsh.model.setCurrent(problem.name)
    # lekérjük az összes csomópontot
    all = gmsh.model.getEntities(0)
    # megadjuk, hogy a csomóponthoz rendelt eleméret mekkora legyen
    gmsh.model.mesh.setSize(all, elemSize)
    # kiválasztjuk a 8-as számú hálózó algoritmust a 2D-s sf1 felülethez
    gmsh.model.mesh.setAlgorithm(2, surf, algorithm)
    # legeneráljuk a hálót a felület kontúrjához (1D-s)
    gmsh.model.mesh.generate(1)
    # legeneráljuk a hálót a felülethez (2D-s)
    gmsh.model.mesh.generate(2)
    # a legenerált háromszög elemeket négyszög elemekké alakítjuk
    if quadrangle
        gmsh.model.mesh.recombine()
    end
    # másodfokú elemekhez:
    # belső csomópontok használata
    if internalNodes
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 0) # 0:belső csomópontokkal 1:csak éleken lévő csomópontokkal
    else
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1) # 0:belső csomópontokkal 1:csak éleken lévő csomópontokkal
    end
    # közelítés fokszáma (1-től 5-ig)
    gmsh.model.mesh.setOrder(approxOrder)
end

# Merevségi mátrix felépítése
function stiffnessMatrixSolid(problem, PhGname; E=1im, ν=1im)
    if E == 1im || ν == 1im
        E = problem.E
        ν = problem.ν
    end
    if problem.dim == 3 && problem.type == "Solid"
        D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0 0 0;
                                        ν 1-ν ν 0 0 0;
                                        ν ν 1-ν 0 0 0;
                                        0 0 0 (1-2ν)/2 0 0;
                                        0 0 0 0 (1-2ν)/2 0;
                                        0 0 0 0 0 (1-2ν)/2]

        dim = 3
        rowsOfB = 6
        b = 1
    elseif problem.dim == 2 && problem.type == "PlaneStress"
        D = E / (1 - ν^2) * [1 ν 0;
                             ν 1 0;
                             0 0 (1-ν)/2] # ÁSF feladat
        dim = 2
        rowsOfB = 3
        b = problem.thickness
    elseif problem.dim == 2 && problem.type == "PlaneStrain"
        D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν 0;
                                        ν 1-ν 0;
                                        0 0 (1-2ν)/2] # SA feladat
        dim = 2
        rowsOfB = 3
        b = 1
    else
        error("stiffnessMatrixSolid: dimension is $(problem.dim), problem type is $(problem.type).")
    end

    # modell kiválasztása
    gmsh.model.setCurrent(problem.name)
    # csomópontok sorszámának lekérése
    #nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(-1, -1)
    # végeselemek típusának, sorszámának és kapcsolati mátrixának (connectivity matrix) lekérése
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    # a lefoglalandó memória méretének becslése
    lengthOfIJV = sum([(div(length(elemNodeTags[i]), length(elemTags[i])) * dim)^2 * length(elemTags[i]) for i in 1:length(elemTags)])
    nn = []
    I = []
    J = []
    V = []
    V = convert(Vector{Float64}, V)
    sizehint!(I, lengthOfIJV)
    sizehint!(J, lengthOfIJV)
    sizehint!(V, lengthOfIJV)
    for i in 1:length(elemTypes)
        et = elemTypes[i]
        elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
        intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order))
        numIntPoints = length(intWeights)
        comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "GradLagrange")
        ∇h = reshape(dfun, :, numIntPoints)
        nnet = zeros(Int, length(elemTags[i]), numNodes)
        invJac = zeros(3, 3numIntPoints)
        Iidx = zeros(Int, numNodes * dim, numNodes * dim)
        Jidx = zeros(Int, numNodes * dim, numNodes * dim)
        for k in 1:numNodes*dim, l in 1:numNodes*dim
            Iidx[k, l] = l
            Jidx[k, l] = k
        end
        ∂h = zeros(dim, numNodes * numIntPoints) # ∂h-t mindig csak felül kellene írni, nem kell újra és újra memóriát foglalni neki.
        B = zeros(rowsOfB * numIntPoints, dim * numNodes) # B-t mindig csak felül kellene írni?
        K1 = zeros(dim * numNodes, dim * numNodes)
        nn2 = zeros(Int, dim * numNodes)
        for j in 1:length(elemTags[i])
            elem = elemTags[i][j]
            for k in 1:numNodes
                nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
            end
            jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
            Jac = reshape(jac, 3, :)
            for k in 1:numIntPoints
                invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
            end
            ∂h *= 0
            for k in 1:numIntPoints, l in 1:numNodes
                ∂h[1:dim, (k-1)*numNodes+l] = invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k] #??????????????????
            end
            B *= 0
            if dim == 2 && rowsOfB == 3
                for k in 1:numIntPoints, l in 1:numNodes
                    B[k*3-0, l*2-0] = B[k*3-2, l*2-1] = ∂h[1, (k-1)*numNodes+l]
                    B[k*3-0, l*2-1] = B[k*3-1, l*2-0] = ∂h[2, (k-1)*numNodes+l]
                end
            elseif dim == 3 && rowsOfB == 6
                for k in 1:numIntPoints, l in 1:numNodes
                    B[k*6-5, l*3-2] = B[k*6-2, l*3-1] = B[k*6-0, l*3-0] = ∂h[1, (k-1)*numNodes+l]
                    B[k*6-4, l*3-1] = B[k*6-2, l*3-2] = B[k*6-1, l*3-0] = ∂h[2, (k-1)*numNodes+l]
                    B[k*6-3, l*3-0] = B[k*6-1, l*3-1] = B[k*6-0, l*3-2] = ∂h[3, (k-1)*numNodes+l]
                end
            else
                error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
            end
            K1 *= 0
            for k in 1:numIntPoints
                B1 = B[k*rowsOfB-(rowsOfB-1):k*rowsOfB, 1:dim*numNodes]
                K1 += B1' * D * B1 * jacDet[k] * intWeights[k]
            end
            for k in 1:dim
                nn2[k:dim:dim*numNodes] = dim * nnet[j, 1:numNodes] .- (dim - k)
            end
            #nn2[1:3:3*numNodes] = 3 * nnet[j, 1:numNodes] .- 2
            #nn2[2:3:3*numNodes] = 3 * nnet[j, 1:numNodes] .- 1
            #nn2[3:3:3*numNodes] = 3 * nnet[j, 1:numNodes]
            append!(I, nn2[Iidx[:]])
            append!(J, nn2[Jidx[:]])
            append!(V, K1[:])
        end
        push!(nn, nnet)
    end
    K = sparse(I, J, V)
    return K
end

# Tömeg mátrix felépítése
function massMatrixPlaneStress(problem, PhGname; ρ=1im)
    if ρ == 1im
        ρ = problem.ρ
    end
    if problem.dim == 3 && problem.type == "Solid"
        dim = 3
        rowsOfH = 3
        b = 1
    elseif problem.dim == 2 && problem.type == "PlaneStress"
        dim = 2
        rowsOfH = 2
        b = problem.thickness
    elseif problem.dim == 2 && problem.type == "PlaneStrain"
        dim = 2
        rowsOfH = 2
        b = 1
    else
        error("stiffnessMatrixSolid: dimension is $(problem.dim), problem type is $(problem.type).")
    end

    gmsh.model.setCurrent(problem.name)
    # csomópontok sorszámának lekérése
    #nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(-1, -1)
    # végeselemek típusának, sorszámának és kapcsolati mátrixának (connectivity matrix) lekérése
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    # a lefoglalandó memória méretének becslése
    lengthOfIJV = sum([(div(length(elemNodeTags[i]), length(elemTags[i])) * dim)^2 * length(elemTags[i]) for i in 1:length(elemTags)])
    nn = []
    I = []
    J = []
    V = [] # Ezt vajon nem kellene átnevezni másnak? Ebben voltak a K elemei is... A 'sparse' parancs készített róla másolatot?
    V = convert(Vector{Float64}, V)
    sizehint!(I, lengthOfIJV)
    sizehint!(J, lengthOfIJV)
    sizehint!(V, lengthOfIJV)
    for i in 1:length(elemTypes)
        et = elemTypes[i]
        elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
        intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(et, "Gauss" * string(2order))
        numIntPoints = length(intWeights)
        comp, fun, ori = gmsh.model.mesh.getBasisFunctions(et, intPoints, "Lagrange")
        h = reshape(fun, :, numIntPoints)
        nnet = zeros(Int, length(elemTags[i]), numNodes)
        Iidx = zeros(Int, numNodes * dim, numNodes * dim)
        Jidx = zeros(Int, numNodes * dim, numNodes * dim)
        for k in 1:numNodes*dim, l in 1:numNodes*dim
            Iidx[k, l] = l
            Jidx[k, l] = k
        end
        nn2 = zeros(Int, 2 * numNodes)
        H = zeros(rowsOfH * numIntPoints, dim * numNodes)
        for k in 1:numIntPoints, l in 1:numNodes
            for kk in 1:dim
                H[k*dim-(dim-kk), l*2-(dim-kk)] = h[(k-1)*numNodes+l]
            end
        end
        M1 = zeros(dim * numNodes, dim * numNodes)
        for j in 1:length(elemTags[i])
            elem = elemTags[i][j]
            for k in 1:numNodes
                nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
            end
            jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
            M1 *= 0
            for k in 1:numIntPoints
                H1 = H[k*dim-(dim-1):k*dim, 1:dim*numNodes]
                M1 += H1' * H1 * jacDet[k] * intWeights[k]
            end
            M1 *= ρ
            for k in 1:dim
                nn2[k:dim:dim*numNodes] = dim * nnet[j, 1:numNodes] .- (dim - k)
            end
            #nn2[1:2:2*numNodes] = 2 * nnet[j, 1:numNodes] .- 1
            #nn2[2:2:2*numNodes] = 2 * nnet[j, 1:numNodes]
            append!(I, nn2[Iidx[:]])
            append!(J, nn2[Jidx[:]])
            append!(V, M1[:])
        end
        push!(nn, nnet)
    end
    M = sparse(I, J, V)
    M = spdiagm(vec(sum(M, dims=2))) # lumped mass matrix
    return M
end

function applyBoundaryConditions!(problem, stiffMat, supports, tractions)
    dof, dof = size(stiffMat)
    massMat = spzeros(dof, dof)
    dampMat = spzeros(dof, dof)
    stiffMat0, massMat, dampMat, fp = applyBoundaryConditions!(problem, stiffMat, massMat, dampMat, supports, tractions)
    massMat = []
    dampMat = []
    return stiffMat0, fp
end

function getTagForPhysicalName(name)
    dimTags = gmsh.model.getPhysicalGroups(-1)
    i = 1
    while gmsh.model.getPhysicalName(dimTags[i][1], dimTags[i][2]) != name
        i += 1
        if i > length(dimTags)
            error("Physical name '$name' does not exist.")
        end
    end
    #display("$name $(dimTags[i])")
    return dimTags[i][2]
end

function applyBoundaryConditions!(problem, stiffMat, massMat, dampMat, supports, tractions)
    gmsh.model.setCurrent(problem.name)
    dof, dof = size(stiffMat)
    pdim = problem.dim
    fp = zeros(dof)
    for n in 1:length(tractions)
        name, fx, fy, fz = tractions[n]
        if problem.dim == 3
            f = [fx, fy, fz]
        elseif problem.dim == 2
            f = [fx, fy]
        else
            error("applyBoundaryConditions: dimension of the problem is $(problem.dim).")
        end
        #tags = gmsh.model.getEntitiesForPhysicalGroup(2, phg)
        dimTags = gmsh.model.getEntitiesForPhysicalName(name)
        for i ∈ 1:length(dimTags)
            dimTag = dimTags[i]
            dim = dimTag[1]
            tag = dimTag[2]
            elementTypes, elementTags, elemNodeTags = gmsh.model.mesh.getElements(dim, tag)
            for ii in 1:length(elementTypes)
                elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elementTypes[ii])
                nnoe = reshape(elemNodeTags[ii], numNodes, :)'
                intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(elementTypes[ii], "Gauss" * string(order))
                numIntPoints = length(intWeights)
                comp, fun, ori = gmsh.model.mesh.getBasisFunctions(elementTypes[ii], intPoints, "Lagrange")
                h = reshape(fun, :, numIntPoints)
                #comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(elementTypes[1], intPoints, "GradLagrange")
                #∇h = reshape(dfun, :, numIntPoints)
                H = zeros(pdim * numIntPoints, pdim * numNodes)
                for j in 1:numIntPoints
                    for k in 1:numNodes
                        for l in 1:pdim
                            H[j*pdim-(pdim-l), k*pdim-(pdim-l)] = h[k, j]
                            #H[j*3-2, k*3-2] = h[k, j]
                            #H[j*3-1, k*3-1] = h[k, j]
                            #H[j*3-0, k*3-0] = h[k, j]
                        end
                    end
                end
                f1 = zeros(pdim * numNodes)
                nn2 = zeros(Int, pdim * numNodes)
                for l in 1:length(elementTags[ii])
                    elem = elementTags[ii][l]
                    jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                    Jac = reshape(jac, 3, :)
                    f1 *= 0
                    for j in 1:numIntPoints
                        H1 = H[j*pdim-(pdim-1):j*pdim, 1:pdim*numNodes]
                        ############### NANSON ###########################################
                        if pdim == 3 && dim == 2
                            xy = Jac[1, 3*j-2] * Jac[2, 3*j-1] - Jac[2, 3*j-2] * Jac[1, 3*j-1]
                            yz = Jac[2, 3*j-2] * Jac[3, 3*j-1] - Jac[3, 3*j-2] * Jac[2, 3*j-1]
                            zx = Jac[3, 3*j-2] * Jac[1, 3*j-1] - Jac[1, 3*j-2] * Jac[3, 3*j-1]
                            Ja = √(xy^2 + yz^2 + zx^2)
                        elseif pdim == 2 && dim == 1
                            Ja = √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2)
                            # Ide lehetne rakni még a térfogati terhelést is?
                        else
                            error("applyBoundaryConditions: dimension of the problem is $(problem.dim), dimension of load is $dim.")
                        end
                        f1 += H1' * f * Ja * intWeights[j]
                    end
                    for k in 1:pdim
                        nn2[k:pdim:pdim*numNodes] = pdim * nnoe[l, 1:numNodes] .- (pdim - k)
                        #nn2[1:3:3*numNodes] = 3 * nnoe[l, 1:numNodes] .- 2
                        #nn2[2:3:3*numNodes] = 3 * nnoe[l, 1:numNodes] .- 1
                        #nn2[3:3:3*numNodes] = 3 * nnoe[l, 1:numNodes]
                    end
                    fp[nn2] += f1
                end
            end
        end
    end

    for i in 1:length(supports)
        name, ux, uy, uz = supports[i]
        #display(name)
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        #phg, ux, uy = supports[i]
        #nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(1, phg)
        if ux != 1im
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= 3
            nodeTagsX .-= 2
            f0 = spzeros(dof, length(nodeTagsX))
            f0 = stiffMat[:, nodeTagsX] * ux
            f0 = sum(f0, dims=2)
            fp -= f0
        end
        if uy != 1im
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= 3
            nodeTagsX .-= 1
            f0 = spzeros(dof, length(nodeTagsX))
            f0 = stiffMat[:, nodeTagsX] * uy
            f0 = sum(f0, dims=2)
            fp -= f0
        end
        if pdim == 3 && uz != 1im
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= 3
            f0 = spzeros(dof, length(nodeTagsY))
            f0 = stiffMat[:, nodeTagsY] * uz
            f0 = sum(f0, dims=2)
            fp -= f0
        end
    end

    for i in 1:length(supports)
        #phg, ux, uy = supports[i]
        #nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(1, phg)
        name, ux, uy, uz = supports[i]
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(-1, phg)
        if ux != 1im
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim-1)
            for j ∈ nodeTagsX
                stiffMat[j, :] .= 0
                stiffMat[:, j] .= 0
                stiffMat[j, j] = 1
                massMat[j, :] .= 0
                massMat[:, j] .= 0
                massMat[j, j] = 1
                dampMat[j, :] .= 0
                dampMat[:, j] .= 0
                dampMat[j, j] = 1
                fp[j] = ux
            end
        end
        if uy != 1im
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= pdim
            nodeTagsX .-= (pdim-2)
            for j ∈ nodeTagsX
                stiffMat[j, :] .= 0
                stiffMat[:, j] .= 0
                stiffMat[j, j] = 1
                massMat[j, :] .= 0
                massMat[:, j] .= 0
                massMat[j, j] = 1
                dampMat[j, :] .= 0
                dampMat[:, j] .= 0
                dampMat[j, j] = 1
                fp[j] = uy
            end
        end
        if pdim == 3 && uz != 1im
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= 3
            for j ∈ nodeTagsY
                stiffMat[j, :] .= 0
                stiffMat[:, j] .= 0
                stiffMat[j, j] = 1
                massMat[j, :] .= 0
                massMat[:, j] .= 0
                massMat[j, j] = 1
                dampMat[j, :] .= 0
                dampMat[:, j] .= 0
                dampMat[j, j] = 1
                fp[j] = uz
            end
        end
    end

    dropzeros!(stiffMat)
    return stiffMat, massMat, dampMat, fp
end

function smallestPeriodTime(K, M)
    #using SymRCM
    #perm = symrcm(K)
    #Kp = K[perm, perm]
    #Mp = M[perm, perm]
    #ω², ϕ = Arpack.eigs(Kp, Mp, nev=1, which=:LM)
    #α = 1e10
    #Ks = K - α * M
    ω², ϕ = Arpack.eigs(K, M, nev=1, which=:LM)

    #if norm(Kp * ϕ - λ²[1] * Mp * ϕ) > 1e-6
    err = norm(K * ϕ[:,1] - ω²[1] * M * ϕ[:,1])
    if err > 1#1e-6
        error("Túl nagy a hiba a legnagyobb sajátérték számításánál: $err")
    end
    Δt = 2π / √(real(ω²[1]))
    return Δt
end

function CDM(K, M, C, T, Δt)
    invM = spdiagm(1 ./ diag(M))
    nsteps = ceil(Int64, T / Δt)
    dof, dof = size(K)

    u = zeros(dof, nsteps)
    v = zeros(dof, nsteps)
    #p = zeros(nsteps)
    t = zeros(nsteps)
    kene = zeros(nsteps)
    sene = zeros(nsteps)
    diss = zeros(nsteps)

    f = zeros(dof)
    u0 = zeros(dof)
    v0 = zeros(dof)
    v0[1:2:dof] .= 1000

    nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(1, 1) #dimTag
    nodeTags *= 2
    nodeTags .-= 1
    v0[nodeTags] .= 0
    a0 = M \ (f - K * u0 - C * v0)
    u00 = u0 - v0 * Δt + a0 * Δt^2 / 2

    u[:, 1] = u0
    v[:, 1] = v0
    t[1] = 0
    kene[1] = dot(v0' * M, v0) / 2
    sene[1] = dot(u0' * K, u0) / 2

    for i in 2:nsteps
        u1 = 2.0 * u0 - u00 + Δt * Δt * invM * (f - K * u0) - Δt * invM * (C * (u0 - u00))
        u[:, i] = u1
        v1 = (u1 - u0) / Δt
        v[:, i] = v1
        t[i] = t[i-1] + Δt
        kene[i] = dot(v1' * M, v1) / 2
        sene[i] = dot(u1' * K, u1) / 2
        #diss[i] = dot(v1' * C, v1)
        u00 = u0
        u0 = u1
    end
    return u, v, t
end

function solveDisplacement(K, f)
    return K \ f
end

# Feszültségek számítása
function solveStress(problem, q)
    E = problem.E
    ν = problem.ν
    if problem.dim == 3 && problem.type == "Solid"
        D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0 0 0;
                                        ν 1-ν ν 0 0 0;
                                        ν ν 1-ν 0 0 0;
                                        0 0 0 (1-2ν)/2 0 0;
                                        0 0 0 0 (1-2ν)/2 0;
                                        0 0 0 0 0 (1-2ν)/2]

        dim = 3
        rowsOfB = 6
        b = 1
    elseif problem.dim == 2 && problem.type == "PlaneStress"
        D = E / (1 - ν^2) * [1 ν 0;
                             ν 1 0;
                             0 0 (1-ν)/2] # ÁSF feladat
        dim = 2
        rowsOfB = 3
        b = problem.thickness
    elseif problem.dim == 2 && problem.type == "PlaneStrain"
        D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν 0;
                                        ν 1-ν 0;
                                        0 0 (1-2ν)/2] # SA feladat
        dim = 2
        rowsOfB = 3
        b = 1
    else
        error("stiffnessMatrixSolid: dimension is $(problem.dim), problem type is $(problem.type).")
    end

    gmsh.model.setCurrent(problem.name)
    # csomópontok sorszámának lekérése
    #nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(-1, -1)
    # végeselemek típusának, sorszámának és kapcsolati mátrixának (connectivity matrix) lekérése
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    numElem = Int[]
    σ = Vector{Float64}[]
    #σx = Vector{Float64}[]
    #σy = Vector{Float64}[]
    #σxy = Vector{Float64}[]
    for i in 1:length(elemTypes)
        et = elemTypes[i]
        elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
        s0 = zeros(rowsOfB * numNodes) # csak SA és ÁSF feladatnál, FSZ-nél már 4 kell
        nodeCoord = zeros(numNodes * 3)
        #display("size of nodeCoord: $(size(nodeCoord))")
        #display("size of localNodeCoord: $(size(localNodeCoord))")
        #display("elemType = $et")
        #display("dim = $dim")
        #display("numNodes = $numNodes")
        for k in 1:dim, j = 1:numNodes
            nodeCoord[k+(j-1)*3] = localNodeCoord[k+(j-1)*dim]
        end
        comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "GradLagrange")
        ∇h = reshape(dfun, :, numNodes)
        nnet = zeros(Int, length(elemTags[i]), numNodes)
        invJac = zeros(3, 3numNodes)
        ∂h = zeros(3, numNodes * numNodes)
        B = zeros(rowsOfB * numNodes, dim * numNodes)
        nn2 = zeros(Int, dim * numNodes)
        for j in 1:length(elemTags[i])
            elem = elemTags[i][j]
            for k in 1:numNodes
                nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
            end
            jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, nodeCoord)
            #display("size of jac: $(size(jac))")
            Jac = reshape(jac, 3, :)
            for k in 1:numNodes
                invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
            end
            ∂h *= 0
            for k in 1:numNodes, l in 1:numNodes
                ∂h[1:dim, (k-1)*numNodes+l] = invJac[1:dim, k*3-2:k*3-(3-dim)] * ∇h[l*3-2:l*3-(3-dim), k] #??????????????????
            end
            B *= 0
            if dim == 2 && rowsOfB == 3
                for k in 1:numNodes, l in 1:numNodes
                    B[k*3-0, l*2-0] = B[k*3-2, l*2-1] = ∂h[1, (k-1)*numNodes+l]
                    B[k*3-0, l*2-1] = B[k*3-1, l*2-0] = ∂h[2, (k-1)*numNodes+l]
                end
            elseif dim == 3 && rowsOfB == 6
                for k in 1:numNodes, l in 1:numNodes
                    B[k*6-5, l*3-2] = B[k*6-2, l*3-1] = B[k*6-0, l*3-0] = ∂h[1, (k-1)*numNodes+l]
                    B[k*6-4, l*3-1] = B[k*6-2, l*3-2] = B[k*6-1, l*3-0] = ∂h[2, (k-1)*numNodes+l]
                    B[k*6-3, l*3-0] = B[k*6-1, l*3-1] = B[k*6-0, l*3-2] = ∂h[3, (k-1)*numNodes+l]
                end
            else
                error("stiffnessMatrix: rows of B is $rowsOfB, dimension of the problem is $dim.")
            end
            push!(numElem, elem)
            for k in 1:dim
                nn2[k:dim:dim*numNodes] = dim * nnet[j, 1:numNodes] .- (dim - k)
            end
            #nn2[1:3:3*numNodes] = 3 * nnet[j, 1:numNodes] .- 2
            #nn2[2:3:3*numNodes] = 3 * nnet[j, 1:numNodes] .- 1
            #nn2[3:3:3*numNodes] = 3 * nnet[j, 1:numNodes]
            #nn2[1:2:2*numNodes] = 2 * nnet[j, 1:numNodes] .- 1
            #nn2[2:2:2*numNodes] = 2 * nnet[j, 1:numNodes]
            s = zeros(9numNodes) # tenzornak 9 eleme van
            #sx = zeros(numNodes)
            #sy = zeros(numNodes)
            #sxy = zeros(numNodes)
            for k in 1:numNodes
                if rowsOfB == 6 && dim == 3
                    B1 = B[k*6-5:k*6, 1:3*numNodes]
                    s0 = D * B1 * q[nn2]
                    s[(k-1)*9+1:k*9] = [s0[1], s0[4], s0[6],
                                        s0[4], s0[2], s0[5],
                                        s0[6], s0[5], s0[3]]
                elseif rowsOfB == 3 && dim == 2 && problem.type == "PlaneStress"
                    B1 = B[k*3-2:k*3, 1:2*numNodes]
                    s0 = D * B1 * q[nn2]
                    s[(k-1)*9+1:k*9] = [s0[1], s0[3], 0,
                                        s0[3], s0[2], 0,
                                         0   ,  0   , 0]
                else
                    error("solveStress: rowsOfB is $rowsOfB, dimension of the problem is $dim, problem type is $(problem.type).")
                end
            end
            push!(σ, s)
        end
    end
    return σ, numElem
end

function showResultUvec(problem, q; name="uvec", visible=false)
    gmsh.model.setCurrent(problem.name)
    dim = problem.dim
    #display("dim = $dim")
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(-1, -1)
    non = length(nodeTags)
    uvec = gmsh.view.add(name)
    u = zeros(3 * non)
    for i in 1:length(nodeTags)
        u[3i-2] = q[dim*nodeTags[i]-(dim-1)]
        u[3i-1] = q[dim*nodeTags[i]-(dim-2)]
        u[3i-0] = dim == 3 ? q[dim*nodeTags[i]-0] : 0
    end
    gmsh.view.addHomogeneousModelData(uvec, 0, problem.name, "NodeData", nodeTags, u, 0, 3)

    gmsh.view.option.setNumber(uvec, "AdaptVisualizationGrid", 1)
    gmsh.view.option.setNumber(uvec, "TargetError", -1e-4)
    gmsh.view.option.setNumber(uvec, "MaxRecursionLevel", 1) # order + 1
    if visible == false
        gmsh.view.option.setNumber(uvec, "Visible", 0)
    end
    return uvec
end

function showResultUX(problem, q; name="ux", visible=false)
    gmsh.model.setCurrent(problem.name)
    dim = problem.dim
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(-1, -1)
    non = length(nodeTags)
    ux = gmsh.view.add(name)
    u = zeros(non)
    for i in 1:length(nodeTags)
        u[i] = q[dim*nodeTags[i]-(dim-1)]
    end
    gmsh.view.addHomogeneousModelData(ux, 0, problem.name, "NodeData", nodeTags, u, 0, 1)

    gmsh.view.option.setNumber(ux, "AdaptVisualizationGrid", 1)
    gmsh.view.option.setNumber(ux, "TargetError", -1e-4)
    gmsh.view.option.setNumber(ux, "MaxRecursionLevel", 1)
    if visible == false
        gmsh.view.option.setNumber(ux, "Visible", 0)
    end
    return ux
end

function showResultUY(problem, q; name="uy", visible=false)
    gmsh.model.setCurrent(problem.name)
    dim = problem.dim
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(-1, -1)
    non = length(nodeTags)
    uy = gmsh.view.add(name)
    u = zeros(non)
    for i in 1:length(nodeTags)
        u[i] = q[dim*nodeTags[i]-(dim-2)]
    end
    gmsh.view.addHomogeneousModelData(uy, 0, problem.name, "NodeData", nodeTags, u, 0, 1)

    gmsh.view.option.setNumber(uy, "AdaptVisualizationGrid", 1)
    gmsh.view.option.setNumber(uy, "TargetError", -1e-4)
    gmsh.view.option.setNumber(uy, "MaxRecursionLevel", 1)
    if visible == false
        gmsh.view.option.setNumber(uy, "Visible", 0)
    end
    return uy
end

function showResultUZ(problem, q; name="uz", visible=false)
    gmsh.model.setCurrent(problem.name)
    dim = problem.dim
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(-1, -1)
    non = length(nodeTags)
    uz = gmsh.view.add(name)
    u = zeros(non)
    for i in 1:length(nodeTags)
        u[i] = dim == 3 ? q[dim*nodeTags[i]-0] : 0
    end
    gmsh.view.addHomogeneousModelData(uz, 0, problem.name, "NodeData", nodeTags, u, 0, 1)

    gmsh.view.option.setNumber(uz, "AdaptVisualizationGrid", 1)
    gmsh.view.option.setNumber(uz, "TargetError", -1e-4)
    gmsh.view.option.setNumber(uz, "MaxRecursionLevel", 1)
    if visible == false
        gmsh.view.option.setNumber(uz, "Visible", 0)
    end
    return uz
end

function showStressResults(problem, S, comp; name="σ", visible=false, smooth=true)
    gmsh.model.setCurrent(problem.name)
    dim = problem.dim
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    σ, numElem = S
    S = gmsh.view.add(name)
    if comp == "s"
        σcomp = σ
        nc =9
    end
    k = 1im
    if comp != "s"
        nc = 1
        if comp == "sx"
            k = 8
        elseif comp == "sy"
            k = 4
        elseif comp == "sz"
            k = 0
        elseif comp == "sxy" || comp == "syx"
            k = 7
        elseif comp == "syz" || comp == "szy"
            k = 3
        elseif comp == "szx" || comp == "sxz"
            k = 6
        else
            error("ShowStressResults: component is $comp ????")
        end
        σcomp = []
        sizehint!(σcomp, length(numElem))
        for i in 1:length(σ)
            sx = zeros(div(length(σ[i]), 9))
            for j in 1:(div(length(σ[i]), 9))
                sx[j] = σ[i][9j-k]
            end
            push!(σcomp, sx)
        end
    end

    gmsh.view.addModelData(S, 0, problem.name, "ElementNodeData", numElem, σcomp, 0, nc)

    if smooth == true
        gmsh.plugin.setNumber("Smooth", "View", -1)
        gmsh.plugin.run("Smooth")
    end

    gmsh.view.option.setNumber(S, "AdaptVisualizationGrid", 1)
    gmsh.view.option.setNumber(S, "TargetError", -1e-4)
    gmsh.view.option.setNumber(S, "MaxRecursionLevel", 1)
    if visible == false
        gmsh.view.option.setNumber(S, "Visible", 0)
    end
    display("$comp..ok")
    return S
end

function plotOnPath(problem, pathName, field, points; numOfStep=0, name="path", visible=false)
    gmsh.model.setCurrent(problem.name)
    dimTags = gmsh.model.getEntitiesForPhysicalName(pathName)
    i = 1
    while dimTags[i][1] != 1
        i += 1
        if i > length(dimTags)
            error("Physical name '$name' with dimension ONE does not exist.")
        end
    end
    path = dimTags[i][2]
    dataType, tags, data, time, numComponents = gmsh.view.getModelData(field, numOfStep)
    bounds = gmsh.model.getParametrizationBounds(1, path)
    bound1 = bounds[1][1]
    bound2 = bounds[2][1]
    step0 = (bound2 - bound1) / (points - 1)
    cv = zeros(4)
    CoordValue = []
    pt0 = gmsh.model.getValue(1, path, [bound1])
    for i in 1:points
        pt1 = gmsh.model.getValue(1, path, [bound1 + (i - 1) * step0])
        cv[1:3] = pt1 - pt0
        val, dis = gmsh.view.probe(field, pt1[1], pt1[2], pt1[3])
        if dis == 0
            if numComponents == 1
                v = val[1]
            elseif numComponents == 3
                v = √(val[1]^2 + val[1]^2 + val[1]^2)
            elseif numComponents == 9
                v = √(0.5 * ((val[1] - val[5])^2 + (val[5] - val[9])^2 + (val[9] - val[1])^2 + 6 * (val[2]^2 + val[3]^2 + val[6]^2)))
            else
                error("Vagy skalás vagy vektor vagy tenzor...")
            end
        else
            v = 0
        end
        cv[4] = v
        append!(CoordValue, cv)
    end
    pathView = gmsh.view.add(name)
    gmsh.view.addListData(pathView, "SP", points, CoordValue)

    gmsh.view.option.setNumber(pathView, "Type", 2)
    gmsh.view.option.setNumber(pathView, "Axes", 1)

    if visible == false
        gmsh.view.option.setNumber(pathView, "Visible", 0)
    end
    return pathView
end

function showResultsVTvec(problem, v, t; name="v(t)", visible=false)
    gmsh.model.setCurrent(problem.name)

    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(2, -1)
    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(-1, -1)
    non = length(nodeTags)
    vvec = gmsh.view.add(name)
    gmsh.option.setNumber("View.DisplacementFactor", 0)
    vv = zeros(3 * non)
    dof, nsteps = size(v)
    for j in 1:nsteps
        for i in 1:length(nodeTags)
            vv[3i-2] = v[2*nodeTags[i]-1, j]
            vv[3i-1] = v[2*nodeTags[i]-0, j]
        end
        gmsh.view.addHomogeneousModelData(vvec, j - 1, "rectangle", "NodeData", nodeTags, vv, t[j], 3)
    end
    
    gmsh.view.option.setNumber(vvec, "NormalRaise", 0.03)
    gmsh.view.option.setNumber(vvec, "DisplacementFactor", 0)
    gmsh.view.option.setNumber(vvec, "AdaptVisualizationGrid", 1)
    gmsh.view.option.setNumber(vvec, "TargetError", -1e-4)
    gmsh.view.option.setNumber(vvec, "MaxRecursionLevel", 1)
    if visible == false
        gmsh.view.option.setNumber(vvec, "Visible", 0)
    end
    return vvec
end

end #module