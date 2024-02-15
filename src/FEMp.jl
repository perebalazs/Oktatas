module FEM

using LinearAlgebra, SparseArrays
using Arpack
import gmsh

struct Problem
    name::String
    type::String
    E::Float64
    ν::Float64
    ρ::Float64
    b::Float64
    function Problem(name; E=2e5,ν=0.3,ρ=7.85e-9,b=1,type="PlaneStress")
        gmsh.model.add(name)
        gmsh.model.occ.synchronize()
        return new(name,type,E,ν,ρ,b)
    end
end

function displacementConstraintOnLines(groupOfLines; ux=1im, uy=1im, name="support")
    phg = gmsh.model.addPhysicalGroup(1, groupOfLines)
    gmsh.model.setPhysicalName(1, phg, name)
    bc0 = name, ux, uy
    return bc0
end

function displacementConstraint(name; ux=1im, uy=1im)
    bc0 = name, ux, uy
    return bc0
end

function tractionOnLines(groupOfLines; fx=0, fy=0, thickness=1, name="load")
    phg = gmsh.model.addPhysicalGroup(1, groupOfLines)
    gmsh.model.setPhysicalName(1, phg, name)
    ld0 = name, fx, fy, thickness
    return ld0
end

function traction(name; fx=0, fy=0, thickness=1)
    ld0 = name, fx, fy, thickness
    return ld0
end

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
function stiffnessMatrixPlaneStress(problem)
    # anyagállandók mátrixa
    E = problem.E
    ν = problem.ν
    b = problem.b
    D = E / (1 - ν^2) * [1 ν 0; ν 1 0; 0 0 (1-ν)/2] # ÁSF feladat
    #D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν 0; ν 1-ν 0; 0 0 (1-2ν)/2] # SA feladat
    # modell kiválasztása
    gmsh.model.setCurrent(problem.name)
    # csomópontok sorszámának lekérése
    #nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(-1, -1)
    # végeselemek típusának, sorszámának és kapcsolati mátrixának (connectivity matrix) lekérése
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(2, -1)
    # a lefoglalandó memória méretének becslése
    lengthOfIJV = sum([(div(length(elemNodeTags[i]), length(elemTags[i])) * 2)^2 * length(elemTags[i]) for i in 1:length(elemTags)])
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
        Iidx = zeros(Int, numNodes * 2, numNodes * 2)
        Jidx = zeros(Int, numNodes * 2, numNodes * 2)
        for k in 1:numNodes*2, l in 1:numNodes*2
            Iidx[k, l] = l
            Jidx[k, l] = k
        end
        ∂h = zeros(3, numNodes * numIntPoints) # ∂h-t mindig csak felül kellene írni, nem kell újra és újra memóriát foglalni neki.
        B = zeros(3 * numIntPoints, 2 * numNodes) # B-t mindig csak felül kellene írni?
        K1 = zeros(2 * numNodes, 2 * numNodes)
        nn2 = zeros(Int, 2 * numNodes)
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
                ∂h[1:2, (k-1)*numNodes+l] = invJac[1:2, k*3-2:k*3-1] * ∇h[l*3-2:l*3-1, k] #??????????????????
            end
            B *= 0
            for k in 1:numIntPoints, l in 1:numNodes
                B[k*3-0, l*2-0] = B[k*3-2, l*2-1] = ∂h[1, (k-1)*numNodes+l]
                B[k*3-0, l*2-1] = B[k*3-1, l*2-0] = ∂h[2, (k-1)*numNodes+l]
            end
            K1 *= 0
            for k in 1:numIntPoints
                B1 = B[k*3-2:k*3, 1:2*numNodes]
                K1 += B1' * D * B1 * jacDet[k] * b * intWeights[k]
            end
            nn2[1:2:2*numNodes] = 2 * nnet[j, 1:numNodes] .- 1
            nn2[2:2:2*numNodes] = 2 * nnet[j, 1:numNodes]
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
function massMatrixPlaneStress(problem)
    b = problem.b
    ρ = problem.ρ
    gmsh.model.setCurrent(problem.name)
    # csomópontok sorszámának lekérése
    nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(-1, -1)
    # végeselemek típusának, sorszámának és kapcsolati mátrixának (connectivity matrix) lekérése
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(2, -1)
    # a lefoglalandó memória méretének becslése
    lengthOfIJV = sum([(div(length(elemNodeTags[i]), length(elemTags[i])) * 2)^2 * length(elemTags[i]) for i in 1:length(elemTags)])
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
        Iidx = zeros(Int, numNodes * 2, numNodes * 2)
        Jidx = zeros(Int, numNodes * 2, numNodes * 2)
        for k in 1:numNodes*2, l in 1:numNodes*2
            Iidx[k, l] = l
            Jidx[k, l] = k
        end
        nn2 = zeros(Int, 2 * numNodes)
        H = zeros(2 * numIntPoints, 2 * numNodes)
        for k in 1:numIntPoints, l in 1:numNodes
            H[k*2-1, l*2-1] = H[k*2-0, l*2-0] = h[(k-1)*numNodes+l]
        end
        M1 = zeros(2 * numNodes, 2 * numNodes)
        for j in 1:length(elemTags[i])
            elem = elemTags[i][j]
            for k in 1:numNodes
                nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
            end
            jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
            M1 *= 0
            for k in 1:numIntPoints
                H1 = H[k*2-1:k*2, 1:2*numNodes]
                M1 += H1' * H1 * jacDet[k] * intWeights[k]
            end
            M1 *= ρ
            nn2[1:2:2*numNodes] = 2 * nnet[j, 1:numNodes] .- 1
            nn2[2:2:2*numNodes] = 2 * nnet[j, 1:numNodes]
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
    dimTags = gmsh.model.getPhysicalGroups(1)
    i = 1
    while gmsh.model.getPhysicalName(1, dimTags[i][2]) != name
        i += 1
        if i > length(dimTags)
            error("Physical name '$name' does not exist.")
        end
    end
    return dimTags[i][2]
end


function applyBoundaryConditions!(problem, stiffMat, massMat, dampMat, supports, tractions)
    gmsh.model.setCurrent(problem.name)
    dof, dof = size(stiffMat)
    fp = zeros(dof)
    for n in 1:length(tractions)
        name, fx, fy, b = tractions[n]
        f = [fx, fy]
        dimTags = gmsh.model.getEntitiesForPhysicalName(name)
        for i ∈ 1:length(dimTags)
            dimTag = dimTags[i]
            dim = dimTag[1]
            tag = dimTag[2]
            #display("traction on $tag")
            if dim != 1
                error("Plane Stress: A perem csak egydimenziós lehet!")
            end
            elementTypes, elementTags, elemNodeTags = gmsh.model.mesh.getElements(1, tag)
            if length(elementTypes) != 1
                error("A peremen nem csak egyfajta végeselem van!")
            end
            #display("elementTags: $elementTags")
            elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elementTypes[1])
            nnoe = reshape(elemNodeTags[1], numNodes, :)'
            intPoints, intWeights = gmsh.model.mesh.getIntegrationPoints(elementTypes[1], "Gauss" * string(order))
            numIntPoints = length(intWeights)
            comp, fun, ori = gmsh.model.mesh.getBasisFunctions(elementTypes[1], intPoints, "Lagrange")
            h = reshape(fun, :, numIntPoints)
            #comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(elementTypes[1], intPoints, "GradLagrange")
            #∇h = reshape(dfun, :, numIntPoints)
            H = zeros(2 * numIntPoints, 2 * numNodes)
            for j in 1:numIntPoints
                for k in 1:numNodes
                    H[j*2-1, k*2-1] = h[k, j]
                    H[j*2-0, k*2-0] = h[k, j]
                end
            end
            f1 = zeros(2 * numNodes)
            nn2 = zeros(Int, 2 * numNodes)
            for l in 1:length(elementTags[1])
                elem = elementTags[1][l]
                jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, intPoints)
                Jac = reshape(jac, 3, :)
                f1 *= 0
                for j in 1:numIntPoints
                    H1 = H[j*2-1:j*2, 1:2*numNodes]
                    Ja = √((Jac[1, 3*j-2])^2 + (Jac[2, 3*j-2])^2)
                    f1 += H1' * f * Ja * b * intWeights[j]
                end
                nn2[1:2:2*numNodes] = 2 * nnoe[l, 1:numNodes] .- 1
                nn2[2:2:2*numNodes] = 2 * nnoe[l, 1:numNodes]
                fp[nn2] += f1
            end
        end

    end

    for i in 1:length(supports)
        name, ux, uy = supports[i]
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(1, phg)
        if ux != 1im
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= 2
            nodeTagsX .-= 1
            f0 = spzeros(dof, length(nodeTagsX))
            f0 = stiffMat[:, nodeTagsX] * ux
            f0 = sum(f0, dims=2)
            fp -= f0
        end
        if uy != 1im
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= 2
            f0 = spzeros(dof, length(nodeTagsY))
            f0 = stiffMat[:, nodeTagsY] * uy
            f0 = sum(f0, dims=2)
            fp -= f0
        end
    end

    for i in 1:length(supports)
        name, ux, uy = supports[i]
        phg = getTagForPhysicalName(name)
        nodeTags, coord = gmsh.model.mesh.getNodesForPhysicalGroup(1, phg)
        if ux != 1im
            nodeTagsX = copy(nodeTags)
            nodeTagsX *= 2
            nodeTagsX .-= 1
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
            nodeTagsY = copy(nodeTags)
            nodeTagsY *= 2
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
                fp[j] = uy
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
function solveStressPlaneStress(problem, q)
    E = problem.E
    ν = problem.ν
    D = E / (1 - ν^2) * [1 ν 0; ν 1 0; 0 0 (1-ν)/2] # ÁSF feladat
    gmsh.model.setCurrent(problem.name)
    # csomópontok sorszámának lekérése
    nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(-1, -1)
    # végeselemek típusának, sorszámának és kapcsolati mátrixának (connectivity matrix) lekérése
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(2, -1)
    numElem = Int[]
    σ = Vector{Float64}[]
    #σx = Vector{Float64}[]
    #σy = Vector{Float64}[]
    #σxy = Vector{Float64}[]
    for i in 1:length(elemTypes)
        et = elemTypes[i]
        elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
        s0 = zeros(3numNodes) # csak SA és ÁSF feladatnál, FSZ-nél már 4 kell
        nodeCoord = zeros(numNodes * 3)
        for k in 1:dim, j in 1:numNodes
            nodeCoord[k+(j-1)*3] = localNodeCoord[k+(j-1)*dim]
        end
        comp, dfun, ori = gmsh.model.mesh.getBasisFunctions(et, nodeCoord, "GradLagrange")
        ∇h = reshape(dfun, :, numNodes)
        nnet = zeros(Int, length(elemTags[i]), numNodes)
        invJac = zeros(3, 3numNodes)
        ∂h = zeros(3, numNodes * numNodes)
        B = zeros(3 * numNodes, 2 * numNodes)
        nn2 = zeros(Int, 2 * numNodes)
        for j in 1:length(elemTags[i])
            elem = elemTags[i][j]
            for k in 1:numNodes
                nnet[j, k] = elemNodeTags[i][(j-1)*numNodes+k]
            end
            jac, jacDet, coord = gmsh.model.mesh.getJacobian(elem, nodeCoord)
            Jac = reshape(jac, 3, :)
            for k in 1:numNodes
                invJac[1:3, 3*k-2:3*k] = inv(Jac[1:3, 3*k-2:3*k])'
            end
            ∂h *= 0
            for k in 1:numNodes, l in 1:numNodes
                ∂h[1:2, (k-1)*numNodes+l] = invJac[1:2, k*3-2:k*3-1] * ∇h[l*3-2:l*3-1, k] #??????????????????
            end
            B *= 0
            for k in 1:numNodes, l in 1:numNodes
                B[k*3-0, l*2-0] = B[k*3-2, l*2-1] = ∂h[1, (k-1)*numNodes+l]
                B[k*3-0, l*2-1] = B[k*3-1, l*2-0] = ∂h[2, (k-1)*numNodes+l]
            end
            push!(numElem, elem)
            nn2[1:2:2*numNodes] = 2 * nnet[j, 1:numNodes] .- 1
            nn2[2:2:2*numNodes] = 2 * nnet[j, 1:numNodes]
            s = zeros(9numNodes) # tenzornak 9 eleme van
            #sx = zeros(numNodes)
            #sy = zeros(numNodes)
            #sxy = zeros(numNodes)
            for k in 1:numNodes
                B1 = B[k*3-2:k*3, 1:2*numNodes]
                s0 = D * B1 * q[nn2]
                s[(k-1)*9+1:k*9] = [s0[1], s0[3], 0, s0[3], s0[2], 0, 0, 0, 0]
                #sx[k] = s0[1]
                #sy[k] = s0[2]
                #sxy[k] = s0[3]
            end
            push!(σ, s)
            #push!(σx, sx)
            #push!(σy, sy)
            #push!(σxy, sxy)
        end
    end
    return σ, numElem
end

function showResultUvec(problem, q; name="uvec", visible=false)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(2, -1)
    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(-1, -1)
    non = length(nodeTags)
    uvec = gmsh.view.add(name)
    u = zeros(3 * non)
    for i in 1:length(nodeTags)
        u[3i-2] = q[2*nodeTags[i]-1]
        u[3i-1] = q[2*nodeTags[i]-0]
    end
    gmsh.view.addHomogeneousModelData(uvec, 0, problem.name, "NodeData", nodeTags, u, 0, 3)

    gmsh.view.option.setNumber(uvec, "AdaptVisualizationGrid", 1)
    gmsh.view.option.setNumber(uvec, "TargetError", -1e-4)
    gmsh.view.option.setNumber(uvec, "MaxRecursionLevel", order + 1)
    if visible == false
        gmsh.view.option.setNumber(uvec, "Visible", 0)
    end
    return uvec
end

function showResultUX(problem, q; name="ux", visible=false)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(2, -1)
    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(-1, -1)
    non = length(nodeTags)
    ux = gmsh.view.add(name)
    u = zeros(non)
    for i in 1:length(nodeTags)
        u[i] = q[2*nodeTags[i]-1]
    end
    gmsh.view.addHomogeneousModelData(ux, 0, problem.name, "NodeData", nodeTags, u, 0, 1)

    gmsh.view.option.setNumber(ux, "AdaptVisualizationGrid", 1)
    gmsh.view.option.setNumber(ux, "TargetError", -1e-4)
    gmsh.view.option.setNumber(ux, "MaxRecursionLevel", order + 1)
    if visible == false
        gmsh.view.option.setNumber(ux, "Visible", 0)
    end
    return ux
end

function showResultUY(problem, q; name="uy", visible=false)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(2, -1)
    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(-1, -1)
    non = length(nodeTags)
    uy = gmsh.view.add(name)
    u = zeros(non)
    for i in 1:length(nodeTags)
        u[i] = q[2*nodeTags[i]-0]
    end
    gmsh.view.addHomogeneousModelData(uy, 0, problem.name, "NodeData", nodeTags, u, 0, 1)

    gmsh.view.option.setNumber(uy, "AdaptVisualizationGrid", 1)
    gmsh.view.option.setNumber(uy, "TargetError", -1e-4)
    gmsh.view.option.setNumber(uy, "MaxRecursionLevel", order + 1)
    if visible == false
        gmsh.view.option.setNumber(uy, "Visible", 0)
    end
    return uy
end

function showResultS(problem, S; name="σ", visible=true, smooth=true)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(2, -1)
    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    σ, numElem = S
    S = gmsh.view.add(name)
    gmsh.view.addModelData(S, 0, problem.name, "ElementNodeData", numElem, σ, 0, 9)

    if smooth == true
        gmsh.plugin.setNumber("Smooth", "View", -1)
        gmsh.plugin.run("Smooth")
    end

    gmsh.view.option.setNumber(S, "AdaptVisualizationGrid", 1)
    gmsh.view.option.setNumber(S, "TargetError", -1e-4)
    gmsh.view.option.setNumber(S, "MaxRecursionLevel", order + 1)
    if visible == false
        gmsh.view.option.setNumber(S, "Visible", 0)
    end
    return S
end

function showResultSX(problem, S; name="σx", visible=false, smooth=true)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(2, -1)
    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    σ, numElem = S
    S = gmsh.view.add(name)
    σx = []
    sizehint!(σx, length(numElem))
    for i in 1:length(σ)
        sx = zeros(div(length(σ[i]), 9))
        for j in 1:(div(length(σ[i]), 9))
            sx[j] = σ[i][9j-8]
        end
        push!(σx, sx)
    end

    gmsh.view.addModelData(S, 0, problem.name, "ElementNodeData", numElem, σx, 0, 1)

    if smooth == true
        gmsh.plugin.setNumber("Smooth", "View", -1)
        gmsh.plugin.run("Smooth")
    end

    gmsh.view.option.setNumber(S, "AdaptVisualizationGrid", 1)
    gmsh.view.option.setNumber(S, "TargetError", -1e-4)
    gmsh.view.option.setNumber(S, "MaxRecursionLevel", order + 1)
    if visible == false
        gmsh.view.option.setNumber(S, "Visible", 0)
    end
    return S
end

function showResultSY(problem, S; name="σy", visible=false, smooth=true)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(2, -1)
    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    σ, numElem = S
    S = gmsh.view.add(name)
    σy = []
    sizehint!(σy, length(numElem))
    for i in 1:length(σ)
        sy = zeros(div(length(σ[i]), 9))
        for j in 1:(div(length(σ[i]), 9))
            sy[j] = σ[i][9j-4]
        end
        push!(σy, sy)
    end

    gmsh.view.addModelData(S, 0, problem.name, "ElementNodeData", numElem, σy, 0, 1)

    if smooth == true
        gmsh.plugin.setNumber("Smooth", "View", -1)
        gmsh.plugin.run("Smooth")
    end

    gmsh.view.option.setNumber(S, "AdaptVisualizationGrid", 1)
    gmsh.view.option.setNumber(S, "TargetError", -1e-4)
    gmsh.view.option.setNumber(S, "MaxRecursionLevel", order + 1)
    if visible == false
        gmsh.view.option.setNumber(S, "Visible", 0)
    end
    return S
end

function showResultSXY(problem, S; name="τxy", visible=false, smooth=true)
    gmsh.model.setCurrent(problem.name)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(2, -1)
    elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    σ, numElem = S
    S = gmsh.view.add(name)
    σxy = []
    sizehint!(σxy, length(numElem))
    for i in 1:length(σ)
        sxy = zeros(div(length(σ[i]), 9))
        for j in 1:(div(length(σ[i]), 9))
            sxy[j] = σ[i][9j-7]
        end
        push!(σxy, sxy)
    end

    gmsh.view.addModelData(S, 0, problem.name, "ElementNodeData", numElem, σxy, 0, 1)

    if smooth == true
        gmsh.plugin.setNumber("Smooth", "View", -1)
        gmsh.plugin.run("Smooth")
    end

    gmsh.view.option.setNumber(S, "AdaptVisualizationGrid", 1)
    gmsh.view.option.setNumber(S, "TargetError", -1e-4)
    gmsh.view.option.setNumber(S, "MaxRecursionLevel", order + 1)
    if visible == false
        gmsh.view.option.setNumber(S, "Visible", 0)
    end
    return S
end

function plotOnPath(problem, path, field, points; numOfStep=0, name="path", visible=false)
    gmsh.model.setCurrent(problem.name)
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
    gmsh.view.option.setNumber(vvec, "MaxRecursionLevel", order + 1)
    if visible == false
        gmsh.view.option.setNumber(vvec, "Visible", 0)
    end
    return vvec
end

end #module