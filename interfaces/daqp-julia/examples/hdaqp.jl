using DAQPBase

As = Matrix{Float64}[]
bus = Vector{Float64}[]
bls = Vector{Float64}[]

push!(As,[1.0 0 0; 0 1 0; 0 0 1])
push!(bus,ones(3))
push!(bls,-ones(3))

push!(As,[1.0 1  1;])
push!(bus,ones(1))
push!(bls,-[1e30])

push!(As,[1.0 -1 0;])
push!(bus,[0.5])
push!(bls,-[0.5])

push!(As,[3.0 1 -1;])
push!(bus,[20.0])
push!(bls,[10.0])

x,es,info = DAQPBase.hidaqp(As,bus,bls)
