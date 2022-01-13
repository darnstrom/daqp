# Init
work = Libc.malloc(sizeof(DAQP.workspace))
ccall((:allocate_daqp_iters,libdaqp), Cvoid, (Ref{Ptr{Cvoid}},Cint),n);

# Change LDP
unsafe_store!(work,fieldoffset(DAQP.Workspace,3),m);
unsafe_store!(work,fieldoffset(DAQP.Workspace,4),ms);
unsafe_store!(work,fieldoffset(DAQP.Workspace,5),pointer(A));
unsafe_store!(work,fieldoffset(DAQP.Workspace,6),pointer(bupper));
unsafe_store!(work,fieldoffset(DAQP.Workspace,7),pointer(blower));
unsafe_store!(work,fieldoffset(DAQP.Workspace,10),pointer(sense));

# Solve LDP
exitflag = ccall((:daqp_ldp,libdaqp), Cint, (Ref{Ptr{Cvoid}},:),work);


# Cleanup 
ccall((:free_daqp_iters,libdaqp), Cvoid, (Ref{Ptr{Cvoid}},:),work);
Libc.free(work);
