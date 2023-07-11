module Testing
    data_dir = joinpath(@__DIR__,"data")
    pdb_file = joinpath(data_dir,"Gromacs","system.pdb")
    trajectory_file = joinpath(data_dir,"Gromacs","trajectory.xtc")
end