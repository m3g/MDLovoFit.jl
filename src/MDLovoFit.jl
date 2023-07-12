module MDLovoFit

using TestItems
using DelimitedFiles
using PDBTools
import Chemfiles
import MDLovoFit_jll

# Structure that will contain the output of MDLovoFit
struct MDLovoFitResult
    iframe::Vector{Int}
    rmsd_low::Vector{Float64}
    rmsd_high::Vector{Float64}
    rmsd_all::Vector{Float64}
    rmsf::AbstractVector{Float64}
    aligned_pdb::String
end

"""
    write_frame!(trajectory_pdb_file, atoms, frame)

Append a frame to a PDB file.

"""
function write_frame!(
    trajectory_pdb_file::IO,
    atoms::AbstractVector{<:PDBTools.Atom},
    frame::Chemfiles.Frame;
) 
    coordinates = Chemfiles.positions(frame)
    for (i,col) in enumerate(eachcol(coordinates))
        iatom = findfirst(==(i), atom.index for atom in atoms)
        if !isnothing(iatom)
           atoms[iatom].x = col[1]
           atoms[iatom].y = col[2]
           atoms[iatom].z = col[3]
           println(trajectory_pdb_file, PDBTools.write_atom(atoms[iatom]))
        end
    end
    println(trajectory_pdb_file, "ENDMDL")
    return nothing
end

"""
    write_tmp_pdb_trajectory(atoms, trajectory; first=1, last=length(trajectory), stride=1)

Write a temporary PDB file containing a trajectory. Returns the trajectory file name.

"""
function write_tmp_pdb_trajectory(
    atoms::AbstractVector{<:PDBTools.Atom}, 
    trajectory::Chemfiles.Trajectory; 
    first=1, last=length(trajectory), stride=1
)
    tmp_file = tempname()
    tmp_trajectory = open(tmp_file, "w")
    for (iframe, frame) in enumerate(trajectory)
        if iframe % stride != 0
            continue
        end
        if iframe < first
            continue
        end
        if iframe > last
            break
        end
        write_frame!(tmp_trajectory, atoms, frame)
    end
    close(tmp_trajectory)
    return tmp_file
end

"""
    mdlovofit(atoms, trajectory_file, fraction; output_pdb=nothing, atoms_to_consider=atoms, first=1, last=nothing, stride=1, iref=1)

Run MDLovoFit on a trajectory.

"""
function mdlovofit(
    atoms::AbstractVector{<:PDBTools.Atom}, 
    trajectory_file::String, 
    fraction::AbstractFloat; 
    output_pdb::Union{String,Nothing} = nothing,
    atoms_to_consider::AbstractVector{<:PDBTools.Atom} = atoms,
    first=1, last=nothing, stride=1, iref=1,
)
    trajectory = Chemfiles.Trajectory(trajectory_file)
    if isnothing(last)
        last = length(trajectory)
    else
        if last > length(trajectory)
            throw(ArgumentError("last > length(trajectory)"))
        end
    end
    # Write temporary PDB file with the trajectory
    tmp_trajectory_file = write_tmp_pdb_trajectory(atoms, trajectory; first=first, last=last, stride=stride)
    close(trajectory)
    # Output PDB file
    if isnothing(output_pdb) 
        output_pdb = tempname()
    end
    # Write atoms to consider file, by default all atoms are considered
    tmp_atoms_to_consider_file = tempname()
    for atom in atoms
        i = findfirst(==(atom.index), atom_to_consider.index for atom_to_consider in atoms_to_consider) 
        if !isnothing(i)
            atoms[i].beta = 1.0
        else
            atoms[i].beta = 0.0
        end
    end
    writePDB(atoms, tmp_atoms_to_consider_file)
    # Name of RMSF file
    rmsf_file = tempname()
    # Name of RMSD file
    rmsd_file = tempname()
    # Run MDLovoFit
    try
        MDLovoFit_jll.mdlovofit() do exe
            pipeline(`$exe -f $fraction -iref $iref -rmsf $rmsf_file -t $output_pdb $tmp_trajectory_file`; stdout=rmsd_file)
        end
    catch 
        "ERROR in MDLovoFit execution"
        "Command executed: $command"
    end
    # Read RMSD file
    rmsd_data = readdlm(rmsd_file; comment_char='#')
    frame_index = collect(first:stride:last)
    rmsd_low = rmsd_data[:,2]
    rmsd_high = rmsd_data[:,3]
    rmsd_all = rmsd_data[:,4]
    # Read RMSF file
    rmsf = readdlm(rmsf_file; comment_char='#')
    # If the output PDB is not desired then delete it
    rm(output_pdb)
    # Retun the data structure with the results
    return MDLovoFitResult(frame_index, rmsd_low, rmsd_high, rmsd_all, rmsf, output_pdb)
end

# Testing module
include("../test/Testing.jl")

end # module MDLovoFit
