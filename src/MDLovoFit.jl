module MDLovoFit

using TestItems
using Statistics: mean
using DelimitedFiles
using PDBTools
import Chemfiles
import MDLovoFit_jll

# Functions of the interface
export MDLovoFitResult
export mdlovofit
export MapFractions
export map_fractions

"""
    MapFractionsResult

Structure that will contain the output of map_fractions.

`fraction` contains the fraction of atoms considered in the alignment.

`rmsd_low` contains the RMSD of the fraction of the structure with the lowest RMSD.

`rmsd_high` contains the RMSD of the fraction not considered for the alignment.

`rmsd_all` contains the RMSD of the whole structure.

"""
struct MapFractionsResult
    fraction::Vector{Float64}
    rmsd_low::Vector{Float64}
    rmsd_high::Vector{Float64}
    rmsd_all::Vector{Float64}
end

function show(io::IO, mf::MapFractionsResult)
    print(io, chomp("""
    ------------------
    MapFractionsResult
    ------------------
    Greatest fraction for which the RMSD-low is smaller than 1: $(round(mf.fraction[findlast(<(1.0), mf.rmsd_low)],digits=2))

    fraction contains the fraction of atoms considered in the alignment.
    rmsd_low contains the RMSD of the fraction of the structure with the lowest RMSD.
    rmsd_high contains the RMSD of the fraction not considered for the alignment.
    rmsd_all contains the RMSD of the whole structure.
    """))
    return nothing
end

# Structure that will contain the output of MDLovoFit
"""
    MDLovoFitResult

Structure that will contain the output of MDLovoFit.

`iframe` is the frame index.

`rmsd_low` is the RMSD of the fraction of the structure with the lowest RMSD.

`rmsd_high` is the RMSD of the fraction of the structure with the highest RMSD.

`rmsd_all` is the RMSD of the whole structure.

`rmsf` is the RMSF as a function of the residue or atom index. 

`aligned_pdb` is the name of the PDB file with the aligned structure.

"""
struct MDLovoFitResult
    fraction::Float64
    iframe::Vector{Int}
    rmsd_low::Vector{Float64}
    rmsd_high::Vector{Float64}
    rmsd_all::Vector{Float64}
    rmsf::Vector{Float64}
    aligned_pdb::String
end

import Base: show
function show(io::IO, result::MDLovoFitResult)
    plow = round(100*result.fraction,digits=1)
    phigh = round(100*(1-result.fraction),digits=1)
    av_low = round((mean(result.rmsd_low)), digits=2)
    av_high = round((mean(result.rmsd_high)), digits=2)
    print(io, chomp("""
    ---------------
    MDLovoFitResult
    ---------------

    Aligned pdb file: $(result.aligned_pdb)
    Number of frames considered: $(length(result.iframe))
    Average RMSD of all atoms: $(mean(result.rmsd_all))
    Average RMSD of the $plow% atoms of lowest RMSD: $av_low
    Average RMSD of the $phigh% atoms of highest RMSD: $av_high

    Frame indices availabe in result.iframe
    RMSD data availabe in rmsd_low, rmsd_high, and rmsd_all

    RMSF data availabe in result.rmsf (Number of atoms: $(length(result.rmsf)))
    """))
    return nothing
end

"""
    write_frame!(trajectory_pdb_file, atoms, frame

Write a frame of a trajectory in the temporary PDB trajectory file.

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
    write_tmp_pdb_trajectory(
        atoms::AbstractVector{<:PDBTools.Atom}, 
        trajectory_file::String; 
        first=1, last=nothing, maxframes=100,
    )

Writes a temporary PDB file containing a trajectory. Returns the trajectory file name, the 
last frame to be read, and the `stride` parameter.

Returns the name of the file, the last frame to be read, and the stride parameter,
computed from the `maxframes` input parameter. The `maxframes` parameter defines 
the maximum number of frames (equally spaced) that will be used in the calculation.

"""
function write_tmp_pdb_trajectory(
    atoms::AbstractVector{<:PDBTools.Atom}, 
    trajectory_file::String; 
    first=1, last=nothing, maxframes=100,
)
    Chemfiles.Trajectory(trajectory_file) do trajectory 
        if isnothing(last)
            last = Int(length(trajectory))
        else
            if last > length(trajectory)
                throw(ArgumentError("last > length(trajectory)"))
            end
        end
        # Write temporary PDB file with the trajectory
        stride = max(1,ceil(Int,(last-first+1)/maxframes))
        tmp_file = tempname()
        open(tmp_file, "w") do file_io
            for (iframe, frame) in enumerate(trajectory)
                if iframe < first
                    continue
                end
                if iframe > last
                    break
                end
                if (iframe-first) % stride != 0
                    continue
                end
                write_frame!(file_io, atoms, frame)
            end
        end
        return tmp_file, last, stride
    end
end

"""
    map_fractions(atoms::AbstractVector{<:PDBTools.Atom}, trajectory_file::String)


"""
function map_fractions(
    atoms::AbstractVector{<:PDBTools.Atom}, 
    trajectory_file::String;
    first=1,
    last=nothing,
    maxframes=100,
)
    tmp_trajectory_file, last, _ = 
        write_tmp_pdb_trajectory(atoms, trajectory_file; first=first, last=last, maxframes=maxframes)
    mapfrac_file = tempname()
    MDLovoFit_jll.mdlovofit() do exe
        run(pipeline(`$exe -mapfrac $tmp_trajectory_file`; stdout=mapfrac_file))
    end
    data = readdlm(mapfrac_file, comments=true, comment_char='#')
    range = 1:findlast(<(1), data[:,1])
    fraction = data[range,1]
    rmsd_low = data[range,2]
    rmsd_high = data[range,3]
    rmsd_all = data[range,4]
    return MapFractionsResult(fraction, rmsd_low, rmsd_high, rmsd_all)
end

function map_fractions(
    selection::String, 
    pdb_file::String, 
    trajectory_file::String;
    kargs...
)
    atoms = readPDB(pdb_file, selection)
    return map_fractions(atoms, trajectory_file)
end

"""
    mdlovofit(
        atoms::AbstractVector{<:PDBTools.Atom}, # atoms of the system (usually a "protein and name CA" selection)
        trajectory_file::String, # name of the trajectory file 
        fraction; # fraction of atoms to be aligned
        output_pdb=nothing, # name of the output file. If nothing, no output file is written
        atoms_to_consider=atoms, # may be a different set of atoms
        first=1, 
        last=nothing, # nothing means the last frame of the trajectory
        iref=1, # reference frame
        maxframes=100, # number of frames to be aligned
    )

Run MDLovoFit on a trajectory.

"""
function mdlovofit(
    atoms::AbstractVector{<:PDBTools.Atom}, 
    trajectory_file::String, 
    fraction::AbstractFloat; 
    output_pdb::Union{String,Nothing} = nothing,
    atoms_to_consider::AbstractVector{<:PDBTools.Atom} = atoms,
    first=1, last=nothing, iref=1, maxframes=100,
)
    # Open trajectory and save it in temporary PDB file, for the selected atoms
    tmp_trajectory_file, last, stride = 
        write_tmp_pdb_trajectory(atoms, trajectory_file; first=first, last=last, maxframes=maxframes)
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
        if isnothing(output_pdb)
            MDLovoFit_jll.mdlovofit() do exe
                run(pipeline(`$exe -f $fraction -iref $iref -rmsf $rmsf_file $tmp_trajectory_file`; stdout=rmsd_file))
            end
        else
            MDLovoFit_jll.mdlovofit() do exe
                run(pipeline(`$exe -f $fraction -iref $iref -rmsf $rmsf_file -t $output_pdb $tmp_trajectory_file`; stdout=rmsd_file))
            end
        end
    catch 
        "ERROR in MDLovoFit execution"
        "Command executed: $command"
    end
    # Read RMSD file
    rmsd_data = readdlm(rmsd_file; comments=true, comment_char='#')
    frame_index = collect(first:stride:last)
    rmsd_low = rmsd_data[:,2]
    rmsd_high = rmsd_data[:,3]
    rmsd_all = rmsd_data[:,4]
    # Read RMSF file
    rmsf = readdlm(rmsf_file; comments=true, comment_char='#')[:,2]
    # Retun the data structure with the results
    if isnothing(output_pdb)
        output_pdb = "Not saved."
    end
    return MDLovoFitResult(fraction, frame_index, rmsd_low, rmsd_high, rmsd_all, rmsf, output_pdb)
end

function mdlovofit(
    selection::String,
    pdbfile::String, 
    trajectory_file::String, 
    fraction::AbstractFloat; 
    kargs...
)
    atoms = readPDB(pdbfile, selection)
    return mdlovofit(
        atoms, 
        trajectory_file, 
        fraction; 
        kargs...
    )
end

# Testing module
include("../test/Testing.jl")

end # module MDLovoFit
