using TestItemRunner

@testiem "mdlovofit -f" begin
    using MDLovoFit.Testing
    using PDBTools
    pdbfile = Testing.pdb_file
    trajectory = Testing.trajectory_file

    atoms = readPDB(pdbfile, "name CA")




end
