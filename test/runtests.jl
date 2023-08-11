using TestItemRunner

@testitem "Aqua.test_all" begin
    import Aqua
    Aqua.test_all(MDLovoFit)
end

#@testiem "mdlovofit -f" begin
#    using MDLovoFit.Testing
#    using PDBTools
#    pdbfile = Testing.pdb_file
#    trajectory = Testing.trajectory_file
#    atoms = readPDB(pdbfile, "name CA")
#end
