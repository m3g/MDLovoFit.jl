using TestItemRunner

@testitem "Aqua.test_all" begin
    import Aqua
    Aqua.test_all(MDLovoFit)
end

@testiem "mdlovofit" begin
    using MDLovoFit
    # PDB file of the system
    pdbfile = MDLovoFit.Testing.pdb_file 
    trajectory = MDLovoFit.Testing.trajectory_file # Gromacs trajecotry
    mf = map_fractions("protein and name CA", pdbfile, trajectory)
    @test sum(mf.fraction) ≈ 49.35
    @test sum(mf.rmsd_all) ≈ 122.41398264
    @test sum(mf.rmsd_high) ≈ 211.04589343599997
    @test sum(mf.rmsd_low) ≈  54.55835694400001

    rm("aligned.pdb"; force=true)
    rm("aligned_rmsf.dat"; force=true)
    rm("aligned_rmsd.dat"; force=true)
    r = mdlovofit(
        "protein and name CA", pdbfile, 
        trajectory, 
        fraction = 0.9, 
        output_pdb="aligned.pdb"
    )
    @test isfile("aligned.pdb")
    @test isfile("aligned_rmsd.dat")
    @test isfile("aligned_rmsf.dat")
    @test r.fraction ≈ 0.9
    @test sum(r.rmsd_high) ≈ 80.349800192
    @test sum(r.rmsd_low) ≈ 24.965490834999997
    @test sum(r.rmsd_all) ≈ 30.795418137
end
