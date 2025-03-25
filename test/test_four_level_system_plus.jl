using TestItems

@testmodule FourLevelSystemPlus begin
    using C13NV.Models: make_nv_system
    using C13NV.Amplitudes: ConstantDrive, LinearChirp
    using C13NV.Units: MHz, kHz, ms

    H, ket, labels = make_nv_system(
        δ₋ = 1.5MHz,
        Ω₋ = ConstantDrive(257kHz),
        ω₋ = LinearChirp(t₀ = 0.1568ms, α = 26.24MHz / ms),
    )
end


@testitem "Four-Level-System (+)" setup = [FourLevelSystemPlus] begin

    using LinearAlgebra: norm
    using RecursiveArrayTools
    using QuantumControl.Controls: get_parameters

    H = FourLevelSystemPlus.H
    parameters = get_parameters(H)
    @test length(parameters) == 3
    @test norm(parameters - [1.0, 1.0, 1.0]) ≈ 0.0
    @test FourLevelSystemPlus.labels ==
          [("G", "-1", "↑"), ("G", "-1", "↓"), ("G", "0", "↑"), ("G", "0", "↓"),]

end


using TestItemRunner
@run_package_tests filter = ti -> ti.filename == @__FILE__
