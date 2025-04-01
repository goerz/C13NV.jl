module Models

import LinearAlgebra
using QuantumPropagators: hamiltonian, liouvillian
using QuantumPropagators.Generators: dissipator, Generator

using ..Units: Units

function ⊗(A, B)
    LinearAlgebra.kron(A, B)
end


function ZeroMatrix(N)
    return Array{Float64}(zeros(N) * zeros(N)')
end

function IdMatrix(N)
    return Array{Float64}(LinearAlgebra.I(N))
end


function ⊕(A, B)
    N = size(A, 1)
    M = size(B, 1)
    @assert size(A, 2) == N
    @assert size(B, 2) == M
    C = ZeroMatrix(N + M)
    C[1:N, 1:N] .= A[1:N, 1:N]
    C[N+1:N+M, N+1:N+M] .= B[1:M, 1:M]
    return C
end


"""Construct the system

```julia
    H_or_L, ket, labels = make_nv_system(;
        A_zz,
        A_zx,
        B,
        γ_c,
        δ₋,
        δ₊,
        Ω₊,
        Ω₋,
        Γ_E_G,
        Γ_0_E_S,
        Γ_1_E_S,
        Γ_0_S_G,
        Γ_1_S_G,
        use_dissipation = !isnothing(Λ)
    )
```

return a generator, a function `ket` to construct bare basis states, and a
list of labels.
"""
function make_nv_system(;
    A_zz::Float64 = 1.0 * Units.MHz,
    A_zx::Float64 = 0.3 * Units.MHz,
    B::Float64 = 120 * Units.Gauss,
    γ_c::Float64 = 1.07 * Units.kHz / Units.Gauss,
    δ₋::Float64 = 0.0,
    δ₊::Float64 = 0.0,
    ω₊ = nothing,
    ω₋ = nothing,
    Ω₊ = nothing,
    Ω₋ = nothing,
    Λ = nothing, # incoherent optical excitation (proportional to laser power)
    Γ_E_G::Float64 = (1 / (12 * Units.ns)),
    Γ_0_E_S::Float64 = 0.0,
    Γ_1_E_S::Float64 = (1 / (24 * Units.ns + 0.9 * Units.ns)),
    Γ_0_S_G::Float64 = (1 / (219 * Units.ns)),
    Γ_1_S_G::Float64 = (1 / (219 * Units.ns)),
    coherent_optical_drive = false,
    use_dissipation = !isnothing(Λ),
)

    labels_N = ["G", "E", "S"]
    if isnothing(Λ)
        labels_N = ["G"]
    end
    labels_ms = ["-1", "0", "+1"]  # only "0" for metastable state
    if isnothing(Ω₊) && isnothing(Ω₋)
        error("At least one of Ω₊ and Ω₋ must be given")
    else
        if isnothing(Ω₊)
            labels_ms = ["-1", "0"]
        elseif isnothing(Ω₋)
            labels_ms = ["0", "+1"]
        end
    end
    labels_mI = ["↑", "↓"]

    labels = [  # for Λ ≠
        ("G", ms_label, mI_label) for ms_label in labels_ms for mI_label in labels_mI
    ]
    if !isnothing(Λ)
        labels = [
            (N_label, ms_label, mI_label) for N_label in ["G", "E"] for
            ms_label in labels_ms for mI_label in labels_mI
        ]
        labels = [labels..., ("S", "0", "↑"), ("S", "0", "↓")]
    end

    NS = length(labels_ms)  # dimension of electronic spin subspace (G/E)
    N = length(labels)  # total dimension of Hilbert space

    function ket(N_label::String, ms_label::String, mI_label::String)
        @assert N_label ∈ labels_N
        @assert ms_label ∈ labels_ms
        @assert mI_label ∈ labels_mI
        index = findfirst(lt -> (lt == (N_label, ms_label, mI_label)), labels)
        Ψ = zeros(ComplexF64, N)
        Ψ[index] = one(ComplexF64)
        return Ψ
    end

    ket(ms_label, mI_label) = ket("G", ms_label, mI_label)

    function bra(args...)
        Ψ = ket(args...)
        return Array(Ψ')
    end

    function _truncate(H)
        @assert size(H) == (6, 6)
        if isnothing(Ω₊)  # keep Ω₋
            H = H[1:4, 1:4]
        elseif isnothing(Ω₋)  # keep Ω₊
            H = H[3:6, 3:6]
        end
        return H
    end

    H₀ = _truncate(
        Float64[
            (A_zz/2-B*γ_c/2+δ₋)       (A_zx/2)            0         0               0                   0
            (A_zx/2)            (-A_zz/2+B*γ_c/2+δ₋)      0         0               0                   0
            0                            0            (-B*γ_c/2)    0               0                   0
            0                            0                0      (B*γ_c/2)          0                   0
            0                            0                0         0       (-A_zz/2-B*γ_c/2+δ₊)      (-A_zx/2)
            0                            0                0         0            (-A_zx/2)       (A_zz/2+B*γ_c/2+δ₊)
        ]
    )

    # shift due to chirp
    H₁₋ = _truncate(
        Float64[
            -1 0  0  0  0  0
            0 -1  0  0  0  0
            0  0  0  0  0  0
            0  0  0  0  0  0
            0  0  0  0  0  0
            0  0  0  0  0  0
        ]
    )

    H₁₊ = _truncate(
        Float64[
            0  0  0  0  0  0
            0  0  0  0  0  0
            0  0  0  0  0  0
            0  0  0  0  0  0
            0  0  0  0 -1  0
            0  0  0  0  0 -1
        ]
    )

    # MW coupling
    H₂₋ = _truncate(
        Float64[
            0     0   0.5   0   0   0
            0     0    0   0.5  0   0
            0.5   0    0    0   0   0
            0    0.5   0    0   0   0
            0     0    0    0   0   0
            0     0    0    0   0   0
        ]
    )

    H₂₊ = _truncate(
        Float64[
            0   0   0   0   0   0
            0   0   0   0   0   0
            0   0   0   0  0.5  0
            0   0   0   0   0  0.5
            0   0  0.5  0   0   0
            0   0   0  0.5  0   0
        ]
    )

    c_ops = Any[]

    if isnothing(Λ) # ground-manifold only
        parts = Any[H₀,]
        if !isnothing(ω₋)
            push!(parts, [H₁₋, ω₋])
        end
        if !isnothing(ω₊)
            push!(parts, [H₁₊, ω₊])
        end
        if !isnothing(Ω₋)
            push!(parts, [H₂₋, Ω₋])
        end
        if !isnothing(Ω₊)
            push!(parts, [H₂₊, Ω₊])
        end
        H = hamiltonian(parts...)
        # TODO: use_dissipation of electronic spin levels
    else  # optical driving ⇒ include excited and metastable manifolds
        @assert (2 * NS) == size(H₀, 1)
        H₀ = H₀ ⊕ H₀ ⊕ ZeroMatrix(2)
        H₁₋ = H₁₋ ⊕ ZeroMatrix(2 * NS) ⊕ ZeroMatrix(2)
        H₁₊ = H₁₊ ⊕ ZeroMatrix(2 * NS) ⊕ ZeroMatrix(2)
        H₂₋ = H₂₋ ⊕ ZeroMatrix(2 * NS) ⊕ ZeroMatrix(2)
        H₂₊ = H₂₊ ⊕ ZeroMatrix(2 * NS) ⊕ ZeroMatrix(2)
        parts = Any[H₀,]
        if coherent_optical_drive
            H_Λ = ((0.5 * [0 1; 1 0]) ⊗ IdMatrix(2 * NS)) ⊕ ZeroMatrix(2)
            push!(parts, [H_Λ, Λ])
        else
            # For an incoherent drive, we'll manually add a time-dependent
            # excitation dissipator at the very end. This is because the
            # `liouvillian` function does not support time-dependent `c_ops`
        end
        if !isnothing(ω₋)
            push!(parts, [H₁₋, ω₋])
        end
        if !isnothing(ω₊)
            push!(parts, [H₁₊, ω₊])
        end
        if !isnothing(Ω₋)
            push!(parts, [H₂₋, Ω₋])
        end
        if !isnothing(Ω₊)
            push!(parts, [H₂₊, Ω₊])
        end
        H = hamiltonian(parts...)
        if Γ_E_G > 0.0
            for ms in labels_ms
                for mI in labels_mI
                    push!(c_ops, √Γ_E_G * ket("G", ms, mI) * bra("E", ms, mI))
                end
            end
        end
        if Γ_0_E_S > 0.0
            for mI in labels_mI
                push!(c_ops, √Γ_0_E_S * ket("S", "0", mI) * bra("E", "0", mI))
            end
        end
        if Γ_1_E_S > 0.0
            if !isnothing(Ω₊)
                for mI in labels_mI
                    push!(c_ops, √Γ_1_E_S * ket("S", "0", mI) * bra("E", "+1", mI))
                end
            end
            if !isnothing(Ω₋)
                for mI in labels_mI
                    push!(c_ops, √Γ_1_E_S * ket("S", "0", mI) * bra("E", "-1", mI))
                end
            end
        end
        if Γ_0_S_G > 0.0
            for mI in labels_mI
                push!(c_ops, √Γ_0_S_G * ket("G", "0", mI) * bra("S", "0", mI))
            end
        end
        if Γ_1_S_G > 0.0
            if !isnothing(Ω₊)
                for mI in labels_mI
                    push!(c_ops, √Γ_1_S_G * ket("G", "+1", mI) * bra("S", "0", mI))
                end
            end
            if !isnothing(Ω₋)
                for mI in labels_mI
                    push!(c_ops, √Γ_1_S_G * ket("G", "-1", mI) * bra("S", "0", mI))
                end
            end
        end
    end

    if use_dissipation
        L = liouvillian(H, c_ops; convention = :TDSE)
        if (!coherent_optical_drive) && (Λ ≠ 0)
            exc_ops = []
            for ms in labels_ms
                for mI in labels_mI
                    push!(exc_ops, ket("E", ms, mI) * bra("G", ms, mI))
                end
            end
            L_Λ = liouvillian(nothing, exc_ops; convention = :TDSE)
            L = Generator([L.ops..., L_Λ], [L.amplitudes..., Λ])
        end
        return L, ket, labels
    else
        return H, ket, labels
    end

end

end
