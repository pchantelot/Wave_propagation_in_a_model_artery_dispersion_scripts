# Assemble layers function
"""
    assemblelayer(Lkk,Lk,L0,Md,Lfkk,Lfk,Lf0,Mfd,bc,bcf,D0b,Bfkt,Bf0t)

Returns the matrices `Fkk`, `Fk`, `F0` and `Fd` which define the polynomial eigenproblem 
for a bilayer with PEP `Lkk`, `Lk`, `L0` and `Md` (`Lfkk`, `Lfk`, `Lf0` and `Mfd`, respectively).
The Matrices `D0b`, `Bfkt`, and `Bf0t` couple the normal displacement and stress at the interface.     
"""
function assemblelayer(Lkk,Lk,L0,Md,Lfkk,Lfk,Lf0,Mfd,bc,bcf,D0b,Bfkt,Bf0t)

    # Create extra blocks
    TRB = zeros(ComplexF64, size(Lkk)[1], size(Lfkk)[1])
    BLB = zeros(ComplexF64, size(Lfkk)[1], size(Lkk)[1])
    # Assemble kk term
    tt = hcat(Lkk,TRB)
    tb = hcat(BLB,Lfkk)
    Fkk = vcat(tt,tb)
    # Assemble k term
    TRBk = zeros(ComplexF64, size(Lkk)[1], size(Lfkk)[1])
    TRBk[bc[2],:] .= .- Bfkt
    tt = hcat(Lk,TRBk)
    tb = hcat(BLB,Lfk)
    Fk = vcat(tt,tb)
    # Assemble 0 term
    TRB0 = zeros(ComplexF64, size(Lkk)[1], size(Lfkk)[1])
    BLB0 = zeros(ComplexF64, size(Lfkk)[1], size(Lkk)[1])
    TRB0[bc[2],:] .= .- Bf0t
    BLB0[bcf[1],:] .= .- D0b
    tt = hcat(L0,TRB0)
    tb = hcat(BLB0,Lf0)
    F0 = vcat(tt,tb)
    # Assemble M term
    tt = hcat(Md,TRB)
    tb = hcat(BLB,Mfd)
    Fd = vcat(tt,tb)

    return Fkk, Fk, F0, Fd

end