# Validate approximations to compute guided waves in curved clamped elastic strips
# POnG
path_to_repository = "Documents/GitHub/Wave_propagation_in_a_model_artery_dispersion_scripts/"
push!(LOAD_PATH,joinpath(pwd(),path_to_repository))
using POnG
using CairoMakie
using FileIO, UnPack, MAT, NaturalSort

# Load Data
reffilenames = filter(x -> occursin(".mat", x), readdir(joinpath(@__DIR__,"Comsol")))
filenames = filter(x -> occursin(".jld2", x), readdir(@__DIR__))
filenames = sort(filenames, lt = natural)

with_theme(theme_latexfonts()) do
    fig = Figure(size = (2 * 1.2 * 360,360))

        ax = Axis(fig[1,1])
        ax.limits = (0, 900, 50, 400)
        ax.xlabel = L"k \,\, \mathrm{(1/m)}"
        ax.ylabel = L"f\,\, \mathrm{(Hz)}"
        ax.xlabelsize = 20 
        ax.ylabelsize = 20
        ax.title = "Out-of-plane first mode"

        for name in filenames
            Data = load(joinpath(@__DIR__,name))
            @unpack data = Data
            k = data.k
            f = data.ω
            filter!(!isnan, k)
            filter!(!isnan, f)
            A = unique(f)
            knew = []
            for i in eachindex(A)
                idxtemp = findall(x-> x == A[i],f)
                push!(knew,maximum(real(k[idxtemp])))
            end
            scatter!(ax, knew, A)
        end

        for name in reffilenames
            Data = matread(joinpath(@__DIR__,"Comsol",name))
            @unpack k, eigenfreq, strip = Data
            k = dropdims(real.(k); dims = 1)
            lines!(ax, k, imag(eigenfreq[:,1])./2 ./ pi, linestyle = :solid,
            label = string("1/κ = ", strip["courbure"]," m"))
        end

        fig[1,2] = Legend(fig,ax,"Comsol")
        colsize!(fig.layout, 1, Aspect(1, 1.2))

        resize_to_layout!(fig)
    display(fig)
end
# NB: for the lowest radius of curvature, the low-frequency difference comes from a crossing between the first and second modes.  