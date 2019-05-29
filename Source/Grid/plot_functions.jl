
using Plots
pyplot()

# Plots energy histograms for all the energy series contained in E_by_T and saves the plot in energy_histogram.pdf in the
# current folder unless another filename is specified. It uses automatic determination of bins for each energy-series
# unless nbins is specified.
# uses StatsBase, Plots - pyplot
# -----------------------------------------------------------------------------------------------------------
function plotEnergyHistograms(E_by_T::Vector{Vector{T}}, temps::Vector{T}; 
        nbins=-1, filename="energy_histograms.pdf") where {T}
    
    N₀ = length(E_by_T)
    title = "Energy histograms for $(length(temps)) temperatures"
    # If we are using automatic number of bins
    if nbins==-1
        # Make first histogram
        hist = fit(Histogram, E_by_T[1])
        edges = collect(hist.edges[1])
        xs = [(edges[i] + edges[i+1])/2 for i = 1:length(edges)-1]
        ys = hist.weights
        plt = plot(xs, ys; label="T = $(round(temps[1]; digits=2))", xaxis="Energy pr. site", yaxis="number in bin", title=title);
        for k = 2:N₀
            hist = fit(Histogram, E_by_T[k])
            edges = collect(hist.edges[1])
            xs = [(edges[i] + edges[i+1])/2 for i = 1:length(edges)-1]
            ys = hist.weights
            plot!(plt, xs, ys; label="T = $(round(temps[k]; digits=2))");
        end
    else
        # If we are using fixed number of bins
        # Make first histogram
        hist = fit(Histogram, E_by_T[1]; nbins=nbins)
        edges = collect(hist.edges[1])
        xs = [(edges[i] + edges[i+1])/2 for i = 1:length(edges)-1]
        ys = hist.weights
        plt = plot(xs, ys; label="T = $(round(temps[1]; digits=2))", xaxis="Energy pr. site", yaxis="number in bin", title=title)
        for k = 2:N₀
            hist = fit(Histogram, E_by_T[k]; nbins=nbins)
            edges = collect(hist.edges[1])
            xs = [(edges[i] + edges[i+1])/2 for i = 1:length(edges)-1]
            ys = hist.weights
            plot!(plt, xs, ys; label="T = $(round(temps[k]; digits=2))");
        end
    end
    savefig(plt, filename);
    nothing
end


