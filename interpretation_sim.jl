using Distributions
using Plots
using RCall
using Random
using QuadGK
using LaTeXStrings
using ColorSchemes
using PGFPlotsX
using Test

pgfplotsx()

begin
    pgfplotsx()
    empty!(PGFPlotsX.CUSTOM_PREAMBLE)
end

theme(
    :default;
    background_color_legend = :transparent,
    foreground_color_legend = :transparent,
    grid = nothing,
    frame = :box,
    legend = :topright,
    legendfonthalign = :left,
    thickness_scaling = 2,
    size = (600, 450),
)



#------------------------------
# Figure 1
#------------------------------
σ = 0.5
π1_fig1 = 0.4
signal1_fig1 = MixtureModel([Cauchy(0, σ); Normal(0, 0)], [π1_fig1; 1 - π1_fig1])
y_grid = -5:0.01:5
alternative_pdfs1 =
    [quadgk(μ -> pdf(Normal(μ), y) * pdf(Cauchy(0, σ), μ), -Inf, +Inf)[1] for y in y_grid]
pdfs1_fig1 = (1 - π1_fig1) * pdf.(Normal(), y_grid) + π1_fig1 .* alternative_pdfs1
@test y_grid[501] == 0.0
f0_1_fig1 = pdfs1_fig1[501]
falt_0_1 = alternative_pdfs1[501]

ribbon_two_groups1 = zero(y_grid), alternative_pdfs1 .* π1_fig1

@test alternative_pdfs1 .* π1_fig1 ≈ pdfs1_fig1 .- (1 - π1_fig1) .* pdf.(Normal(), y_grid)

plot(
    y_grid,
    pdfs1_fig1,
    legend = :outertop,
    xlim = (-5, 5),
    ylim = (0, 0.37),
    xlabel = L"y",
    color = :black,
    ylabel = "Density",
    label = L"m_{\pi_1}(y) = (1-\pi_1) \phi(y) + \pi_1 m_1(y)",
)
plot!(
    y_grid,
    zero(y_grid),
    ribbon = ribbon_two_groups1,
    fillalpha = 0.4,
    fillcolor = "#018AC4",
    linealpha = 0,
    label = L"\pi_1 m_1(y)",
)

savefig("null_nonnull_two_groups.pdf")

ρ1_fig1 = π1_fig1 * (1 - falt_0_1 / pdf(Normal(), 0))

@test ρ1_fig1 ≈ 1 - f0_1_fig1/ pdf(Normal(),0)

ribbon_activity1 = zero(y_grid), (pdfs1_fig1 .- (1 - ρ1_fig1) .* pdf.(Normal(), y_grid))

plot(
    y_grid,
    pdfs1_fig1,
    legend = :outertop,
    xlim = (-5, 5),
    ylim = (0, 0.37),
    xlabel = L"y",
    color = :black,
    ylabel = "Density",
    label = L"\psi_{\rho}(y) = (1-\rho) \phi(y) + \rho \psi_1(y)",
)
plot!(
    y_grid,
    zero(y_grid),
    ribbon = ribbon_activity1,
    fillalpha = 0.7,
    fillcolor = :darkorange,
    linealpha = 0,
    label = L"\rho \psi_1(y)",
)
plot!(y_grid, (1 - ρ1) .* pdf.(Normal(), y_grid), label = L"\phi(y)")

savefig("inactive_active_two_groups.pdf")


#------------------------------
# Figure 3
#------------------------------

π1 = 0.2
signal1 = MixtureModel([Cauchy(0, σ); Normal(0, 0)], [π1; 1 - π1])
signal2 = Cauchy(0, σ / 5)


Gs_signal1 = cdf.(signal1, y_grid)
Gs_signal2 = cdf.(signal2, y_grid)

plot(
    y_grid,
    [Gs_signal1 Gs_signal2],
    label = [L"P_1 \; \textrm{(Dirac-Cauchy)}" L"P_2 \; \textrm{(Cauchy)}"],
    xlabel = L"\phantom{y}x\phantom{y}",
    ylabel = L"P(x)",
    xlim = (-5, 5),
    legend = :topleft,
    color = [ColorSchemes.seaborn_colorblind[10] ColorSchemes.seaborn_colorblind[8]],
    linewidth = 1.6,
    alpha = 1,
    linestyle = [:solid :dot],
    size = (700, 350),
)

savefig("illustration_sim_signal_cdfs.pdf")


pdfs1 = (1 - π1) * pdf.(Normal(), y_grid) + π1 .* alternative_pdfs1


pdfs2 = [
    quadgk(μ -> pdf(Normal(μ), y) * pdf(signal2, μ), -Inf, +Inf)[1] for y in y_grid
]

f0_1 = pdfs1[501]
f0_2 = pdfs2[501]

lfdr1 = (1-π1) .* pdf.(Normal(), y_grid) ./ pdfs1
lfdr2 = fill(0.0, length(y_grid))

sum(lfdr1 .* pdfs1) / sum(pdfs1)

liar1 = (f0_1 / pdf(Normal(), 0)) * pdf.(Normal(), y_grid) ./ pdfs1
liar2 = (f0_2 / pdf(Normal(), 0)) * pdf.(Normal(), y_grid) ./ pdfs2



plot(y_grid, [pdfs1 pdfs2])


Random.seed!(1)

n = 10_000
εs = randn(n)
μs2 = rand(signal2, n)
μs1 = rand(signal1, n)

Zs_2 = μs2 .+ εs
Zs_1 = μs1 .+ εs

Zs_2 = sort(Zs_2)
Zs_1 = sort(Zs_1)


@rput Zs_1
@rput Zs_2


R"""
library(locfdr)
locfdr_1 <- locfdr(Zs_1, nulltype =0, bre=seq(-7,7,by=0.05), plot=0)$fdr
locfdr_2 <- locfdr(Zs_2, nulltype =0,  bre=seq(-7,7,by=0.05), plot=0)$fdr
"""

@rget locfdr_1
@rget locfdr_2

R"""
library(fdrtool)
ps_1 <- 2*(1 - pnorm(abs(Zs_1)))
ps_2 <- 2*(1 - pnorm(abs(Zs_2)))
fdrtool_lfdr1 <- fdrtool(ps_1, statistic="pvalue")$lfdr
fdrtool_lfdr2 <- fdrtool(ps_2, statistic="pvalue")$lfdr
"""


@rget fdrtool_lfdr1
@rget fdrtool_lfdr2

R"""
library(qvalue)
qvalue_lfdr_1 <- qvalue(ps_1)$lfdr 
qvalue_lfdr_2 <- qvalue(ps_2)$lfdr 
"""

@rget qvalue_lfdr_1
@rget qvalue_lfdr_2


idx1_in_05 = findall(0 .<= Zs_1 .<= 5)
idx2_in_05 = findall(0 .<= Zs_2 .<= 5)

cols = ColorSchemes.seaborn_colorblind6#ColorSchemes.Set2_6

plot(
    y_grid,
    [lfdr1 liar1],
    label = ["lnsr (population)" "clar (population)"],
    color = [cols[1] cols[2]],
    xlabel = L"y",
    ylabel = L"\widehat{\mathrm{lfdr}}(y)",
    xlim = (0, 5),
    linewidth = 2.3,
)


plot!(
    Zs_1[idx1_in_05],
    [locfdr_1[idx1_in_05] fdrtool_lfdr1[idx1_in_05] qvalue_lfdr_1[idx1_in_05]],
    label = ["locfdr (estimate)" "fdrtool (estimate)" "qvalue (estimate)"],
    color = [cols[3] cols[4] cols[6]],
    linestyle = :solid,
    linewidth = 1.4,
)

savefig("illustration_sim_1_lfdr.pdf")

plot(
    y_grid,
    [lfdr2 liar2],
    label = ["lnsr (population)" "clar (population)"],
    color = [cols[1] cols[2]],
    xlabel = L"y",
    ylabel = L"\widehat{\mathrm{lfdr}}(y)",
    xlim = (0, 5),
    linewidth = 2.3,
)

plot!(
    Zs_2[idx2_in_05],
    [locfdr_2[idx2_in_05] fdrtool_lfdr2[idx2_in_05] qvalue_lfdr_2[idx2_in_05]],
    label = ["locfdr (estimate)" "fdrtool (estimate)" "qvalue (estimate)"],
    color = [cols[3] cols[4] cols[6]],
    linestyle = :solid,
    linewidth = 1.4,
)

savefig("illustration_sim_2_lfdr.pdf")


R"""
library(sessioninfo)
session_info(to_file = "R_session_info.txt")
"""