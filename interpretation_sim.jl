using Distributions
using Plots
using RCall
using Random
using QuadGK
using LaTeXStrings
using ColorSchemes
using PGFPlotsX
using Test
using Roots

pgfplotsx()
empty!(PGFPlotsX.CUSTOM_PREAMBLE)

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
# Figure 2
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
    label = L"\psi_{\rho_1}(y) = (1-\rho_1) \phi(y) + \rho_1 \psi_1(y)",
)
plot!(
    y_grid,
    zero(y_grid),
    ribbon = ribbon_activity1,
    fillalpha = 0.7,
    fillcolor = :darkorange,
    linealpha = 0,
    label = L"\rho_1 \psi_1(y)",
)

savefig("inactive_active_two_groups.pdf")


#------------------------------
# Figure 6
#------------------------------
# Activity decomposition for asymmetric distribution

# 
μ₂ = 2.0
π₂ = 0.6
postmean_numerator(x) = x*pdf(Normal(), x)*(1-π₂) + μ₂*pdf(Normal(),μ₂)*π₂


μ₁ = fzero( x-> x*pdf(Normal(), x)*(1-π₂) + μ₂*pdf(Normal(),μ₂)*π₂, -1, 0)

pdfs_asymm =  0.8*pdf.(Normal(), y_grid) .+ 0.2*((1 - π₂) * pdf.(Normal(μ₁), y_grid) + π₂ * pdf.(Normal(μ₂), y_grid))

ρ1_asymm = 1-pdfs_asymm[501] / pdf(Normal(), 0)

tmp_activity = pdfs_asymm .- (1 - ρ1_asymm) .* pdf.(Normal(), y_grid)
ribbon_activity_asymm = zero(y_grid), tmp_activity


plot(
    y_grid,
    pdfs_asymm,
    legend = :outertop,
    xlim = (-5, 5),
    ylim = (0, 0.37),
    xlabel = L"y",
    color = :black,
    ylabel = "Density",
    label = L"\psi_{\rho_1}(y) = (1-\rho_1) \phi(y) + \rho_1 \psi_1(y)",
)
plot!(
    y_grid,
    zero(y_grid),
    ribbon = ribbon_activity_asymm,
    fillalpha = 0.7,
    fillcolor = :darkorange,
    linealpha = 0,
    label = L"\rho_1 \psi_1(y)",
)

savefig("asymmetric_active_two_groups.pdf")

# temp

loc_fdr = (1-ρ1_asymm) .* pdf.(Normal(), y_grid) ./ pdfs_asymm
# q: is sech identity still true?
using Empirikos
prior = DiscreteNonParametric([0; 2.0; μ₁], [0.8; 0.12; 0.08])
function exp_formula(y)
    posterior = Empirikos.posterior(StandardNormalSample(y), prior)
    #sum( sech.(y .* support(posterior)) .* probs(posterior))
    sum( exp.(-y .* support(posterior)) .* probs(posterior))
end 

function sech_formula(y)
    posterior = Empirikos.posterior(StandardNormalSample(y), prior)
    #sum( sech.(y .* support(posterior)) .* probs(posterior))
    #myf(x) = 1/(exp(x)-x)
    myf(x) = sech(x)
    sum( myf.(y .* support(posterior)) .* probs(posterior))
end 


function sign_rate_formula(y)
    posterior = Empirikos.posterior(StandardNormalSample(y), prior)
    #sum( sech.(y .* support(posterior)) .* probs(posterior))
    y <= 0 ? sum(probs(posterior)[2:3]) : sum(probs(posterior)[1:2])
end

plot(y_grid, loc_fdr)
plot!(y_grid, sign_rate_formula.(y_grid))

plot(y_grid, loc_fdr)
plot!(y_grid, sech_formula.(y_grid))
plot!(y_grid, sech_formula.(y_grid))

plot(y_grid, sech_formula.(y_grid) .- loc_fdr)

f(y) =  sum( [0.8; 0.12; 0.08] .* [(tanh(x*y)+1)*exp(-x^2/2) for x in [0; 2.0; μ₁]])
f(1.0)
f(-5.0)


# plot(y_grid, (1-ρ1_asymm).* pdf.(Normal(), y_grid) ./ pdfs_asymm)
# plot(y_grid, log.(pdfs_asymm ./ pdf.(Normal(), y_grid)))
# plot(y_grid, pdfs_asymm ./ pdf.(Normal(), y_grid))




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
    label = [L"\mathrm{P}_1 \; \textrm{(Dirac-Cauchy)}" L"\mathrm{P}_2 \; \textrm{(Cauchy)}"],
    xlabel = L"\phantom{y}x\phantom{y}",
    ylabel = L"\mathrm{P}(x)",
    xlim = (-5, 5),
    legend = :topleft,
    color = [ColorSchemes.seaborn_colorblind[10] ColorSchemes.seaborn_colorblind[8]],
    linewidth = 1.6,
    alpha = 1,
    linestyle = [:solid :dot],
    size = (600,350) #(700, 350),
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

# Locfdr

R"""
library(locfdr)
locfdr_1 <- locfdr(Zs_1, nulltype =0, bre=seq(-7,7,by=0.05), plot=0)$fdr
locfdr_2 <- locfdr(Zs_2, nulltype =0,  bre=seq(-7,7,by=0.05), plot=0)$fdr
"""

@rget locfdr_1
@rget locfdr_2

# Fdrtool

R"""
library(fdrtool)
ps_1 <- 2*(1 - pnorm(abs(Zs_1)))
ps_2 <- 2*(1 - pnorm(abs(Zs_2)))
fdrtool_lfdr1 <- fdrtool(ps_1, statistic="pvalue")$lfdr
fdrtool_lfdr2 <- fdrtool(ps_2, statistic="pvalue")$lfdr
"""


@rget fdrtool_lfdr1
@rget fdrtool_lfdr2

# Qvalue

R"""
library(qvalue)
qvalue_lfdr_1 <- qvalue(ps_1)$lfdr 
qvalue_lfdr_2 <- qvalue(ps_2)$lfdr 
"""

@rget qvalue_lfdr_1
@rget qvalue_lfdr_2

# TailInflation

R"""
source('TailInflation/TailInflationA.R')
"""

R"""
tail_inflation1 <- TailInflation.A(Zs_1[Zs_1 <= 7 & Zs_1 >= -7])
"""
R"""
tail_inflation2 <- TailInflation.A(Zs_2[Zs_2 <= 7 & Zs_1 >= -7])
"""

R"""
D_tt_1 <- LocalLinearSplines1.2A(tail_inflation1$tau,Zs_1)
thh_tt_1 <- D_tt_1 %*% tail_inflation1$theta
thh_tt_at_0_1 <-  (LocalLinearSplines1.2A(tail_inflation1$tau,0) %*% tail_inflation1$theta)[1]
tailinflation_lfdr_1 <- exp(thh_tt_at_0_1)/exp(thh_tt_1)
"""

R"""
D_tt_2 <- LocalLinearSplines1.2A(tail_inflation2$tau,Zs_2)
thh_tt_2 <- D_tt_2 %*% tail_inflation2$theta
thh_tt_at_0_2 <-  (LocalLinearSplines1.2A(tail_inflation2$tau,0) %*% tail_inflation2$theta)[1]
tailinflation_lfdr_2 <- exp(thh_tt_at_0_2)/exp(thh_tt_2)
"""

@rget tailinflation_lfdr_1
@rget tailinflation_lfdr_2



idx1_in_05 = findall(0 .<= Zs_1 .<= 5)
idx2_in_05 = findall(0 .<= Zs_2 .<= 5)


colors = [
    RGB(0.00, 0.45, 0.70),  # Blue
    RGB(0.90, 0.20, 0.00),  # Bright Red
    RGB(0.00, 0.65, 0.45),  # Bright Teal
    RGB(0.55, 0.25, 0.15),  # Rich Brown
    RGB(0.50, 0.00, 0.50),  # Deep Purple
    RGB(0.70, 0.50, 0.00),  # Golden Brown
]

plot(
    y_grid,
    [lfdr1 liar1],
    label = ["lnsr (population)" "clar (population)"],
    color = [colors[1] colors[2]],
    xlabel = L"y",
    ylabel = L"\widehat{\mathrm{lfdr}}(y)",
    xlim = (0, 5),
    linewidth = 3.3,
    size = (650, 450),
    alpha = [0.5 0.5]
)

plot!(
    Zs_1[idx1_in_05],
    [locfdr_1[idx1_in_05] fdrtool_lfdr1[idx1_in_05] qvalue_lfdr_1[idx1_in_05] tailinflation_lfdr_1[idx1_in_05]],
    label = ["locfdr (estimate)" "fdrtool (estimate)" "qvalue (estimate)" "tailinflation (estimate)"],
    color = [colors[6] colors[4] colors[3] colors[5]],
    linewidth = 1.2,
    alpha = 0.8,
    linestyle = [:solid :dash :solid :dot]  
)

savefig("illustration_sim_1_lfdr.pdf")



plot(
    y_grid,
    [lfdr2 liar2],
    label = ["lnsr (population)" "clar (population)"],
    color = [colors[1] colors[2]],
    xlabel = L"y",
    ylabel = L"\widehat{\mathrm{lfdr}}(y)",
    xlim = (0, 5),
    linewidth = 3.3,
    size = (650, 450),
    alpha = [0.5 0.5]
)

plot!(
    Zs_2[idx2_in_05],
    [locfdr_2[idx2_in_05] fdrtool_lfdr2[idx2_in_05] qvalue_lfdr_2[idx2_in_05] tailinflation_lfdr_2[idx2_in_05]],
    label = ["locfdr (estimate)" "fdrtool (estimate)" "qvalue (estimate)" "tailinflation (estimate)"],
    color = [colors[6] colors[4] colors[3] colors[5]],
    linewidth = 1.2,
    alpha = 0.8,
    linestyle = [:solid :dash :solid :dot] 
)


savefig("illustration_sim_2_lfdr.pdf")


R"""
library(sessioninfo)
session_info(to_file = "R_session_info.txt")
"""