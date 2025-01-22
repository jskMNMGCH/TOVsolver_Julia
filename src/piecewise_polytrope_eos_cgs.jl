"""
This code is written in CGS unit system!
In accordance with convention, the argument of 'log_p1', where p1 has units of dyn/cm^2.
"""
module PiecePoly

"""
Following construction of EoS is based on the piecewise polytrope framework.
Ref. ↓
Constraints on a phenomenologically parameterized neutron-star equation of state
Jocelyn S. Read (Wisconsin U., Milwaukee), Benjamin D. Lackey (Wisconsin U., Milwaukee), 
Benjamin J. Owen (Granada U., Theor. Phys. Astrophys. and Penn State U.), John L. Friedman (Wisconsin U., Milwaukee)
e-Print: 0812.2163 [astro-ph]
DOI: 10.1103/PhysRevD.79.124032
Published in: Phys.Rev.D 79 (2009), 124032
"""

# include files
include("constants.jl")

# Parameters of "Crust EoS" are fixed.
# ("crust name", upper bound of mass density 'rho' [g/cm^3], Polytrope index 'Gamma', Polytrope Coefficient 'K')
fixed_crust = [
    ("crust1", 2.44034e7, 1.58425, 6.80110e-9),
    ("crust2", 3.78358e11, 1.28733, 1.06186e-6),
    ("crust3", 2.62780e12, 0.62223, 5.32697e1),
    ("crust4", nothing, 1.35692, 3.99874e-8)
]

fixed_ρ_bounds = [10^14.7, 10^15]

# export fixed_crust, fixed_ρ_bounds

function p_rho_def(rho, K, Gamma)
"""
    Computes the pressure as a function of density, based on the polytropic relation.
    Inputs:
    - `rho`: Density.
    - `K`: Polytropic constant.
    - `Gamma`: Polytropic index.
    Returns:
    - Pressure corresponding to the given density.
"""
    return K * (rho^Gamma)
end

function eps_rho_def(rho, K, Gamma, a)
"""
    Calculates the energy density as a function of density, using a specific polytropic model.
    Inputs:
    - `rho`: Density.
    - `K`: Polytropic constant.
    - `Gamma`: Polytropic index.
    - `a`: Constant related to energy contributions.
    Returns:
    - Energy density corresponding to the given density.
"""
    return (1.0 + a) * rho + K * (rho^Gamma) / (Gamma - 1.0)
end

function a_def(pre_eps_lim, pre_rho_lim, K, Gamma; param_crust=fixed_crust)
"""
    Determines the parameter `a`, which adjusts the energy density at given density and pressure limits.
    Inputs:
    - `pre_eps_lim`: Energy density at the lower limit.
    - `pre_rho_lim`: Density at the lower limit.
    - `K`: Polytropic constant.
    - `Gamma`: Polytropic index.
    Returns:
    - Parameter `a`.
"""
    if pre_eps_lim <= eps_rho_def(param_crust[1][2], param_crust[1][4], param_crust[1][3], 0.0)
        return 0.0
    else
        return pre_eps_lim / pre_rho_lim - 1.0 - K * (pre_rho_lim^(Gamma - 1.0)) / (Gamma - 1.0)
    end
end

function next_K_def(p_lim, rho_lim, next_Gamma)
"""
    Calculates the polytropic constant `K` for the next density segment.
    Inputs:
    - `p_lim`: Pressure at the upper boundary of the current segment.
    - `rho_lim`: Density at the upper boundary of the current segment.
    - `next_Gamma`: Polytropic index for the next segment.
    Returns:
    - Polytropic constant `K` for the next segment.
"""
    return p_lim / (rho_lim^next_Gamma)
end

function calc_K1_rhob(log_p1, Gamma1; K_crust=fixed_crust[4][4], Gamma_crust=fixed_crust[4][3], param_bound=fixed_ρ_bounds)
"""
    Computes the initial polytropic constant `K1` and baryon density `rho_b` for the first inner segment.
    Inputs:
    - `log_p1`: Logarithm of the pressure at the first deviding density. [dyn/cm^2] = [c^2 g/cm^3]
    - `Gamma1`: Polytropic index for the first inner segment.
    - `K_crust`: Polytropic constant for the crust (default value provided).
    - `Gamma_crust`: Polytropic index for the crust (default value provided).
    Returns:
    - `K1`: Initial polytropic constant.
    - `rho_b`: Baryon density at the transition.
"""
    K1 = 10^log_p1 / c^2 / param_bound[1]^Gamma1  # [dyn/cm^2] = [c^2 g/cm^3]
    rho_b = (K_crust / K1)^(1.0 / (Gamma1 - Gamma_crust))
    return K1, rho_b
end

function calc_outerEoS(;param_crust=fixed_crust)
"""
    Calculates the list of `a` parameters for the outer crust equations of state.
    Inputs:
    - `param_crust`: Parameters for the crust segments (default to `fixed_crust`).
    Returns:
    - List of `a` parameters for the outer crust.
"""
    a_crust_list = [0.0]
    eps_pre_crust = eps_rho_def(param_crust[1][2], param_crust[1][4], param_crust[1][3], 0.0)

    for i in 2:length(param_crust)
        rho_temp = param_crust[i][2]
        a_temp = a_def(eps_pre_crust, param_crust[i - 1][2], param_crust[i][4], param_crust[i][3])
        push!(a_crust_list, a_temp)
        if i < length(param_crust)
            eps_pre_crust = eps_rho_def(rho_temp, param_crust[i][4], param_crust[i][3], a_crust_list[i])
        end
    end
    return a_crust_list
end

function param_of_innerEoS(log_p1, Gamma; param_c=fixed_crust, param_bound=fixed_ρ_bounds)
"""
    Defines the parameters for the inner equations of state based on polytropic segments.
    Inputs:
    - `log_p1`: Logarithm of the transition pressure.
    - `Gamma`: List of polytropic indices for the inner segments.
    - `param_c`: Parameters for the crust (default to `fixed_crust`).
    Returns:
    - Lists of density limits, `a` parameters, and polytropic constants for the inner EoS.
"""
    K1, rhob = calc_K1_rhob(log_p1, Gamma[1])
    rho_lim_list = [rhob, param_bound[1], param_bound[2], Inf]
    K_list = [K1]
    a_list = []

    K_crust_end = param_c[end][4]
    Gamma_crust_end = param_c[end][3]
    eps_temp = eps_rho_def(rho_lim_list[1], K_crust_end, Gamma_crust_end, calc_outerEoS(param_crust=param_c)[end])

    for i in 1:length(Gamma)
        a_temp = a_def(eps_temp, rho_lim_list[i], K_list[i], Gamma[i])
        push!(a_list, a_temp)
        # println(rho_lim_list[i+1])
        p_temp = p_rho_def(rho_lim_list[i+1], K_list[i], Gamma[i])
        eps_temp = eps_rho_def(rho_lim_list[i+1], K_list[i], Gamma[i], a_list[i])

        if i < length(Gamma)
            K_temp = next_K_def(p_temp, rho_lim_list[i+1], Gamma[i + 1])
            push!(K_list, K_temp)
        end
    end
    return rho_lim_list, a_list, K_list
end


function joint_params(rho_lim_l, a_l, K_l; par_crust=fixed_crust)
"""
    Combines the crust and inner equations of state parameters into a single set of lists.
    Inputs:
    - `rho_lim_l`: Density limits from the inner EoS.
    - `a_l`: `a` parameters from the inner EoS.
    - `K_l`: Polytropic constants from the inner EoS.
    - `par_crust`: Crust parameters (default to `fixed_crust`).
    Returns:
    - Combined density limits, `a` parameters, and polytropic constants.
"""
    rho_lim_all = vcat([c[2] for c in par_crust[1:end-1]], rho_lim_l)
    a_all = vcat(calc_outerEoS(param_crust=par_crust), a_l)
    K_all = vcat([c[4] for c in par_crust], K_l)
    return rho_lim_all, a_all, K_all
end


function get_all_params(log_p1, Gamma; p_c=fixed_crust)
"""
    Computes all parameters required to describe the equations of state.
    Inputs:
    - `log_p1`: Logarithm of the transition pressure.
    - `Gamma`: List of polytropic indices.
    - `p_c`: Parameters for the crust (default to `fixed_crust`).
    Returns:
    - Density limits, `a` parameters, polytropic constants, and polytropic indices for the full EoS.
"""
    rho, a, k = param_of_innerEoS(log_p1, Gamma)
    rho_lim_all, a_all, K_all = joint_params(rho, a, k)
    Gamma_all = vcat([c[3] for c in p_c], Gamma)
    return rho_lim_all, a_all, K_all, Gamma_all
end

function make_polyEos(rho_lim_arr, a_arr, K_arr, Gamma_arr; initial_rho=10, final_rho=1e18)
"""
    Function to create a piecewise polytropic equation of state (EoS).
    
    Arguments:
    - `rho_lim_arr`: Array of mass density defining the upper boundaries of each polytropic segment. [g/cm^3]
    - `a_arr`: Array of coefficients for energy density calculations.
    - `K_arr`: Array of polytropic constants for each segment.
    - `Gamma_arr`: Array of polytropic indices (Γ) for each segment.
    - `initial_rho`: Starting density in CGS units (default: 10). [g/cm^3]
    - `final_rho: Final density in CGS units (default: 1e18). [g/cm^3]
    
    # Returns:
    - A tuple `(ε, p)`:
      - `ε`: Array containing the energy density values for the entire density range. [g/cm^3]
      - `p`: Array containing the pressure values for the entire density range. [g/cm^3]
    
    # Notes:
    This function constructs energy density (`ε`) and pressure (`p`) values for a range of densities
    using piecewise polytropic equations. Each density segment is processed individually based on
    the provided parameters.
"""
    ε = []
    p = []
    if rho_lim_arr[end] == Inf
        vs = [[initial_rho], rho_lim_arr[1:end-1], [final_rho]]  # 10^19 = -8446744073709551616 overflow!
        rho_boundary = reduce(vcat, vs)
    else
        vs = [[initial_rho], rho_lim_arr]
        rho_boundary = reduce(vcat, vs)
    end
    for i in 1:length(rho_boundary)-1
        logrange_len = ceil(Int, (rho_boundary[i+1]-rho_boundary[i])/rho_boundary[i+1]*2^9)  # 少数点以下を切り上げて，　Intに変換
        rho_piece = logrange(rho_boundary[i], rho_boundary[i+1], length=logrange_len)
        ε = vcat(ε, PiecePoly.eps_rho_def.(rho_piece, K_arr[i], Gamma_arr[i], a_arr[i]))
        p = vcat(p, PiecePoly.p_rho_def.(rho_piece, K_arr[i], Gamma_arr[i]))
    end
    return ε, p
end

end # end of PiecePoly