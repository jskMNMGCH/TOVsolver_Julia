"""
This code is written in Geometrized unit system, where distances are expressed in centimeters [cm].
If energy_density and pressure are given in [g/cm^3] unit, 
"""
module MainModule

# include files
using Base.Threads
using ProgressMeter
include("solver_code.jl") # includeにより solver_code.jl の内容を現在のスコープに取り込む.
include("constants.jl") # includeにより solver_code.jl の内容を現在のスコープに取り込む.
using .SolverCode

# functions
function make_eos_monotonic(e, P)
"""
    This function processes the given arrays `e` (energy density) and `P` (pressure) to remove points where the energy density and the pressure do not increase monotonically.
    
    Args:
        e (array): Array of energy densities.
        P (array): Array of pressures.
    
    Returns:
        mono_e (array): Monotonic energy density values.
        mono_P (array): Monotonic pressure values.
"""
    mono_e = [e[1]]
    mono_P = [P[1]]
    for i in 2:length(e)
        dP = P[i] - P[i-1]
        de = e[i] - e[i-1]
        # The region where dp = 0 corresponds to the area during a first-order phase transition, so it should be removed.
        # 1次相転移領域が十分にソフトな領域として与えられるようにEoSの離散点を与えないといけないことに注意．
        # 対数グリッドでそのまま離散点を与えると，高密度にある1次相転移は誤差が大きくなりがちなので注意
        if dP > 0.0 && de > 0.0
            push!(mono_e, e[i]) 
            push!(mono_P, P[i])
        end
    end
    return mono_e, mono_P
end


# function out_RMT(ε, pres; 
#     debug::Bool=false, 
#     parallel::Bool=false,
#     min_pc=1.0*MeVfm3_to_gcm3/unit_g, 
#     max_pc=1.3e3*MeVfm3_to_gcm3/unit_g, 
#     num_pc::Int=100, 
#     dr::Float64=1e2
#     )
# """ 
#     This function calculates the radius, mass, and tidal deformability by solving the Tolman–Oppenheimer–Volkoff (TOV) equations for a given set of energy densities (`ε`) and pressures (`pres`).
    
#     Args:
#         ε (array): Array of energy densities ([1/cm^2], [g/cm^3] is scaled by eps_ref [g/cm]).
#         pres (array): Array of pressures ([1/cm^2], [g/cm^3] is scaled by eps_ref [g/cm]).
#         debug (bool, optional): Flag to enable debugging mode. Default is `false`.
    
#     Returns:
#         Tuple: A tuple containing:
#             - `R` (1D array): Array of radii corresponding to the input densities.
#             - `M` (1D array): Array of masses corresponding to the input densities.
#             - `Λ` (1D array): Array of tidal deformabilities corresponding to the input densities.
#             - `sol_list` (1D array): List of solutions for each central density.
# """
#     # The dimensionless pressure 2.9e-6 corresponds to about 3.0 MeV/fm^3 in QHC18_gv100H16.
#     # The dimensionless pressure 3.8e-3 corresponds to about 1300 MeV/fm^3 in QHC18_gv100H16.
#     min_pc_id = searchsortedfirst(pres, min_pc)
#     max_pc_id = min(searchsortedfirst(pres, max_pc), length(pres)-1)
#     len = min(max_pc_id-min_pc_id, num_pc)
#     pc_indices = round.(Int, range(min_pc_id, max_pc_id, length=len))

#     R        = Vector{Float64}(undef, len)
#     M        = Vector{Float64}(undef, len)
#     Λ        = Vector{Float64}(undef, len)
#     k_mag    = Vector{Float64}(undef, len)
#     sol_list = Vector{Any}(undef, len)

#     if parallel
#         # progress = Progress(len, 1, "Total number of processes = $(len)")
#         @threads for j in 1:len
#             i = pc_indices[j]
#             debug_flag = debug && (i == pc_indices[1])
#             try
#                 result, sol = SolverCode.solveTOV_RMT(i, ε, pres, debug_flag, max_dr=dr)
#                 R[j]        = result[1]
#                 M[j]        = result[2]
#                 Λ[j]        = result[3]
#                 k_mag[j]    = result[4]
#                 sol_list[j] = sol
#             catch err
#                 println("Error in thread $(threadid())  (j=$j, i=$i) : ", err)
#             end
#             # next!(progress)
#         end
#     else
#         for j in 1:len
#             i = pc_indices[j]
#             debug_flag = debug && (i == pc_indices[1])
#             try
#                 result, sol = SolverCode.solveTOV_RMT(i, ε, pres, debug_flag, max_dr=dr)
#                 R[j]        = result[1]
#                 M[j]        = result[2]
#                 Λ[j]        = result[3]
#                 k_mag[j]    = result[4]
#                 sol_list[j] = sol
#             catch err
#                 println("Error (j=$j, i=$i) : ", err)
#             end
#         end
#     end
#     return [R, M, Λ, k_mag], sol_list
# end

function out_RMT(ε, pres;
    debug::Bool=false,
    parallel::Bool=false,
    min_pc = 1.0*MeVfm3_to_gcm3/unit_g,
    max_pc = 1.3e3*MeVfm3_to_gcm3/unit_g,
    num_pc::Int=100,
    dr::Float64=1e2,
    keep_sol::Bool=false          # ★sol を保持するかどうかだけ切替
)
"""
    out_RMT(ε, pres; keep_sol=false, ...)

- keep_sol=false: sol は保持せず、各反復で破棄する
- keep_sol=true : sol_list に保持して返す

Returns:
  keep_sol=false -> (RMT, nothing)
  keep_sol=true  -> (RMT, sol_list)

RMT = [R, M, Λ, k_mag]
"""
    min_pc_id = searchsortedfirst(pres, min_pc)
    max_pc_id = min(searchsortedfirst(pres, max_pc), length(pres)-1)
    len = min(max_pc_id - min_pc_id, num_pc)
    pc_indices = round.(Int, range(min_pc_id, max_pc_id, length=len))

    R     = Vector{Float64}(undef, len)
    M     = Vector{Float64}(undef, len)
    Λ     = Vector{Float64}(undef, len)
    k_mag = Vector{Float64}(undef, len)          # ★常に必要
    sol_list = keep_sol ? Vector{Any}(undef, len) : nothing

    function do_one!(j::Int)
        i = pc_indices[j]
        debug_flag = debug && (i == pc_indices[1])
        try
            if keep_sol
                result, sol = SolverCode.solveTOV_RMT(i, ε, pres, debug_flag, max_dr=dr)
                R[j]     = result[1]
                M[j]     = result[2]
                Λ[j]     = result[3]
                k_mag[j] = result[4]
                sol_list[j] = sol
            else
                # sol を保存しない（このスコープを抜けたら参照が消える）
                result, _ = SolverCode.solveTOV_RMT(i, ε, pres, debug_flag, max_dr=dr)
                R[j]     = result[1]
                M[j]     = result[2]
                Λ[j]     = result[3]
                k_mag[j] = result[4]
            end
        catch err
            if parallel
                println("Error in thread $(threadid()) (j=$j, i=$i): ", err)
            else
                println("Error (j=$j, i=$i): ", err)
            end
            R[j] = NaN; M[j] = NaN; Λ[j] = NaN; k_mag[j] = NaN
            if keep_sol
                sol_list[j] = nothing
            end
        end
        return nothing
    end

    if parallel
        Threads.@threads for j in 1:len
            do_one!(j)
        end
    else
        for j in 1:len
            do_one!(j)
        end
    end

    return [R, M, Λ, k_mag], sol_list
end

# function truncate_rows_based_on_second_row(array::Vector{Vector{Float64}})
#     second_row = array[2]
#     for i in 2:length(second_row)-1
#         if (second_row[i-1] < second_row[i] && second_row[i] > second_row[i+1]) 
#             if second_row[1] < second_row[end]
#                 return [row[1:i] for row in array]
#             elseif second_row[end] < second_row[1]
#                 return [row[i:end] for row in array]
#             end
#         end
#     end
#     println("\nNo local maxima found in the second row.")
#     return array
# end

function truncate_rows_based_on_second_row(array::Vector{Vector{Float64}};
        Mmin::Float64 = -Inf,
        Mmax::Float64 =  Inf)

    M = array[2]
    n = length(M)

    # 探索対象のインデックス（M が範囲内）
    idxs = findall(m -> isfinite(m) && (Mmin <= m <= Mmax), M)
    length(idxs) < 3 && return array  # 局所最大に必要な点数がない

    # 方向判定（先頭→末尾が増加なら左側を残す、減少なら右側を残す）
    increasing = M[first(idxs)] < M[last(idxs)]

    # idxs 内で局所最大を探す（隣接点も idxs 内にあるところだけ）
    for t in 2:length(idxs)-1
        i_prev = idxs[t-1]
        i      = idxs[t]
        i_next = idxs[t+1]

        if (M[i_prev] < M[i] && M[i] > M[i_next])
            cut = i
            return increasing ? [row[1:cut] for row in array] : [row[cut:end] for row in array]
        end
    end

    println("No local maxima found in M within [$(Mmin), $(Mmax)].")
    return array
end

function out_RMT_point(ε, pres, center_pres::Float64; debug::Bool=false)
    i = searchsortedfirst(pres, center_pres)
    result, sol = SolverCode.solveTOV_RMT(i, ε, pres, debug)
    R = result[1]
    M = result[2]  
    Λ = result[3]
    return [R, M, Λ], sol
end

function cs(energy_density, pressure)
"""
    This function calculates the speed of sound (c_s) for a given array of energy densities (`energy_density`) 
    and pressures (`pressure`).

    Args:
        energy_density (array): Array of energy density values (ε) at discrete points.
        pressure (array): Array of pressure values (P) corresponding to `energy_density`.

    Returns:
        Tuple: A tuple containing:
            - `energy_density[2:end]` (1D array): Truncated array of energy densities, starting from the second element.
            - `cs` (1D array): Array of speed of sound values calculated using the formula: cs^2 = dp/dε.
"""
    cs = [sqrt( (pressure[i]-pressure[i-1])/(energy_density[i]-energy_density[i-1]) ) for i in 2:length(pressure)]
    return energy_density[2:end], cs
end


function cs2(energy_density, pressure)
"""
    This function calculates the square of the speed of sound (c_s^2) for a given array of energy densities (`energy_density`) 
    and pressures (`pressure`).

    Args:
        energy_density (array): Array of energy density values (ε) at discrete points.
        pressure (array): Array of pressure values (P) corresponding to `energy_density`.

    Returns:
        Tuple: A tuple containing:
            - `energy_density[2:end]` (1D array): Truncated array of energy densities, starting from the second element.
            - `cs2` (1D array): Array of speed of sound values calculated using the formula: cs^2 = dp/dε.
"""
    cs2 = [(pressure[i]-pressure[i-1])/(energy_density[i]-energy_density[i-1]) for i in 2:length(pressure)]
    return energy_density[2:end], cs2
end


function check_tidal_constraint(RMT; M_interval::Float64=1e-2)
    """
    This function checks if at least one tidal deformability value (`T`) at `M = 1.4` satisfies the specified constraint.

    Args:
        RMT (tuple or array): A data structure containing:
            - `RMT[2]` (array): Array of mass values corresponding to compact stars.
            - `RMT[3]` (array): Array of tidal deformability values corresponding to the masses in `RMT[2]`.

    Returns:
        Bool: `true` if at least one value satisfies the constraint; otherwise, `false`.
    """
    # Sort masses and get sorting indices
    sorted_indices = sortperm(RMT[2])  # Indices to sort RMT[2]
    sorted_M = RMT[2][sorted_indices]
    sorted_T = RMT[3][sorted_indices]  # Sort T accordingly

    # Define the valid range for tidal deformability (90% Confidence Interval)
    max_threshold = 190 + 390
    min_threshold = 190 - 120

    # Find all indices where M is close to 1.4
    indices_close = findall(x -> isapprox(x, 1.4; atol=M_interval), sorted_M)

    # If no matching masses, return false
    if isempty(indices_close)
        return false
    end

    # Check if any corresponding T value satisfies the constraint
    for idx in indices_close
        if min_threshold <= sorted_T[idx] <= max_threshold
            return true  # At least one value satisfies the constraint
        end
    end

    return false  # No value satisfies the constraint
end

end  # End of the MainModule
