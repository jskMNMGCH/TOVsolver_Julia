module MainModule

# include files
include("solver_code.jl") # includeにより solver_code.jl の内容を現在のスコープに取り込む.
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
    mono_e = []
    mono_P = []
    for i in 2:length(e)
        dP = P[i] - P[i-1]
        de = e[i] - e[i-1]
        # The region where dp = 0 corresponds to the area during a first-order phase transition, so it should be removed.
        if dP > 0.0 && de > 0.0 
            push!(mono_e, e[i]) 
            push!(mono_P, P[i])
        end
    end
    return mono_e, mono_P
end

function out_RMT(ε, pres; debug=false, min_pc=2.9e-6, max_pc=3.8e-3, num_pc=100)
    """ 
    This function calculates the radius, mass, and tidal deformability by solving the Tolman–Oppenheimer–Volkoff (TOV) equations 
    for a given set of energy densities (`ε`) and pressures (`pres`).    
    Args:
        ε (array): Array of energy densities (dimensionless, scaled by ε_ref).
        pres (array): Array of pressures (dimensionless, scaled by ε_ref).
        debug (bool, optional): Flag to enable debugging mode. Default is `false`.
    
    Returns:
        Tuple: A tuple containing:
            - `R` (1D array): Array of radii corresponding to the input densities.
            - `M` (1D array): Array of masses corresponding to the input densities.
            - `Λ` (1D array): Array of tidal deformabilities corresponding to the input densities.
            - `sol_list` (1D array): List of solutions for each central density.
    """
    R = []
    M = []
    Λ = []
    sol_list = []
    # The dimensionless pressure 2.9e-6 corresponds to about 1.0 MeV/fm^3 in QHC18_gv100H16.
    # The dimensionless pressure 3.8e-3 corresponds to about 1300 MeV/fm^3 in QHC18_gv100H16.
    min_pc_id = searchsortedfirst(pres, min_pc)
    max_pc_id = min(searchsortedfirst(pres, max_pc), length(pres))
    len = min(max_pc_id-min_pc_id, num_pc)
    pc_indices = round.(Int, range(min_pc_id, max_pc_id, length=len))
    
    for i in pc_indices[end:-1:1]
        if i == pc_indices[end] && debug
            debug_flag = true
        else
            debug_flag = false
        end
        
        try
            result, sol = SolverCode.solveTOV_RMT(i, ε, pres, debug_flag)
            push!(R, result[1]) 
            push!(M, result[2])  
            push!(Λ, result[3])
            push!(sol_list, sol)
        catch e
            println("Error at i = $i")
        end
    end
    return [R, M, Λ], sol_list
end

end  # end of MainModule