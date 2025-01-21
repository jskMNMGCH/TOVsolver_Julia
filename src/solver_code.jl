"""
This code is written in Geometrized unit system, where distances are expressed in centimeters [cm].
If energy_density and pressure are given in [g/cm^3] unit, 
the inputs shoulds be energy_density/unit_g and pressure/unit_g [1/cm^3].
"""
module SolverCode

# Julia packages
using DifferentialEquations
using Dates
# include files
include("constants.jl")

# functions
function Debug(input_list; base_filename="debug/debug_data.txt")
"""
    This function appends the provided data to a timestamped log file.
    The file name will be dynamically generated by appending the current date and time to the provided `base_filename`.
    Args:
        input_list: A list of values (of any type) that you want to log.
        base_filename: The base name for the log file. The default is "debug/debug_data.txt". 
                       The current date and time will be appended to this base name.
"""
    now_str = Dates.format(now(), "yyyy_mm_dd_HH_MM")
    # 日付と時刻をファイル名に追加
    filename_with_date = replace(base_filename, ".txt" => "_$now_str.txt")
    open(filename_with_date, "a") do io
        for val in input_list
            write(io, string(val))
            write(io, ' ')
        end
        write(io, '\n')
    end
end

function tidal_deformability(y, M, R)
"""
    Calculates the tidal deformability for a given compactness and other parameters
    Args:
        y: A dimensionless parameter.
        M: The mass of the compact star (dimensionless).
        R: The radius of the compact star (dimensionless).
    Return:
        tidal_def: The tidal deformability.
"""
    C = M / R
    # println("Compactness ", C)
    Eps = (2.0 * C * (6.0 - 3.0 * y + 3.0 * C * (5.0 * y - 8.0)) + 4.0 * C^3 * (13.0 - 11.0 * y + C * (3.0 * y - 2.0) + 2.0 * C^2 * (1.0 + y)) + 3.0 * (1.0 - 2.0 * C)^2 * (2.0 - y + 2.0 * C * (y - 1.0))*log(1.0 - 2.0 * C))
    
    tidal_def = 16.0 / (15.0 * Eps) * (1.0 - 2.0 * C)^2 * (2.0 + 2.0 * C * (y - 1.0) - y)
    #println("(Eps, y, tidal) = ($(Eps), $(y), $(tidal_def))")
    return tidal_def
end


function TOV_def!(du, u, p, t)
"""
    This function defines the system of differential equations for the Tolman–Oppenheimer–Volkoff (TOV) equations. 
    The function performs linear interpolation between the nearest data points for energy density (`E`) and pressure (`Pres`) based on the current pressure value (`pres`). 
    The derivative (`f`) is calculated as the inverse of the sound speed squared (`f = 1 / c_s^2`).   
    Args:
        du (array): The computed derivatives.
        u (array): The current values.
        p (tuple): A tuple containing the energy density (`E`), pressure array (`Pres`), and a debug flag (`debug_flag`).
        t (float): The current radial coordinate.
    Returns:
        Nothing. The `du` array is updated in-place with the calculated derivatives.
    Side effects:
        If `debug_flag` is `true`, a debug output is generated by calling the `Debug` function with the current values of 'eps', 'pres' and 'f'.
"""
    E, Pres, debug_flag = p
    pres, m, h, b = u
    
    idx_now = searchsortedfirst(Pres, pres)
    if idx_now >= length(Pres)
        println("idx_end_point appear: ", idx_now, ", Pressure: ", pres*eps_ref*gcm3_to_MeVfm3, " [MeV/fm^3]")
        idx_now = length(Pres)-1
    elseif idx_now <= 1
        println("idx_end_point appear: ", idx_now, ", Pressure: ", pres*eps_ref*gcm3_to_MeVfm3, " [MeV/fm^3]")
        idx_now = 2
    end
    s = (pres - Pres[idx_now-1])/(Pres[idx_now]-Pres[idx_now-1])
    eps = (s*E[idx_now] + (1.0-s)*E[idx_now-1])  # [cm^-2]
    f = (E[idx_now]-E[idx_now-1])/(Pres[idx_now]-Pres[idx_now-1]) # 後進差分でfが計算されることを前提に上2行で線形補完している. 
    # f = ()
    
    dpdr = -(eps + pres) * (m + 4.0 * π * t^3 * pres) / (t*m *(t/m - 2.0))  # [cm^-3]
    dmdr = 4.0 * π * t^2 * eps  # []
    dhdr = b  # [cm]
    dbdr = ( 2.0*(1.0-2.0*m/t)^-1*h*(-2.0*π*(5.0*eps + 9.0*pres + f*(eps+pres)) 
        + 3.0/t^2 + 2.0*(1.0-2.0*m/t)^-1*(1/t^2)*(m/t + 4.0*π*t^2*pres)^2) 
        + 2.0*b/t * (1.0-2.0*m/t)^-1*(-1.0 + m/t + 2.0*π*t^2*eps-pres) )# []

    if debug_flag
        Debug([eps, pres, f])
    end
    
    du[1] = dpdr
    du[2] = dmdr
    du[3] = dhdr
    du[4] = dbdr
end


function solveTOV_RMT(center_idx, ε, pres, debug_flag)
"""
    This function solves the Tolman–Oppenheimer–Volkoff (TOV) equations for a given central density and pressure 
    to calculate the radius, mass, and tidal deformability of a compact object. 
    The function iterates over a radial grid, solving the TOV equations numerically with an ODE solver and returns the radius, mass, and tidal deformability.
    Args:
        center_idx (int): The index in the density and pressure arrays, corresponding to the central density and pressure.
        ε (array): Array of energy densities ([1/cm^2], [g/cm^3] is scaled by eps_ref [g/cm]).
        pres (array): Array of pressures ([1/cm^2], [g/cm^3] is scaled by eps_ref [g/cm]).
        debug_flag (bool): A flag to enable debugging output.
    Returns:
        Tuple:
            - `R` (float): The radius of the object in the unit of [km].
            - `M` (float): The mass of the object in the solar masses unit (dimensionless).
            - `Λ` (float): The tidal deformability of the object (dimensionless).
            - `sol` (array): The solution of the ODE at each step.
""" 
    # set initial state and parameters
    r = 10^-17/unit_l  # dimless length. initial radius.
    dr = 100.0/unit_l  # step value
    dhdr = dr
    
    p0 = pres[center_idx]
    p_surface = max(p0/10^12, pres[2]) # max で データより小さくならないようにするべき
    # println("Pres at the surface ", p_surface)
    ε0 = ε[center_idx]
    m0 = 4.0/3.0*π*ε0*r^3
    h0 = r^2
    b0 = 2.0 * r
    u0 = [p0, m0, h0, b0]
    p = (ε, pres, debug_flag)
    tspan = (r, Inf)

    # コールバックの設定
    condition(u, t, integrator) = u[1] > p_surface # 条件を満たす時のみ積分を進める
    affect!(integrator) = terminate!(integrator)  # 終了命令
    cb = ContinuousCallback(condition, affect!)
    # 問題を定義
    prob = ODEProblem(TOV_def!, u0, tspan, p)
    # 数値積分を実行
    sol = solve(prob, DP5(); dt=dr, callback=cb, dtmax=dr, maxiters=1e6) # maxiters=1e5 (default)

    # 最終状態を取得
    final_state = sol[end]
    R = sol.t[end]
    M = final_state[2]

    # calculation of the ε_surface by the linear interpolation.
    idx_surface = searchsortedfirst(pres, final_state[1])
    if idx_surface <= 1
        idx_surface = 2
    end
    ratio = (final_state[1] - pres[idx_surface-1])/(pres[idx_surface]-pres[idx_surface-1])
    ε_surface = (ratio*ε[idx_surface] + (1.0-ratio)*ε[idx_surface-1])
    # println("effect of the nonzero ε at the surface: ",  (4.0*π*R^3*ε_surface/M)/(R*final_state[4]/final_state[3]))
    
    # Tidal deformability の計算
    y = R * final_state[4] / final_state[3] - 4.0*π*R^3*ε_surface/M  # if EoS has non-zero ε, the 2nd term of y is important!
    Λ = tidal_deformability(y, M, R)
    return [R*unit_l/10^5, M*unit_g/Msun, Λ], sol
end

end  # end of SolverCode