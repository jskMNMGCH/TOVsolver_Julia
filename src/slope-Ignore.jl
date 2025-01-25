module Slope

include("constants.jl") # includeにより solver_code.jl の内容を現在のスコープに取り込む.

function get_id_eps_j(ε_joint, energy_density)
    id = searchsortedfirst(energy_density, ε_joint)
    if id < 1
        id =  1
        println("id < 1")
    elseif id > length(energy_density)
        id = length(energy_density)
        println("id > length(energy_density)")
    end
    return id, energy_density[id]
end

function cs2_at_id(id::Int, energy_density, pressure)
    cs2 = (pressure[id] - pressure[id-1])/(energy_density[id]-energy_density[id-1])
    return cs2
end

function alpha_K_of_slope(ε_joint_id::Int, ε_cut::Float64, energy_density, pressure)
    cs2 = cs2_at_id(ε_joint_id, energy_density, pressure)
    ε_joint = energy_density[ε_joint_id]
    α = cs2*(ε_joint-ε_cut)/pressure[ε_joint_id]
    K = pressure[ε_joint_id]*(ε_joint-ε_cut)^(-α)
    return α, K
end

function eos_with_slope(id::Int, ε_cut::Float64, energy_density, pressure; num_eps_slope::Float64=2.0^9)
    α, K = alpha_K_of_slope(id, ε_cut, energy_density, pressure)
    ε_joint = energy_density[id]
    ε_slope_len = ceil(Int, (1-ε_cut/ε_joint)*num_eps_slope)
    ε_slope_part = logrange(ε_cut, ε_joint, length=ε_slope_len)[1:end-1]
    p_slope_part = K.*(ε_slope_part.-ε_cut).^α
    whole_ε = vcat(ε_slope_part, energy_density[id:end])
    whole_p = vcat(p_slope_part, pressure[id:end])
    return whole_ε, whole_p
end

function closest_to_zero_indices(arr::AbstractVector{<:Number}, k::Int)
    """
    Returns the indices of the `k` elements in the array `arr` that are closest to 0.

    Args:
        arr: A vector of numbers (can be non-monotonic).
        k: The number of indices to return.

    Returns:
        Vector{Int}: A vector containing the indices of the `k` elements closest to 0.
    """
    # Check if the input array is empty
    if isempty(arr)
        error("The input array is empty.")
    end

    # Check if k is valid
    if k < 1 || k > length(arr)
        error("`k` must be between 1 and the length of the array.")
    end

    # Compute the absolute differences from 0
    diffs = abs.(arr)

    # Get the indices of the `k` smallest elements
    return partialsortperm(diffs, 1:k)
end

function shortest_distance_with_indices(curve1::Vector{Vector{Float64}}, curve2::Vector{Vector{Float64}})
    """
    Calculate the shortest distance between two 2D curves and return the indices of the points
    that achieve this shortest distance.

    Args:
        curve1::Vector{Vector{Float64}}: [[x1, x2, ...], [y1, y2, ...]] representing the first curve.
        curve2::Vector{Vector{Float64}}: [[x1, x2, ...], [y1, y2, ...]] representing the second curve.

    Returns:
        Tuple{Float64, Int, Int}: A tuple containing:
            - The shortest distance (Float64),
            - The index of the point on curve1,
            - The index of the point on curve2.
    """
    # Validate input dimensions
    if length(curve1) != 2 || length(curve2) != 2
        throw(ArgumentError("Both curves must be in the form [[x1, x2, ...], [y1, y2, ...]]."))
    end

    # Extract points
    x1, y1 = curve1[1], curve1[2]
    x2, y2 = curve2[1], curve2[2]

    # Ensure x and y have the same length
    if length(x1) != length(y1) || length(x2) != length(y2)
        throw(ArgumentError("The number of x and y coordinates must match for each curve."))
    end

    # Initialize variables to store the shortest distance and indices
    min_dist_squared = Inf  # Start with a large value
    min_index1 = -1         # Index for curve1
    min_index2 = -1         # Index for curve2

    # Compute the squared distances between all pairs of points
    for i in 1:length(x1)
        for j in 1:length(x2)
            dist_squared = (x1[i] - x2[j])^2 + (y1[i] - y2[j])^2
            if dist_squared < min_dist_squared
                min_dist_squared = dist_squared
                min_index1 = i  # Update index for curve1
                min_index2 = j  # Update index for curve2
            end
        end
    end

    # Return the shortest distance and indices
    return sqrt(min_dist_squared), min_index1, min_index2
end





end  # End of the Slope