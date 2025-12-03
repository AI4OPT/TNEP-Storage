using Base.Threads

include("../structs/ExpansionPlanner.jl")

function get_benders_lb(superdir::String)
    """
    Computes a valid global lower bound by summing first-stage and second-stage LBs

    First-stage LB: computed by achieving zero load-shed feasibility on the "most challenging" representative day
    Second-stage LB: computed by taking the maximum parameters of all first-stage variables and optimizing operational cost
    """
    # TODO: check if calculation is necessary (already precomputed?)

    first_stage_lb = get_first_stage_lb(superdir)

    return
end

function get_first_stage_lb(superdir)
    """
    First-stage LB: computed by achieving zero load-shed feasibility on the "most challenging" representative day
    Solves <only_feasibility = true> for all representative days
    """

    

end