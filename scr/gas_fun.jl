# gas functions
function load_gas_data(set)
    """
    extracts, transforms and loads gas network data
    """
    prod_data = CSV.read("data/$(set[:gas_case])/gas_prod.csv", DataFrame; header=1)
    node_data = CSV.read("data/$(set[:gas_case])/gas_node.csv", DataFrame; header=1)
    pipe_data = CSV.read("data/$(set[:gas_case])/gas_edge.csv", DataFrame; header=1)
    load_fact = CSV.read("data/demand_factor.csv", DataFrame; header=1, limit=set[:T])

    # sets
    E  = collect(1:length(pipe_data[!,:w]))
    E_a = vcat(findall(x->x!=0, pipe_data[!,:b]))
    N  = collect(1:length(node_data[!,:demand]))
    # gas producer data
    c1 = Array(prod_data[!,:c1])
    c2 = Array(prod_data[!,:c2])
    c̀2 = sqrt.(c2)
    ϑ̅ = Array(prod_data[!,:p_max])
    ϑ̲ = Array(prod_data[!,:p_min])
    # node data
    δ = Array(node_data[!,:demand])
    ρ̅ = Array(node_data[!,:presh_max])
    ρ̲ = Array(node_data[!,:presh_min])
    # edge data
    n_s = Array(pipe_data[!,:n_s])
    n_r = Array(pipe_data[!,:n_r])
    w = Array(pipe_data[!,:w])
    s = Array(pipe_data[!,:K_h]) * .5
    b = Array(pipe_data[!,:b])
    κ̅ = Array(pipe_data[!,:kappa_max])
    κ̲ = Array(pipe_data[!,:kappa_min])
    # reference pressure node
    ref = 26
    # node-edge incidence matrix
    A = zeros(size(E,1),size(N,1))
    A_pl = zeros(size(E,1),size(N,1))
    A_mn = zeros(size(E,1),size(N,1))
    for l in E
        A[l,n_s[l]] = 1
        A_pl[l,n_s[l]] = 1
        A[l,n_r[l]] = -1
        A_mn[l,n_r[l]] = -1
    end
    # gas - pressure conversion matrix
    B = zeros(size(N,1),size(E,1))
    for j in E
        B[n_s[j],j] = b[j]
    end
    # initial linepack
    # ψ_o = Array(pipe_data[!,:psi_o])
    ψ_o = zeros(length(E))
    for i in E
        ψ_o[i] = s[i]*(node_data[n_s[i],:presh_init] + node_data[n_r[i],:presh_init]) * 0.05
    end
    # save gas system data
    net_data = Dict(:c1 => c1, :c2 => c2, :c̀2 => c̀2, :ϑ̅ => ϑ̅, :ϑ̲ => ϑ̲,
                    :δ => δ, :k_δ => load_fact[:,1], :ρ̅ => ρ̅, :ρ̲ => ρ̲,
                    :n_s => n_s, :n_r => n_r, :w => w, :ref => ref, :s => s./5,
                    :ψ_o => ψ_o, :κ̅ => κ̅, :κ̲ => κ̲, :b => b, :B => B,
                    :A => A', :A⁺ => A_pl', :A⁻ => A_mn',
                    :E => E, :E_a => E_a, :N => N)
    @info("done loading gas network data")
    return net_data
end

function non_convex_opt(net)
    """
    solves non-convex gas network optimization
    """
    # build model
    model = Model(optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 15000, "print_level" => 0))
    # variable declaration
    @variable(model, ϱ[1:length(net[:N]),t=1:set[:T]])
    @variable(model, φ[1:length(net[:E]),t=1:set[:T]])
    @variable(model, κ[1:length(net[:E]),t=1:set[:T]])
    @variable(model, ϑ[1:length(net[:N]),t=1:set[:T]])
    @variable(model, φ⁻[1:length(net[:E]),t=1:set[:T]])
    @variable(model, φ⁺[1:length(net[:E]),t=1:set[:T]])
    @variable(model, ψ[1:length(net[:E]),t=1:set[:T]])
    # minimize gas production cost
    @objective(model, Min, sum(sum(net[:c1][i]*ϑ[i,t] + net[:c2][i]*ϑ[i,t]^2 for i in net[:N]) for t in 1:set[:T]))
    # gas variable limits
    @constraint(model, inj_lim_max[i=net[:N],t=1:set[:T]], ϑ[i,t] <= net[:ϑ̅][i])
    @constraint(model, inj_lim_min[i=net[:N],t=1:set[:T]], ϑ[i,t] >= net[:ϑ̲][i])
    @constraint(model, pre_lim_max[i=net[:N],t=1:set[:T]], ϱ[i,t] <= net[:ρ̅][i])
    @constraint(model, pre_lim_min[i=net[:N],t=1:set[:T]], ϱ[i,t] >= net[:ρ̲][i])
    @constraint(model, com_lim_max[i=net[:E],t=1:set[:T]], κ[i,t] <= net[:κ̅][i])
    @constraint(model, com_lim_min[i=net[:E],t=1:set[:T]], κ[i,t] >= net[:κ̲][i])
    @constraint(model, SoC_lim_min[i=net[:E],t=set[:T]  ], ψ[i,t] >= net[:ψ_o][i])
    # gas flow equations
    @NLconstraint(model, w_eq[l=net[:E],t=1:set[:T]],  φ[l,t]*abs(φ[l,t]) - net[:w][l]^2 *((ϱ[ns(l),t] + κ[l,t])^2 - ϱ[nr(l),t]^2) == 0)
    @constraint(model, gas_bal[i=net[:N],t=1:set[:T]],
        sum(net[:A⁺][i,j]*φ⁺[j,t] + net[:A⁻][i,j]*φ⁻[j,t] for j in net[:E])
        == ϑ[i,t] - sum(net[:b][l]*κ[l,t] for l in net[:E] if ns(l) == i) - (Δ(t)*S(t)*ξ[:μ̂])[i])
    @constraint(model, φ_pl[l=net[:E_a],t=1:set[:T]],  φ[l,t] >= 0)
    @constraint(model, φ_avr[l=net[:E],t=1:set[:T]],   φ[l,t] == 1/2*(φ⁺[l,t] + φ⁻[l,t]))
    @constraint(model, ψ_avr[l=net[:E],t=1:set[:T]],   ψ[l,t] == 1/2*net[:s][l]*(ϱ[ns(l),t] + κ[l,t] + ϱ[nr(l),t]))
    @constraint(model, LP_b_t[l=net[:E],t=2:set[:T]],  ψ[l,t] == ψ[l,t-1] +  φ⁺[l,t] - φ⁻[l,t])
    @constraint(model, LP_b_1[l=net[:E],t=1],                   ψ[l,t] == net[:ψ_o][l] +  φ⁺[l,t] - φ⁻[l,t])
    @info("solving the non-convex gas optimization model takes time, wait ...")
    # solve model
    optimize!(model)
    @info("done solving the non-convex gas optimization model in $(solve_time(model)) sec. with status: $(termination_status(model))")
    # return solution
    solution = Dict(
    :ϑ => JuMP.value.(ϑ),
    :ϱ => JuMP.value.(ϱ),
    :φ => JuMP.value.(φ),
    :κ => JuMP.value.(κ),
    :ψ => JuMP.value.(ψ),
    :cost => JuMP.objective_value.(model),
    :status => termination_status(model),
    :CPUtime => solve_time(model),
    :model => model)
    return solution
end

function linearization(net,model)
    """
    extracts the senstivity coefficients from the non-convex Weymouth equation
    """
    con_nl = model[:w_eq]
    var_nl = [model[:ϱ];model[:φ];model[:κ]]
    # aux functions
    raw_index(var::MOI.VariableIndex) = var.value
    raw_index(con::NonlinearConstraintIndex) = con.value
    # extract variables
    variables = all_variables(model)
    # extract Jacobian structure
    d = NLPEvaluator(model)
    MOI.initialize(d, [:Jac])
    jac_str = MOI.jacobian_structure(d)
    # extract operating point
    opertaing_point = value.(variables)
    # evaluate Weymouth eq. Jacobian at the operating point
    V = zeros(length(jac_str))
    MOI.eval_constraint_jacobian(d, V, opertaing_point)
    # prepare the result
    I = [js[1] for js in jac_str] # rows, equations
    J = [js[2] for js in jac_str] # cols, variables
    jac = sparse(I, J, V)
    rows = raw_index.(index.(con_nl))
    cols = raw_index.(index.(var_nl))
    # return evaluated Jacobians for each time period (as dictionaries)
    Jacobian = Dict(:ϱ => [], :φ => [], :κ => [])
    function return_jac(net,jac,rows,cols)
        Jac = Matrix(jac[rows, cols])
        Jac_ϱ = Jac[net[:E],net[:N]]
        Jac_φ = Jac[net[:E],length(net[:N]) .+ net[:E]]
        Jac_κ = Jac[net[:E],(length(net[:N]) .+ length(net[:E])) .+ net[:E]]
        return Jac_ϱ, Jac_φ, Jac_κ
    end
    for t in 1:set[:T]
          J_ϱ, J_φ, J_κ  = return_jac(net,jac,rows[:,t],cols[:,t])
          push!(Jacobian[:ϱ],J_ϱ); push!(Jacobian[:φ],J_φ); push!(Jacobian[:κ],J_κ)
    end
    # return operating points for each time period (as dictionaries)
    point = Dict(:ϱ => [], :φ => [], :κ => [])
    for t in 1:set[:T]
          push!(point[:ϱ],opertaing_point[raw_index.(index.(var_nl))][net[:N],t])
          push!(point[:φ],opertaing_point[raw_index.(index.(var_nl))][length(net[:N]) .+ net[:E],t])
          push!(point[:κ],opertaing_point[raw_index.(index.(var_nl))][length(net[:N]) .+ length(net[:E]) .+ net[:E],t])
    end
    # compute linearization coefficients
    lin_res = Dict(:W_0 => [], :W_1 => [], :W_2 => [], :point => point)
    for t in 1:set[:T]
          push!(lin_res[:W_0], inv(Jacobian[:φ][t]) * (Jacobian[:ϱ][t] * point[:ϱ][t] + Jacobian[:κ][t] * point[:κ][t] + Jacobian[:φ][t] * point[:φ][t]) )
          push!(lin_res[:W_1], -inv(Jacobian[:φ][t]) * Jacobian[:ϱ][t])
          push!(lin_res[:W_2], -inv(Jacobian[:φ][t]) * Jacobian[:κ][t])
    end
    @info("done extracting sensetivities and reference pressures")
    return  lin_res
end

function linearized_opt(net,lin_res)
    """
    solves non-convex gas network optimization
    """
    # build model
    model = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "max_iter" => 15000))
    # variable declaration
    @variable(model, ϱ[1:length(net[:N]),t=1:set[:T]])
    @variable(model, φ[1:length(net[:E]),t=1:set[:T]])
    @variable(model, κ[1:length(net[:E]),t=1:set[:T]])
    @variable(model, ϑ[1:length(net[:N]),t=1:set[:T]])
    @variable(model, φ⁻[1:length(net[:E]),t=1:set[:T]])
    @variable(model, φ⁺[1:length(net[:E]),t=1:set[:T]])
    @variable(model, ψ[1:length(net[:E]),t=1:set[:T]])
    # minimize gas production cost
    @objective(model, Min, sum(sum(net[:c1][i]*ϑ[i,t] + net[:c2][i]*ϑ[i,t]^2 for i in net[:N]) for t in 1:set[:T]))
    # gas variable limits
    @constraint(model, inj_lim_max[i=net[:N],t=1:set[:T]], ϑ[i,t] <= net[:ϑ̅][i])
    @constraint(model, inj_lim_min[i=net[:N],t=1:set[:T]], ϑ[i,t] >= net[:ϑ̲][i])
    @constraint(model, pre_lim_max[i=net[:N],t=1:set[:T]], ϱ[i,t] <= net[:ρ̅][i])
    @constraint(model, pre_lim_min[i=net[:N],t=1:set[:T]], ϱ[i,t] >= net[:ρ̲][i])
    @constraint(model, com_lim_max[i=net[:E],t=1:set[:T]], κ[i,t] <= net[:κ̅][i])
    @constraint(model, com_lim_min[i=net[:E],t=1:set[:T]], κ[i,t] >= net[:κ̲][i])
    @constraint(model, SoC_lim_min[i=net[:E],t=set[:T]  ], ψ[i,t] >= net[:ψ_o][i])
    # gas flow equations
    @constraint(model, w_eq[l=net[:E], t=1:set[:T]],   φ[l,t] .== lin_res[:W_0][t][l] .+ sum(lin_res[:W_1][t][l,n]*ϱ[n,t] for n in net[:N]) .+ sum(lin_res[:W_2][t][l,k]*κ[k,t] for k in net[:E]))
    # @constraint(model, gas_bal[i=net[:N],t=1:set[:T]], sum(net[:A⁺][i,j]*φ⁺[j,t] + net[:A⁻][i,j]*φ⁻[j,t] for j in net[:E]) == ϑ[i,t] - sum(net[:b][l]*κ[l,t] for l in net[:E] if ns(l) == i) - net[:δ][i]*net[:k_δ][t])
    @constraint(model, gas_bal[i=net[:N],t=1:set[:T]],
        sum(net[:A⁺][i,j]*φ⁺[j,t] + net[:A⁻][i,j]*φ⁻[j,t] for j in net[:E])
        == ϑ[i,t] - sum(net[:b][l]*κ[l,t] for l in net[:E] if ns(l) == i) - (Δ(t)*S(t)*ξ[:μ̂])[i])
    @constraint(model, φ_pl[l=net[:E_a],t=1:set[:T]],  φ[l,t] >= 0)
    @constraint(model, φ_avr[l=net[:E],t=1:set[:T]],   φ[l,t] == 1/2*(φ⁺[l,t] + φ⁻[l,t]))
    @constraint(model, ψ_avr[l=net[:E],t=1:set[:T]],   ψ[l,t] == 1/2*net[:s][l]*(ϱ[ns(l),t] + κ[l,t] + ϱ[nr(l),t]))
    @constraint(model, LP_b_t[l=net[:E],t=2:set[:T]],  ψ[l,t] == ψ[l,t-1] +  φ⁺[l,t] - φ⁻[l,t])
    @constraint(model, LP_b_1[l=net[:E],t=1],                   ψ[l,t] == net[:ψ_o][l] +  φ⁺[l,t] - φ⁻[l,t])
    @constraint(model, ref_ϱ[t=1:set[:T]],                  ϱ[net[:ref],t] == lin_res[:point][:ϱ][t][net[:ref]])
    # solve model
    optimize!(model)
    @info("done solving the linearized gas optimization model in $(solve_time(model)) sec. with status: $(termination_status(model))")
    # return solution
    solution = Dict(
    :ϑ => JuMP.value.(ϑ),
    :ϱ => JuMP.value.(ϱ),
    :φ => JuMP.value.(φ),
    :κ => JuMP.value.(κ),
    :ψ => JuMP.value.(ψ),
    :φ⁺ => JuMP.value.(φ⁺),
    :φ⁻ => JuMP.value.(φ⁻),
    :cost => JuMP.objective_value.(model),
    :status => termination_status(model),
    :model => model)
    return solution
end

function chance_constrained_gas_opt(net,lin_res,set)
    """
    solves the chance-constrained gas network optimization
    """
    # build model
    model = Model(optimizer_with_attributes(Mosek.Optimizer, "LOG" => 0))
    # gas network variable declaration
    P =  [@variable(model, [1:length(net[:N]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    Φ =  [@variable(model, [1:length(net[:E]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    K =  [@variable(model, [1:length(net[:E]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    Θ =  [@variable(model, [1:length(net[:N]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    Φ⁻ = [@variable(model, [1:length(net[:E]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    Φ⁺ = [@variable(model, [1:length(net[:E]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    Ψ =  [@variable(model, [1:length(net[:E]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    # auxilary variables for variability control
    s_ϱ =  [@variable(model, [1:length(net[:N])]) for t = 1:set[:T]]
    # auxilary variables for DRO constraints
    y_ϑ =  [@variable(model, [1:length(net[:N])]) for t = 1:set[:T]]
    x_ϑ =  [@variable(model, [1:length(net[:N])]) for t = 1:set[:T]]
    y_ρ =  [@variable(model, [1:length(net[:N])]) for t = 1:set[:T]]
    x_ρ =  [@variable(model, [1:length(net[:N])]) for t = 1:set[:T]]
    y_κ =  [@variable(model, [1:length(net[:E])]) for t = 1:set[:T]]
    x_κ =  [@variable(model, [1:length(net[:E])]) for t = 1:set[:T]]
    # minimize gas injection cost + pressure variability measure
    @objective(model, Min, sum(net[:c1]' * Θ[t] * S(t) * ξ[:μ̂]
                                + tr(Θ[t]' * diagm(net[:c2]) * Θ[t] * S(t) * ξ[:Σ̂] * S(t)')
                           for t in 1:set[:T])
                           # + set[:α_ϱ]*sum(tr((P[t]*S(t)-P[t-1]*S(t-1))*ξ[:Σ̂]*(P[t]*S(t)-P[t-1]*S(t-1))') for t in 2:set[:T])
                           + set[:α_ϱ] * sum(sum(s_ϱ[t]) for t in 2:set[:T])
            )
    # variance constraints
    @constraint(model, pre_var[t=2:set[:T], i=net[:N]],
        [s_ϱ[t][i];(cholesky(ξ[:Σ̂]).U .+ zeros(set[:k_t][set[:T]],set[:k_t][set[:T]]))*(P[t]*S(t) - P[t-1]*S(t-1))[i,:]] in SecondOrderCone())
    # injection deviation constraints
    @constraint(model, inj_lim[t=1:set[:T], i=net[:N]],
        [set[:α_ϑ]*(Θ[t]*S(t)*ξ[:μ̂])[i];
            (cholesky(ξ[:Σ̂]).U .+ zeros(set[:k_t][set[:T]],set[:k_t][set[:T]])) * (Θ[t]*S(t) - Θ[t]*S(t)*ξ[:μ̂]*S(1))[i,:]]
                in SecondOrderCone())
    # linepack flexibility limit constraints
    @constraint(model, lp_lim[t=1:set[:T], i=net[:E]],
        [set[:α_ψ]*(Ψ[t]*S(t)*ξ[:μ̂])[i];
            (cholesky(ξ[:Σ̂]).U .+ zeros(set[:k_t][set[:T]],set[:k_t][set[:T]])) * (Ψ[t]*S(t) - Ψ[t]*S(t)*ξ[:μ̂]*S(1))[i,:]]
                in SecondOrderCone())
    # equivalents of stochastic equations
    @constraint(model, gas_bal[t=1:set[:T]], (net[:A⁺]*Φ⁺[t]*S(t) + net[:A⁻]*Φ⁻[t]*S(t) - Θ[t]*S(t) + net[:B]*K[t]*S(t) + Δ(t)*S(t)) .== 0)
    @constraint(model, lin_wey[t=1:set[:T]], (Φ[t]*S(t) - [lin_res[:W_0][t] zeros(length(net[:E]),set[:k_t][t]-1)]*S(t) - lin_res[:W_1][t]*P[t]*S(t) - lin_res[:W_2][t]*K[t]*S(t)) .== 0)
    @constraint(model, ref_pre[t=1:set[:T]], (P[t][net[:ref],:]' - [lin_res[:point][:ϱ][t][net[:ref]] zeros(1,set[:k_t][t]-1)])*S(t) .== 0)
    @constraint(model, av_flow[t=1:set[:T]], (Φ[t]*S(t) - 1/2*Φ⁺[t]*S(t) - 1/2*Φ⁻[t]*S(t)) .== 0)
    @constraint(model, l_p_def[t=1:set[:T]], (Ψ[t]*S(t) - 1/2*diagm(net[:s])*(K[t] + abs.(net[:A])'*P[t])*S(t)) .== 0)
    @constraint(model, l_p_bal_1[t=1], (Ψ[t]*S(t) - Φ⁺[t]*S(t) + Φ⁻[t]*S(t)) .== 0)
    @constraint(model, l_p_bal_t[t=2:set[:T]], (Ψ[t]*S(t) - Ψ[t-1]*S(t-1) - Φ⁺[t]*S(t) + Φ⁻[t]*S(t)) .== 0)
    # reformulated chance constraints
    # double-sided gas injection chance constraints
    @constraint(model, ϑ_1[t=1:set[:T], i=1:length(net[:N])],
        [sqrt(set[:ε̅_g]) * (1/2*(net[:ϑ̅][i] - net[:ϑ̲][i]) - x_ϑ[t][i]);vcat(ξ[:Σ̅][2:end,2:end] * (Θ[t]*S(t))[i,2:end],y_ϑ[t][i])] in SecondOrderCone())
    @constraint(model, ϑ_2[t=1:set[:T], i=1:length(net[:N])],
        (Θ[t]*S(t))[i,1] - 1/2*(net[:ϑ̅][i] + net[:ϑ̲][i]) <= y_ϑ[t][i] + x_ϑ[t][i])
    @constraint(model, ϑ_3[t=1:set[:T], i=1:length(net[:N])],
        (Θ[t]*S(t))[i,1] - 1/2*(net[:ϑ̅][i] + net[:ϑ̲][i]) >= - y_ϑ[t][i] - x_ϑ[t][i])
    @constraint(model, ϑ_4[t=1:set[:T], i=1:length(net[:N])],
        1/2*(net[:ϑ̅][i] - net[:ϑ̲][i]) >= x_ϑ[t][i] >= 0)
    @constraint(model, ϑ_5[t=1:set[:T], i=1:length(net[:N])],
        y_ϑ[t][i] >= 0)
    # double-sided pressure chance constraints
    @constraint(model, ρ_1[t=1:set[:T], i=1:length(net[:N])],
        [sqrt(set[:ε̅_g]) * (1/2*(net[:ρ̅][i] - net[:ρ̲][i]) - x_ρ[t][i]);vcat(ξ[:Σ̅][2:end,2:end] * (P[t]*S(t))[i,2:end],y_ρ[t][i])] in SecondOrderCone())
    @constraint(model, ρ_2[t=1:set[:T], i=1:length(net[:N])],
        (P[t]*S(t))[i,1] - 1/2*(net[:ρ̅][i] + net[:ρ̲][i]) <= y_ρ[t][i] + x_ρ[t][i])
    @constraint(model, ρ_3[t=1:set[:T], i=1:length(net[:N])],
        (P[t]*S(t))[i,1] - 1/2*(net[:ρ̅][i] + net[:ρ̲][i]) >= - y_ρ[t][i] - x_ρ[t][i])
    @constraint(model, ρ_4[t=1:set[:T], i=1:length(net[:N])],
        1/2*(net[:ρ̅][i] - net[:ρ̲][i]) >= x_ρ[t][i] >= 0)
    @constraint(model, ρ_5[t=1:set[:T], i=1:length(net[:N])],
        y_ρ[t][i] >= 0)
    # double-sided pressure regulation chance constraints
    @constraint(model, κ_1[t=1:set[:T], i=1:length(net[:E])],
        [sqrt(set[:ε̅_g]) * (1/2*(net[:κ̅][i] - net[:κ̲][i]) - x_κ[t][i]);vcat(ξ[:Σ̅][2:end,2:end] * (K[t]*S(t))[i,2:end],y_κ[t][i])] in SecondOrderCone())
    @constraint(model, κ_2[t=1:set[:T], i=1:length(net[:E])],
        (K[t]*S(t))[i,1] - 1/2*(net[:κ̅][i] + net[:κ̲][i]) <= y_κ[t][i] + x_κ[t][i])
    @constraint(model, κ_3[t=1:set[:T], i=1:length(net[:E])],
        (K[t]*S(t))[i,1] - 1/2*(net[:κ̅][i] + net[:κ̲][i]) >= - y_κ[t][i] - x_κ[t][i])
    @constraint(model, κ_4[t=1:set[:T], i=1:length(net[:E])],
        1/2*(net[:κ̅][i] - net[:κ̲][i]) >= x_κ[t][i] >= 0)
    @constraint(model, κ_5[t=1:set[:T], i=1:length(net[:E])],
        y_κ[t][i] >= 0)
    # chance constraints on linepack and flow variables
    @constraint(model, ψ_T[t=set[:T], i=1:length(net[:E])],
        [(Ψ[t]*S(t)*ξ[:μ̂] - net[:ψ_o])[i]; z(set[:ε̅_g])*ξ[:Σ̅]*(Ψ[t]*S(t))[i,:]]
            in SecondOrderCone())
    @constraint(model, φ_pl[t=1:set[:T], i=net[:E_a]],
        [(Φ[t]*S(t)*ξ[:μ̂] .- 0)[i]; z(set[:ε̅_g])*ξ[:Σ̅]*(Φ[t]*S(t))[i,:]]
            in SecondOrderCone())
    # zero recourse for nodes and edges with no controllable components
    @constraint(model, zero_inj[t=1:set[:T], i=findall(x->x==0, net[:ϑ̅])], Θ[t][i,2:end] .== 0)
    @constraint(model, zero_reg[t=1:set[:T], i=setdiff(findall(x->x==0, net[:κ̅]),findall(x->x<0, net[:κ̲]))],K[t][i,2:end] .== 0)
    # solve model
    optimize!(model)
        @info("done solving chance-constrained gas network optimization in $(solve_time(model)) sec. with status $(termination_status(model))")
    # return solution
    solution = Dict(
    :P => [JuMP.value.(P[t]) for t in 1:set[:T]],
    :Φ => [JuMP.value.(Φ[t]) for t in 1:set[:T]],
    :K => [JuMP.value.(K[t]) for t in 1:set[:T]],
    :Θ => [JuMP.value.(Θ[t]) for t in 1:set[:T]],
    :Φ⁻ => [JuMP.value.(Φ⁻[t]) for t in 1:set[:T]],
    :Φ⁺ => [JuMP.value.(Φ⁺[t]) for t in 1:set[:T]],
    :Ψ => [JuMP.value.(Ψ[t]) for t in 1:set[:T]],
    :cost => sum(net[:c1]' * JuMP.value.(Θ[t]) * S(t) * ξ[:μ̂] + tr(JuMP.value.(Θ[t])' * diagm(net[:c2]) * JuMP.value.(Θ[t]) * S(t) * ξ[:Σ̂] * S(t)') for t in 1:set[:T]),
    :ϱ_std_tot => sum(sum(JuMP.value.(s_ϱ[t])) for t in 2:1:set[:T]) ,
    :status => termination_status(model),
    :obj => JuMP.objective_value(model),
    :model => model)
    return solution
end

function export_results(set,ξ,sol_lin,sol_sto_wo_var,sol_sto_wt_var)
    """
    exports results
    """
    ## wind data
    wind_results = zeros(1+set[:Ns],set[:T])
    wind_results[1,:] =  [sum(W(t)*S(t)*ξ[:μ̂]) for t in 1:set[:T]]
    for i in 1:set[:Ns]
        wind_results[i+1,:] = [sum(W(t)*S(t)*ξ[:ξ̂][:,i]) for t in 1:set[:T]]
    end
    ## gas consumption data
    gas_coms_results = zeros(1+set[:Ns],set[:T])
    gas_coms_results[1,:] =  [sum(Δ(t)*S(t)*ξ[:μ̂]) for t in 1:set[:T]]
    for i in 1:set[:Ns]
        gas_coms_results[i+1,:] = [sum(Δ(t)*S(t)*ξ[:ξ̂][:,i]) for t in 1:set[:T]]
    end
    ## first-stage linepack difference between det and stoch solutions
    δ_Ψ = (sol_sto_wo_var[:Ψ][1] .- sol_lin[:ψ][:,1])./sol_lin[:ψ][:,1] * 100
    δ_Ψ_n = zeros(length(δ_Ψ)); α = 10;
    for i in 1:51 # normalization
        δ_Ψ_n[i] = 1 - 1/2 * (1 - δ_Ψ[i] / sqrt(δ_Ψ[i]^2 + α))
    end
    ## pressure results
    det_pressure = [sum(sol_lin[:ϱ][:,t]) for t in 1:set[:T]]
    sto_pressure_wo_var = zeros(1+set[:Ns],set[:T])
    sto_pressure_wt_var = zeros(1+set[:Ns],set[:T])
    sto_pressure_wo_var[1,:] = [sum(sol_sto_wo_var[:P][t]*S(t)*ξ[:μ̂]) for t in 1:set[:T]]
    sto_pressure_wt_var[1,:] = [sum(sol_sto_wt_var[:P][t]*S(t)*ξ[:μ̂]) for t in 1:set[:T]]
    for i in 1:set[:Ns]
        sto_pressure_wo_var[i+1,:] = [sum(sol_sto_wo_var[:P][t]*S(t)*ξ[:ξ̂][:,i]) for t in 1:set[:T]]
        sto_pressure_wt_var[i+1,:] = [sum(sol_sto_wt_var[:P][t]*S(t)*ξ[:ξ̂][:,i]) for t in 1:set[:T]]
    end
    ## linepack results
    det_linepack = [sum(sol_lin[:ψ][:,t]) for t in 1:set[:T]]
    sto_linepack_wo_var = zeros(1+set[:Ns],set[:T])
    sto_linepack_wt_var = zeros(1+set[:Ns],set[:T])
    sto_linepack_wo_var[1,:] = [sum(sol_sto_wo_var[:Ψ][t]*S(t)*ξ[:μ̂]) for t in 1:set[:T]]
    sto_linepack_wt_var[1,:] = [sum(sol_sto_wt_var[:Ψ][t]*S(t)*ξ[:μ̂]) for t in 1:set[:T]]
    for i in 1:set[:Ns]
        sto_linepack_wo_var[i+1,:] = [sum(sol_sto_wo_var[:Ψ][t]*S(t)*ξ[:ξ̂][:,i]) for t in 1:set[:T]]
        sto_linepack_wt_var[i+1,:] = [sum(sol_sto_wt_var[:Ψ][t]*S(t)*ξ[:ξ̂][:,i]) for t in 1:set[:T]]
    end
    ## presusre regulation
    det_pres_reg = [sum(abs.(sol_lin[:κ][:,t])) for t in 1:set[:T]]
    sto_pres_reg_wo_var = zeros(1+set[:Ns],set[:T])
    sto_pres_reg_wt_var = zeros(1+set[:Ns],set[:T])
    sto_pres_reg_wo_var[1,:] = [sum(abs.(sol_sto_wo_var[:K][t]*S(t)*ξ[:μ̂])) for t in 1:set[:T]]
    sto_pres_reg_wt_var[1,:] = [sum(abs.(sol_sto_wt_var[:K][t]*S(t)*ξ[:μ̂])) for t in 1:set[:T]]
    for i in 1:set[:Ns]
        sto_pres_reg_wo_var[i+1,:] = [sum(abs.(sol_sto_wo_var[:K][t]*S(t)*ξ[:ξ̂][:,i])) for t in 1:set[:T]]
        sto_pres_reg_wt_var[i+1,:] = [sum(abs.(sol_sto_wt_var[:K][t]*S(t)*ξ[:ξ̂][:,i])) for t in 1:set[:T]]
    end

    CSV.write("results/wind_results.csv",  Tables.table([[i for i in 1:set[:T]]' ; wind_results./maximum(wind_results)]'), writeheader=false)
    CSV.write("results/gas_coms_results.csv",  Tables.table([[i for i in 1:set[:T]]' ; gas_coms_results./maximum(gas_coms_results)]'), writeheader=false)
    CSV.write("results/linepack_diff.csv",  Tables.table(δ_Ψ_n'), writeheader=false)

    base_value = max(maximum(det_pressure),maximum(sto_pressure_wo_var),maximum(sto_pressure_wt_var))
    CSV.write("results/det_pressure.csv",  Tables.table([[i for i in 1:set[:T]]' ; det_pressure'./base_value]'), writeheader=false)
    CSV.write("results/sto_pressure_wo_var.csv",  Tables.table([[i for i in 1:set[:T]]' ; sto_pressure_wo_var./base_value]'), writeheader=false)
    CSV.write("results/sto_pressure_wt_var.csv",  Tables.table([[i for i in 1:set[:T]]' ; sto_pressure_wt_var./base_value]'), writeheader=false)

    base_value = max(maximum(det_linepack),maximum(sto_linepack_wo_var),maximum(sto_linepack_wt_var))
    CSV.write("results/det_linepack.csv",  Tables.table([[i for i in 1:set[:T]]' ; det_linepack'./base_value]'), writeheader=false)
    CSV.write("results/sto_linepack_wo_var.csv",  Tables.table([[i for i in 1:set[:T]]' ; sto_linepack_wo_var./base_value]'), writeheader=false)
    CSV.write("results/sto_linepack_wt_var.csv",  Tables.table([[i for i in 1:set[:T]]' ; sto_linepack_wt_var./base_value]'), writeheader=false)

    base_value = max(maximum(det_pres_reg),maximum(sto_pres_reg_wo_var),maximum(sto_pres_reg_wt_var))
    CSV.write("results/det_pres_reg.csv",  Tables.table([[i for i in 1:set[:T]]' ; det_pres_reg'./base_value]'), writeheader=false)
    CSV.write("results/sto_pres_reg_wo_var.csv",  Tables.table([[i for i in 1:set[:T]]' ; sto_pres_reg_wo_var./base_value]'), writeheader=false)
    CSV.write("results/sto_pres_reg_wt_var.csv",  Tables.table([[i for i in 1:set[:T]]' ; sto_pres_reg_wt_var./base_value]'), writeheader=false)
end

function topology_conf(set)
    """
    Identifies sensetivity coefficients and reference pressure
    for different gas network topologies
    """
    # number of valves
    v_num = length(set[:b_valve_edge])
    # number of topology configurations
    C = 2^v_num
    # binary configuration list
    con_list = all_config(v_num)
    # compute senstivities and ref. pressure for each configuration
    output = Dict()
    @info("star enumerating candidate topologies")
    for i in 1:C
        # modify network topology
        net_data = load_gas_data(set)
        for j in 1:v_num
            net_data[:w][set[:b_valve_edge][j]] = net_data[:w][set[:b_valve_edge][j]] * con_list[i,j]
        end
        # run linearization
        sol_non_convex  = non_convex_opt(net_data)
        lin_res         = linearization(net_data,sol_non_convex[:model])
        # save results
        conf_data = Dict(:W_0 => [lin_res[:W_0][t] for t in 1:set[:T]],
                         :W_1 => [lin_res[:W_1][t] for t in 1:set[:T]],
                         :W_2 => [lin_res[:W_2][t] for t in 1:set[:T]],
                         :point => Dict(:ϱ => [lin_res[:point][:ϱ][t] for t in 1:set[:T]])
        )
        push!(output, i => conf_data)
    end
    @info("done enumerating candidate topologies")
    return output
end

function all_config(n)
    """
    Finds the list of possible network topologies
    """
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    @variable(model, ξ[1:n])
    @constraint(model, con[i=1:n], 0 <= ξ[i] <= 1)
    # create a polyhedron from the sample set
    poly = polyhedron(model, CDDLib.Library(:exact))
    # obtain V-representation of the polyhedron
    vr = vrep(poly)
    vr = MixedMatVRep(vr)
    p = vr.V
    return float(p)
end

function out_of_sample(solution,net,set,ξ)
    """
    Performs out-of-sample analysis
    """
    ϑ = zeros(length(net[:N]),set[:T],set[:Ns])
    ϱ = zeros(length(net[:N]),set[:T],set[:Ns])
    ψ = zeros(length(net[:E]),set[:T],set[:Ns])
    ϕ = zeros(length(net[:E]),set[:T],set[:Ns])
    κ = zeros(length(net[:E]),set[:T],set[:Ns])

    ϑ_mag =  zeros(set[:T])
    ϱ_mag =  zeros(set[:T])
    ψ_mag =  zeros(set[:T])
    ϕ_mag =  zeros(set[:T])
    κ_mag =  zeros(set[:T])

    ϑ_cvar =  zeros(set[:T],set[:Ns])
    ϱ_cvar =  zeros(set[:T],set[:Ns])
    ψ_cvar =  zeros(set[:T],set[:Ns])
    ϕ_cvar =  zeros(set[:T],set[:Ns])
    κ_cvar =  zeros(set[:T],set[:Ns])
    for t in 1:set[:T], s in 1:set[:Ns]
         ϑ[:,t,s] = solution[:Θ][t]*S(t)*ξ[:ξ̂][:,s]
         ϱ[:,t,s] = solution[:P][t]*S(t)*ξ[:ξ̂][:,s]
         ψ[:,t,s] = solution[:Ψ][t]*S(t)*ξ[:ξ̂][:,s]
         ϕ[:,t,s] = solution[:Φ][t]*S(t)*ξ[:ξ̂][:,s]
         κ[:,t,s] = solution[:K][t]*S(t)*ξ[:ξ̂][:,s]
    end
    for t in 1:set[:T], i in 1:48, s in 1:set[:Ns]
        (ϑ[i,t,s] - net[:ϑ̅][i]) >= 0 ? ϑ_mag[t] += (ϑ[i,t,s] - net[:ϑ̅][i]) : NaN
        (net[:ϑ̲][i] - ϑ[i,t,s]) >= 0 ? ϑ_mag[t] += (net[:ϑ̲][i] - ϑ[i,t,s]) : NaN
        (ϑ[i,t,s] - net[:ϑ̅][i]) >= 0 ? ϑ_cvar[t,s] += (ϑ[i,t,s] - net[:ϑ̅][i]) : NaN
        (net[:ϑ̲][i] - ϑ[i,t,s]) >= 0 ? ϑ_cvar[t,s] += (net[:ϑ̲][i] - ϑ[i,t,s]) : NaN

        (ϱ[i,t,s] - net[:ρ̅][i]) >= 0 ? ϱ_mag[t] += (ϱ[i,t,s] - net[:ρ̅][i]) : NaN
        (net[:ρ̲][i] - ϱ[i,t,s]) >= 0 ? ϱ_mag[t] += (net[:ρ̲][i] - ϱ[i,t,s]) : NaN
        (ϱ[i,t,s] - net[:ρ̅][i]) >= 0 ? ϱ_cvar[t,s] += (ϱ[i,t,s] - net[:ρ̅][i]) : NaN
        (net[:ρ̲][i] - ϱ[i,t,s]) >= 0 ? ϱ_cvar[t,s] += (net[:ρ̲][i] - ϱ[i,t,s]) : NaN
    end
    for t in 1:set[:T], i in net[:E_a], s in 1:set[:Ns]
        ϕ[i,t,s] <= 0  ? ϕ_mag[t] += ϕ[i,t,s] : NaN
        ϕ[i,t,s] <= 0  ? ϕ_cvar[t,s] += ϕ[i,t,s] : NaN

        (κ[i,t,s] - net[:κ̅][i]) >= 0  ? κ_mag[t] += (κ[i,t,s] - net[:κ̅][i]) : NaN
        (net[:κ̲][i] - κ[i,t,s]) >= 0  ? κ_mag[t] += (net[:κ̲][i] - κ[i,t,s]) : NaN
        (κ[i,t,s] - net[:κ̅][i]) >= 0  ? κ_cvar[t,s] += (κ[i,t,s] - net[:κ̅][i]) : NaN
        (net[:κ̲][i] - κ[i,t,s]) >= 0  ? κ_cvar[t,s] += (net[:κ̲][i] - κ[i,t,s]) : NaN
    end
    for i in net[:E], s in 1:set[:Ns]
        net[:ψ_o][i] - ψ[i,set[:T],s] >= 0 ? ψ_mag[set[:T]] += (net[:ψ_o][i] - ψ[i,set[:T],s]) : NaN
        net[:ψ_o][i] - ψ[i,set[:T],s] >= 0 ? ψ_cvar[set[:T],s] += (net[:ψ_o][i] - ψ[i,set[:T],s]) : NaN
    end
    ϑ_mag = ϑ_mag./set[:Ns]
    ϱ_mag = ϱ_mag./set[:Ns]/1000
    ϕ_mag = ϱ_mag./set[:Ns]
    κ_mag = ϱ_mag./set[:Ns]/1000
    ψ_mag = ψ_mag./set[:Ns]

    ϑ_cvar = mean(sort(sum(ϑ_cvar,dims=1),dims=2)[950:end])
    ϱ_cvar = mean(sort(sum(ϱ_cvar,dims=1),dims=2)[950:end])/1000
    ψ_cvar = mean(sort(sum(ψ_cvar,dims=1),dims=2)[950:end])
    ϕ_cvar = mean(sort(sum(ϕ_cvar,dims=1),dims=2)[950:end])
    κ_cvar = mean(sort(sum(κ_cvar,dims=1),dims=2)[950:end])/1000

    ### Expected cost
    cost_exp = [net[:c1]'*solution[:Θ][t]*S(t)*ξ[:μ̂] + tr(solution[:Θ][t]'*diagm(net[:c2])*solution[:Θ][t]*S(t)*ξ[:Σ̂]*S(t)') for t in 1:set[:T]]
    mean_lp_deployment = sum(sum(ψ[:,t,s]) for t in 1:set[:T] for s in 1:set[:Ns])*1/set[:Ns]
    mean_compr_deployment = sum(sum(κ[i,t,s]) for i in setdiff(net[:E_a],[51 50]) for t in 1:set[:T] for s in 1:set[:Ns])*1/set[:Ns]
    mean_valve_deployment = sum(sum(abs.(κ[i,t,s])) for i in [51 50] for t in 1:set[:T] for s in 1:set[:Ns])*1/set[:Ns]
    return Dict(:mean_compr_deployment => mean_compr_deployment,
                :mean_lp_deployment => mean_lp_deployment,
                :mean_valve_deployment => mean_valve_deployment,
                :cost_exp => sum(cost_exp),
                :press_violation_exp => sum(ϱ_mag .+ κ_mag),
                :mass_violation_exp => sum(ϑ_mag .+ ϕ_mag .+ ψ_mag),
                :press_violation_wc => sum(ϱ_cvar .+ κ_cvar),
                :mass_violation_wc => sum(ϑ_cvar .+ ϕ_cvar .+ ψ_cvar),
                :first_stage_inj => sum(solution[:Θ][1]*S(1)*ξ[:μ̂]))
end

function cc_opt_MISOCP(net,set,top_data)
    C = Int(2^length(set[:b_valve_edge]))
    # build model
    model = Model(optimizer_with_attributes(Mosek.Optimizer, "LOG" => 0))
    # variable declaration
    P =  [@variable(model, [1:length(net[:N]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    F =  [[@variable(model, [1:length(net[:E]), 1:set[:k_t][t]]) for c = 1:C] for t in 1:set[:T]]
    Z =  [[@variable(model, [1:length(net[:E]), 1:set[:k_t][t]]) for c = 1:C] for t in 1:set[:T]]
    v =  [@variable(model, [1:1], Bin) for c in 1:C]
    Φ =  [@variable(model, [1:length(net[:E]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    K =  [@variable(model, [1:length(net[:E]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    Θ =  [@variable(model, [1:length(net[:N]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    Φ⁻ = [@variable(model, [1:length(net[:E]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    Φ⁺ = [@variable(model, [1:length(net[:E]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    Ψ =  [@variable(model, [1:length(net[:E]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    s_ϱ =  [@variable(model, [1:length(net[:N])]) for t = 1:set[:T]]
    # minimize gas production cost
    @objective(model, Min, sum(net[:c1]' * Θ[t] * S(t) * ξ[:μ̂]
                                + tr(Θ[t]' * diagm(net[:c2]) * Θ[t] * S(t) * ξ[:Σ̂] * S(t)')
                           for t in 1:set[:T])
                           + set[:α_ϱ] * sum(sum(s_ϱ[t]) for t in 2:set[:T])
            )
    # variance constraints
    @constraint(model, pre_var[t=2:set[:T], i=net[:N]],
        [s_ϱ[t][i];(cholesky(ξ[:Σ̂]).U .+ zeros(set[:k_t][set[:T]],set[:k_t][set[:T]]))*(P[t]*S(t) - P[t-1]*S(t-1))[i,:]] in SecondOrderCone())
    # equivalents of stochastic equations
    @constraint(model, gas_bal[t=1:set[:T]], (net[:A⁺]*Φ⁺[t]*S(t) + net[:A⁻]*Φ⁻[t]*S(t) - Θ[t]*S(t) + net[:B]*K[t]*S(t) + Δ(t)*S(t)) .== 0)
    @constraint(model, lin_wey_1[t=1:set[:T],c=1:C], (F[t][c]*S(t) - [top_data[c][:W_0][t] zeros(length(net[:E]),set[:k_t][t]-1)]*S(t) - top_data[c][:W_1][t]*P[t]*S(t) - top_data[c][:W_2][t]*K[t]*S(t)) .== 0)
    @constraint(model, lin_wey_2[t=1:set[:T]], (Φ[t] .- sum(Z[t][c] for c in 1:C))*S(t) .== 0)
    @constraint(model, bilin_1[t=1:set[:T],c=1:C,i=net[:E],j=1:t], set[:M̲] <= Z[t][c][i,j] <= set[:M̅])
    @constraint(model, bilin_2[t=1:set[:T],c=1:C,i=1,j=1:t], set[:M̲]*v[c] .<= Z[t][c][i,j])
    @constraint(model, bilin_3[t=1:set[:T],c=1:C,i=1,j=1:t], set[:M̅]*v[c] .>= Z[t][c][i,j])
    @constraint(model, bilin_4[t=1:set[:T],c=1:C,i=1,j=1:t], F[t][c][i,j] .- (1 .- v[c]) * set[:M̅] .<= Z[t][c][i,j])
    @constraint(model, bilin_5[t=1:set[:T],c=1:C,i=1,j=1:t], F[t][c][i,j] .- (1 .- v[c]) * set[:M̲] .>= Z[t][c][i,j])
    @constraint(model, bin_con, sum(v) .== 1)
    @constraint(model, ref_pre[t=1:set[:T]], (P[t][net[:ref],:]' - sum(v[c].*[top_data[c][:point][:ϱ][t][net[:ref]] zeros(1,set[:k_t][t]-1)] for c in 1:C))*S(t) .== 0)
    @constraint(model, av_flow[t=1:set[:T]], (Φ[t]*S(t) - 1/2*Φ⁺[t]*S(t) - 1/2*Φ⁻[t]*S(t)) .== 0)
    @constraint(model, l_p_def[t=1:set[:T]], (Ψ[t]*S(t) - 1/2*diagm(net[:s])*(K[t] + abs.(net[:A])'*P[t])*S(t)) .== 0)
    @constraint(model, l_p_bal_1[t=1], (Ψ[t]*S(t) - Φ⁺[t]*S(t) + Φ⁻[t]*S(t)) .== 0)
    @constraint(model, l_p_bal_t[t=2:set[:T]], (Ψ[t]*S(t) - Ψ[t-1]*S(t-1) - Φ⁺[t]*S(t) + Φ⁻[t]*S(t)) .== 0)
    # reformulated chance constraints
    @constraint(model, ϑ_max[t=1:set[:T], i=1:length(net[:N])],
        [(net[:ϑ̅] - Θ[t]*S(t)*ξ[:μ̂])[i]; z(set[:ε̅_g]*100)*ξ[:Σ̅]*(Θ[t]*S(t))[i,:]]
            in SecondOrderCone())
    @constraint(model, ϑ_min[t=1:set[:T], i=1:length(net[:N])],
        [(Θ[t]*S(t)*ξ[:μ̂] .- net[:ϑ̲])[i]; z(set[:ε̅_g]*100)*ξ[:Σ̅]*(Θ[t]*S(t))[i,:]]
            in SecondOrderCone())
    @constraint(model, ρ_max[t=1:set[:T], i=1:length(net[:N])],
        [(net[:ρ̅] - P[t]*S(t)*ξ[:μ̂])[i]; z(set[:ε̅_g]*100)*ξ[:Σ̅]*(P[t]*S(t))[i,:]]
            in SecondOrderCone())
    @constraint(model, ρ_min[t=1:set[:T], i=1:length(net[:N])],
        [(P[t]*S(t)*ξ[:μ̂] - net[:ρ̲])[i]; z(set[:ε̅_g]*100)*ξ[:Σ̅]*(P[t]*S(t))[i,:]]
            in SecondOrderCone())
    @constraint(model, κ_max[t=1:set[:T], i=1:length(net[:E])],
        [(net[:κ̅] - K[t]*S(t)*ξ[:μ̂])[i]; z(set[:ε̅_g]*100)*ξ[:Σ̅]*(K[t]*S(t))[i,:]]
            in SecondOrderCone())
    @constraint(model, κ_min[t=1:set[:T], i=1:length(net[:E])],
        [(K[t]*S(t)*ξ[:μ̂] - net[:κ̲])[i]; z(set[:ε̅_g]*100)*ξ[:Σ̅]*(K[t]*S(t))[i,:]]
            in SecondOrderCone())
    @constraint(model, ψ_T[t=set[:T], i=1:length(net[:E])],
        [(Ψ[t]*S(t)*ξ[:μ̂] - net[:ψ_o])[i]; z(set[:ε̅_g]*100)*ξ[:Σ̅]*(Ψ[t]*S(t))[i,:]]
            in SecondOrderCone())
    @constraint(model, φ_pl[t=1:set[:T], i=net[:E_a]],
        [(Φ[t]*S(t)*ξ[:μ̂] .- 0)[i]; z(set[:ε̅_g]*100)*ξ[:Σ̅]*(Φ[t]*S(t))[i,:]]
            in SecondOrderCone())
    # solve model
    optimize!(model)
    @info("done solving topology optimization in $(solve_time(model)) sec. with status: $(termination_status(model))")
    # return solution
    solution = Dict(
    :v => [JuMP.value.(v[c]) for c in 1:C],
    :P => [JuMP.value.(P[t]) for t in 1:set[:T]],
    :Φ => [JuMP.value.(Φ[t]) for t in 1:set[:T]],
    :K => [JuMP.value.(K[t]) for t in 1:set[:T]],
    :Θ => [JuMP.value.(Θ[t]) for t in 1:set[:T]],
    :Φ⁻ => [JuMP.value.(Φ⁻[t]) for t in 1:set[:T]],
    :Φ⁺ => [JuMP.value.(Φ⁺[t]) for t in 1:set[:T]],
    :Ψ => [JuMP.value.(Ψ[t]) for t in 1:set[:T]],
    :cost => sum(net[:c1]' * JuMP.value.(Θ[t]) * S(t) * ξ[:μ̂] + tr(JuMP.value.(Θ[t])' * diagm(net[:c2]) * JuMP.value.(Θ[t]) * S(t) * ξ[:Σ̂] * S(t)') for t in 1:set[:T]),
    :ϱ_std_tot => sum(sum(JuMP.value.(α_ϱ[t])) for t in 2:1:set[:T]) ,
    :status => termination_status(model),
    :obj => JuMP.objective_value(model),
    :model => model)
    return solution
end



function cc_opt_MISOCP_2(net,set,top_data)
    C = Int(2^length(set[:b_valve_edge]))
    # build model
    model = Model(optimizer_with_attributes(Mosek.Optimizer, "MSK_DPAR_MIO_TOL_REL_GAP" => 1e-5))
    # model = Model(optimizer_with_attributes(Gurobi.Optimizer))
    # set_optimizer_attribute(model, "MIQCPMethod", 1)
    # variable declaration
    P =  [@variable(model, [1:length(net[:N]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    F =  [[@variable(model, [1:length(net[:E]), 1:set[:k_t][t]]) for c = 1:C] for t in 1:set[:T]]
    Z =  [[@variable(model, [1:length(net[:E]), 1:set[:k_t][t]]) for c = 1:C] for t in 1:set[:T]]
    v =  [@variable(model, [1:1], Bin) for c in 1:C]
    Φ =  [@variable(model, [1:length(net[:E]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    K =  [@variable(model, [1:length(net[:E]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    Θ =  [@variable(model, [1:length(net[:N]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    Φ⁻ = [@variable(model, [1:length(net[:E]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    Φ⁺ = [@variable(model, [1:length(net[:E]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    Ψ =  [@variable(model, [1:length(net[:E]), 1:set[:k_t][t]]) for t = 1:set[:T]]
    α_ϱ =  [@variable(model, [1:length(net[:N])]) for t = 1:set[:T]]
    # minimize gas production cost
    @objective(model, Min, sum(net[:c1]' * Θ[t] * S(t) * ξ[:μ̂]
                                + tr(Θ[t]' * diagm(net[:c2]) * Θ[t] * S(t) * ξ[:Σ̂] * S(t)')
                           for t in 1:set[:T])
                           + set[:α_ϱ] * sum(sum(α_ϱ[t]) for t in 2:set[:T])
            )
    # variance constraints
    @constraint(model, pre_var[t=2:set[:T], i=net[:N]],
        [α_ϱ[t][i];(cholesky(ξ[:Σ̂]).U .+ zeros(set[:k_t][set[:T]],set[:k_t][set[:T]]))*(P[t]*S(t) - P[t-1]*S(t-1))[i,:]] in SecondOrderCone())
    # equivalents of stochastic equations
    @constraint(model, gas_bal[t=1:set[:T]],
        (net[:A⁺]*Φ⁺[t]*S(t) + net[:A⁻]*Φ⁻[t]*S(t) - Θ[t]*S(t) + net[:B]*K[t]*S(t) + Δ(t)*S(t)) .== 0)
    @constraint(model, lin_wey[t=1:set[:T],c=1:C],
        (F[t][c]*S(t) - [top_data[c][:W_0][t] zeros(length(net[:E]),set[:k_t][t]-1)]*S(t) - top_data[c][:W_1][t]*P[t]*S(t) - top_data[c][:W_2][t]*K[t]*S(t)) .== 0)
    @constraint(model, lin_wey_2[t=1:set[:T]],
        (Φ[t] .- sum(Z[t][c] for c in 1:C))*S(t) .== 0)
    @constraint(model, lin_wey_3[t=1:set[:T]],
        (P[t][net[:ref],:]' -  sum(v[c].*[top_data[c][:point][:ϱ][t][net[:ref]] zeros(1,set[:k_t][t]-1)] for c in 1:C)) * S(t) .== 0)
    @constraint(model, av_flow[t=1:set[:T]],
        (Φ[t]*S(t) - 1/2*Φ⁺[t]*S(t) - 1/2*Φ⁻[t]*S(t)) .== 0)
    @constraint(model, l_p_def[t=1:set[:T]],
        (Ψ[t]*S(t) - 1/2*diagm(net[:s])*(K[t] + abs.(net[:A])'*P[t])*S(t)) .== 0)
    @constraint(model, l_p_bal_1[t=1],
        (Ψ[t]*S(t) - Φ⁺[t]*S(t) + Φ⁻[t]*S(t)) .== 0)
    @constraint(model, l_p_bal_t[t=2:set[:T]],
        (Ψ[t]*S(t) - Ψ[t-1]*S(t-1) - Φ⁺[t]*S(t) + Φ⁻[t]*S(t)) .== 0)
    # linearization of the bilinear term
    @constraint(model, bilin_1[t=1:set[:T],c=1:C,i=net[:E],j=1:t],
        Z[t][c][i,j] <= set[:M̅])
    @constraint(model, bilin_2[t=1:set[:T],c=1:C,i=net[:E],j=1:t],
        Z[t][c][i,j] >= set[:M̲])
    @constraint(model, bilin_3[t=1:set[:T],c=1:C,i=net[:E],j=1:t],
        Z[t][c][i,j] .<= set[:M̅] * v[c])
    @constraint(model, bilin_4[t=1:set[:T],c=1:C,i=net[:E],j=1:t],
        Z[t][c][i,j] .>= set[:M̲] * v[c])
    @constraint(model, bilin_5[t=1:set[:T],c=1:C,i=net[:E],j=1:t],
        Z[t][c][i,j] .<=  F[t][c][i,j] .- (1 .- v[c]) * set[:M̲])
    @constraint(model, bilin_6[t=1:set[:T],c=1:C,i=net[:E],j=1:t],
        Z[t][c][i,j] .>=  F[t][c][i,j] .- (1 .- v[c]) * set[:M̅])
    # sum up all binaries to one
    @constraint(model, sum(v) .== 1)
    # reformulated chance constraints
    @constraint(model, ϑ_max[t=1:set[:T], i=1:length(net[:N])],
        [(net[:ϑ̅] - Θ[t]*S(t)*ξ[:μ̂])[i]; z(set[:ε̅_g])*ξ[:Σ̅]*(Θ[t]*S(t))[i,:]]
            in SecondOrderCone())
    @constraint(model, ϑ_min[t=1:set[:T], i=1:length(net[:N])],
        [(Θ[t]*S(t)*ξ[:μ̂] .- net[:ϑ̲])[i]; z(set[:ε̅_g])*ξ[:Σ̅]*(Θ[t]*S(t))[i,:]]
            in SecondOrderCone())
    @constraint(model, ρ_max[t=1:set[:T], i=1:length(net[:N])],
        [(net[:ρ̅] - P[t]*S(t)*ξ[:μ̂])[i]; z(set[:ε̅_g])*ξ[:Σ̅]*(P[t]*S(t))[i,:]]
            in SecondOrderCone())
    @constraint(model, ρ_min[t=1:set[:T], i=1:length(net[:N])],
        [(P[t]*S(t)*ξ[:μ̂] - net[:ρ̲])[i]; z(set[:ε̅_g])*ξ[:Σ̅]*(P[t]*S(t))[i,:]]
            in SecondOrderCone())
    @constraint(model, κ_max[t=1:set[:T], i=1:length(net[:E])],
        [(net[:κ̅] - K[t]*S(t)*ξ[:μ̂])[i]; z(set[:ε̅_g])*ξ[:Σ̅]*(K[t]*S(t))[i,:]]
            in SecondOrderCone())
    @constraint(model, κ_min[t=1:set[:T], i=1:length(net[:E])],
        [(K[t]*S(t)*ξ[:μ̂] - net[:κ̲])[i]; z(set[:ε̅_g])*ξ[:Σ̅]*(K[t]*S(t))[i,:]]
            in SecondOrderCone())
    @constraint(model, ψ_T[t=set[:T], i=1:length(net[:E])],
        [(Ψ[t]*S(t)*ξ[:μ̂] - net[:ψ_o])[i]; z(set[:ε̅_g])*ξ[:Σ̅]*(Ψ[t]*S(t))[i,:]]
            in SecondOrderCone())
    @constraint(model, φ_pl[t=1:set[:T], i=net[:E_a]],
        [(Φ[t]*S(t)*ξ[:μ̂] .- 0)[i]; z(set[:ε̅_g])*ξ[:Σ̅]*(Φ[t]*S(t))[i,:]]
            in SecondOrderCone())
    # solve model
    optimize!(model)
    @info("cc model terminates with status: $(termination_status(model))")
    # return solution
    solution = Dict(
    :v => [JuMP.value.(v[c]) for c in 1:C],
    :P => [JuMP.value.(P[t]) for t in 1:set[:T]],
    :Φ => [JuMP.value.(Φ[t]) for t in 1:set[:T]],
    :K => [JuMP.value.(K[t]) for t in 1:set[:T]],
    :Θ => [JuMP.value.(Θ[t]) for t in 1:set[:T]],
    :Φ⁻ => [JuMP.value.(Φ⁻[t]) for t in 1:set[:T]],
    :Φ⁺ => [JuMP.value.(Φ⁺[t]) for t in 1:set[:T]],
    :Ψ => [JuMP.value.(Ψ[t]) for t in 1:set[:T]],
    :cost => sum(net[:c1]' * JuMP.value.(Θ[t]) * S(t) * ξ[:μ̂] + tr(JuMP.value.(Θ[t])' * diagm(net[:c2]) * JuMP.value.(Θ[t]) * S(t) * ξ[:Σ̂] * S(t)') for t in 1:set[:T]),
    :ϱ_std_tot => sum(sum(JuMP.value.(α_ϱ[t])) for t in 2:1:set[:T]) ,
    :status => termination_status(model),
    :obj => JuMP.objective_value(model),
    :model => model)
    return solution
end
