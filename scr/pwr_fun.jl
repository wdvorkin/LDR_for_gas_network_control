# power functions
function load_pwr_data(set)
    data_net = PowerModels.parse_file("data/$(set[:pwr_case])")
    gen_data = CSV.read("data/gen_data.csv", DataFrame; header=1)
    load_dyn = CSV.read("data/load_factor.csv", DataFrame; header=1)[:,1]
    # Network size
    G = size(gen_data,1)
    N = length(data_net["bus"])
    E = length(data_net["branch"])
    D = length(data_net["load"])

    # order bus indexing
    bus_keys=collect(keys(data_net["bus"]))
    bus_key_dict = Dict()
    for i in 1:N
        push!(bus_key_dict, i => bus_keys[i])
    end
    node(key) = [k for (k,v) in bus_key_dict if v == key][1]

    # Load generation data
    p̅ = zeros(G); p̲ = zeros(G); c1 = zeros(G); c2 = zeros(G); M_p = zeros(N,G)
    for g in 1:G
        p̅[g] = gen_data[g,:P_max]
        p̲[g] = gen_data[g,:P_min] * 0
        c1[g] = gen_data[g,:C]
        c2[g] = gen_data[g,:C_2]
        M_p[gen_data[g,:node],g] = 1
    end

    # Load demand data
    load_key=collect(keys(data_net["load"]))
    d = zeros(D); M_d = zeros(N,D)
    for h in load_key
        d[parse(Int64,h)] = data_net["load"][h]["pd"]*data_net["baseMVA"] + 1e-3
        M_d[node(string(data_net["load"][h]["load_bus"])),parse(Int64,h)] = 1
    end
    load_dyn = ones(set[:T])

    # Load transmission data
    line_key=collect(keys(data_net["branch"]))
    β = zeros(E); f̅ = zeros(E); n_s = trunc.(Int64,zeros(E)); n_r = trunc.(Int64,zeros(E))
    for l in line_key
        β[data_net["branch"][l]["index"]] = -imag(1/(data_net["branch"][l]["br_r"] + data_net["branch"][l]["br_x"]im))
        n_s[data_net["branch"][l]["index"]] = data_net["branch"][l]["f_bus"]
        n_r[data_net["branch"][l]["index"]] = data_net["branch"][l]["t_bus"]
        f̅[data_net["branch"][l]["index"]] = data_net["branch"][l]["rate_a"]*data_net["baseMVA"]
    end
    # merge parallel lines
    ff = zeros(N,N); ββ = zeros(N,N)
    for l in line_key
        ff[node(string(n_s[data_net["branch"][l]["index"]])),node(string(n_r[data_net["branch"][l]["index"]]))] += f̅[data_net["branch"][l]["index"]]
        ff[node(string(n_r[data_net["branch"][l]["index"]])),node(string(n_s[data_net["branch"][l]["index"]]))] += f̅[data_net["branch"][l]["index"]]
        ββ[node(string(n_s[data_net["branch"][l]["index"]])),node(string(n_r[data_net["branch"][l]["index"]]))]  = β[data_net["branch"][l]["index"]]
        ββ[node(string(n_r[data_net["branch"][l]["index"]])),node(string(n_s[data_net["branch"][l]["index"]]))]  = β[data_net["branch"][l]["index"]]
    end
    # find all parallel lines
    parallel_lines = []
    for l in line_key, e in line_key
        if l != e && node(string(n_s[data_net["branch"][l]["index"]])) == node(string(n_s[data_net["branch"][e]["index"]])) && node(string(n_r[data_net["branch"][l]["index"]])) == node(string(n_r[data_net["branch"][e]["index"]]))
            push!(parallel_lines,l)
        end
    end
    # update number of edges
    E = E - Int(length(parallel_lines)/2)
    # get s and r ends of all edge
    n_s = trunc.(Int64,zeros(E)); n_r = trunc.(Int64,zeros(E))
    ff = LowerTriangular(ff)
    for l in 1:E
        n_s[l] = findall(!iszero, ff)[l][1]
        n_r[l] = findall(!iszero, ff)[l][2]
    end
    β = zeros(E); f̅ = zeros(E);
    likely_congested_lines = [2 25 34 76 77 112 116 174]
    for l in 1:E
        β[l] = ββ[n_s[l],n_r[l]]
        l in likely_congested_lines ? f̅[l] = ff[n_s[l],n_r[l]]*3 : f̅[l] = ff[n_s[l],n_r[l]]*3
    end

    # Find reference node
    ref = 0
    for n in 1:N
        if sum(M_p[n,:]) == 0 &&  sum(M_d[n,:]) == 0 == 0
            ref = n
        end
    end

    # Compute PTDF matrix
    B_line = zeros(E,N); B̃_bus = zeros(N,N); B = zeros(N,N)
    for n in 1:N
        for l in 1:E
            if n_s[l] == n
                B[n,n] += β[l]
                B_line[l,n] = β[l]
            end
            if n_r[l] == n
                B[n,n] += β[l]
                B_line[l,n] = -β[l]
            end
        end
    end
    for l in 1:E
        B[Int(n_s[l]),Int(n_r[l])] = - β[l]
        B[Int(n_r[l]),Int(n_s[l])] = - β[l]
    end
    B̃_bus = remove_col_and_row(B,ref)
    B̃_bus = inv(B̃_bus)
    B̃_bus = build_B̆(B̃_bus,ref)
    PTDF = B_line*B̃_bus

    # power-gas data
    gas_bus = [10 12 15 25 26 61 65 66 69]
    power_node = [45 33 36 31 30 40 39 38 25]
    Λ = zeros(48,N)
    Λ = zeros(48,N)
    for i in 1:length(power_node)
        Λ[power_node[i],gas_bus[i]] = 0.25
    end

    @info("done loading power system data")

    # safe network data
    net = Dict(
    # transmission data
    :f̅ => f̅, :n_s => n_s, :n_r => n_r, :T => round.(PTDF,digits=8),
    # load data
    :d => round.(d,digits=5), :M_d => M_d, :k_l => load_dyn,
    # generation data
    :p̅ => round.(p̅,digits=5), :p̲ => round.(p̲,digits=5), :M_p => M_p,
    :c1 => round.(c1,digits=5), :c2 => round.(c2,digits=5), :c̀2 => round.(sqrt.(c2),digits=5),
    # graph data
    :N => N, :E => E, :G => G, :D => D, :ref => ref,
    :gas_bus => [10 12 15 25 26 61 65 66 69],
    :Λ => Λ
    )
    return net
end
function solve_chance_constrained_OPF(net,set,ξ)
    # DC-OPF definition
    model = Model(optimizer_with_attributes(Mosek.Optimizer, "LOG" => 0))
    # model variables
    P =  [@variable(model, [1:net[:G], 1:set[:k_t][t]]) for t = 1:set[:T]]
    v =  [@variable(model, [1:net[:N]]) for t = 1:set[:T]]
    # auxiliary variables
    A =  [@variable(model, [1:net[:N], 1:set[:k_t][t]]) for t = 1:set[:T]]
    B =  [@variable(model, [1:net[:E], 1:set[:k_t][t]]) for t = 1:set[:T]]
    y_f =  [@variable(model, [1:net[:E]]) for t = 1:set[:T]]
    x_f =  [@variable(model, [1:net[:E]]) for t = 1:set[:T]]
    y_p =  [@variable(model, [1:net[:G]]) for t = 1:set[:T]]
    x_p =  [@variable(model, [1:net[:G]]) for t = 1:set[:T]]
    # model objective (exp cost + var of non-CCGT units)
    @objective(model, Min, sum(net[:c1]'P[t]*S(t)*ξ[:μ̂] + 10000 * sum(v[t][i] for i in setdiff(1:net[:N],net[:gas_bus])) for t in 1:set[:T]))
    # prioritize CCGT plants to balance renewables
    @constraint(model, ϕ[t=1:set[:T],i=setdiff(1:net[:N],net[:gas_bus])],
        [v[t][i];ξ[:Σ̅]*(net[:M_p]*P[t]*S(t))[i,:]] in SecondOrderCone())
    # OPF equations
    @constraint(model, α[t=1:set[:T]], A[t] .== (net[:M_p]*P[t] .+ ξ[:M_w]*W(t) .- net[:M_d]*D(t)*net[:k_l][t]))
    @constraint(model, β[t=1:set[:T]], B[t] .== net[:T] * A[t])
    @constraint(model, λ[t=1:set[:T]], ones(net[:N])' * A[t] .== 0)
    # double-sided reformulation of power flow chance constraints
    @constraint(model, μ_1[t=1:set[:T], l=1:net[:E]], [sqrt(set[:ε̅_p]) * (net[:f̅] .- x_f[t])[l]; vcat(ξ[:Σ̅]*(B[t]*S(t))[l,:],y_f[t][l])] in SecondOrderCone())
    @constraint(model, μ_2[t=1:set[:T], l=1:net[:E]], (B[t] * S(t) * ξ[:μ̂])[l] <= y_f[t][l] + x_f[t][l])
    @constraint(model, μ_3[t=1:set[:T], l=1:net[:E]], (B[t] * S(t) * ξ[:μ̂])[l] >= - y_f[t][l] - x_f[t][l])
    @constraint(model, μ_4[t=1:set[:T], l=1:net[:E]], net[:f̅][l] >= x_f[t][l] >= 0)
    @constraint(model, μ_5[t=1:set[:T], l=1:net[:E]], y_f[t][l] >= 0)
    # double-sided generation (reformulated) chance constraints
    @constraint(model, κ_1[t=1:set[:T], i=1:net[:G]], [sqrt(set[:ε̅_p]) * (1/2*(net[:p̅] - net[:p̲]) - x_p[t])[i]; vcat(ξ[:Σ̅]*(P[t]*S(t))[i,:],y_p[t][i])] in SecondOrderCone())
    @constraint(model, κ_2[t=1:set[:T], i=1:net[:G]], (P[t]*S(t)*ξ[:μ̂])[i] - 1/2*(net[:p̅] - net[:p̲])[i] <= y_p[t][i] + x_p[t][i])
    @constraint(model, κ_3[t=1:set[:T], i=1:net[:G]], (P[t]*S(t)*ξ[:μ̂])[i] - 1/2*(net[:p̅] - net[:p̲])[i] >= - y_p[t][i] - x_p[t][i])
    @constraint(model, κ_4[t=1:set[:T], i=1:net[:G]], 1/2*(net[:p̅] - net[:p̲])[i] >= x_p[t][i] >= 0)
    @constraint(model, κ_5[t=1:set[:T], i=1:net[:G]], y_p[t][i] >= 0)
    # solve model
    optimize!(model)
    @info("done solving chance-constrained OPF in $(solve_time(model)) sec. with status $(termination_status(model))")
    sol = Dict(:status => termination_status(model),
                :obj => JuMP.objective_value(model),
                :P => [JuMP.value.(P[t]) for t = 1:set[:T]],
                :CPUtime => solve_time(model))
    return sol
end
function wind_forecast(set)
    # power system data
    w_f = CSV.read("data/wind.csv", DataFrame; header=0)[:,1]
    bus = [3 8 11 20 24 26 31 38 43 49 53]
    W_ = length(bus)
    M_w = zeros(net[:N],W_)
    for i in 1:W_
        M_w[bus[i],i] = 1
    end
    # uncertainty data
    k = set[:k_t][set[:T]]
    # mean and covariance
    μ = ones(k-1)
    Σ = zeros(k-1,k-1)
    for i in 1:k-1
        Σ[i,i] = set[:σ]
    end
    # first two moments of altered dimentions
    μ̂ = vcat(1,μ)
    # E[ξξᵀ] = Cov + μμᵀ
    Σ̂ = zeros(k,k)
    Σ̂[1,1] = 1; Σ̂[2:end,1] .= μ; Σ̂[1,2:end] .= μ;
    Σ̂[2:end,2:end] .=  Σ .+ μ*μ'
    # factorization of E[ξξᵀ]
    Σ̅ = cholesky(Σ̂).L .+ zeros(k,k)
    Σ̅[:,1] .= 0
    ξ̂ = ones(k,set[:Ns])
    ξ̂[2:end,:] = rand(MvNormal(μ,Σ),set[:Ns])
    return Dict(:w_f => w_f, :bus => bus, :W => W_, :M_w => M_w, :μ => μ, :Σ => Σ, :μ̂ => μ̂, :Σ̂ => Σ̂, :Σ̅ => Σ̅, :ξ̂ => ξ̂)
end
