# aux functions
ns(l) = Int(net[:n_s][l])
nr(l) = Int(net[:n_r][l])
z(x) = quantile(Normal(0,1),1-x)

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--α_ϱ", "-p"
            help = "pressure penalty α_ϱ"
            arg_type = Float64
            default = 0.0
        "--α_ϑ", "-i"
            help = "injection variability coefficient α_ϑ"
            arg_type = Float64
            default = 0.025
        "--α_ψ", "-l"
            help = "linepack variability coefficient α_ψ"
            arg_type = Float64
            default = 100.0
        "--σ²", "-s"
            help = "wind covariance"
            arg_type = Float64
            default = 0.15
        "--ε̅_p", "-e"
            help = "electric power constraint violation tolerance"
            arg_type = Float64
            default = 0.01
        "--ε̅_g", "-g"
            help = "gas constraint violation tolerance"
            arg_type = Float64
            default = 0.005
        "--Ns", "-n"
            help = "number of samples for out-of-sample analysis"
            arg_type = Int64
            default = 1000
    end
    return parse_args(s)
end
args = parse_commandline()
outdir = "output"
mkpath(outdir)

function remove_col_and_row(B,refbus)
    @assert size(B,1) == size(B,2)
    n = size(B,1)
    return B[1:n .!= refbus, 1:n .!= refbus]
end

function full_matrix(A,ref)
    Nb = size(A,1)+1
    V = zeros(Nb,Nb)
    for i in 1:Nb, j in 1:Nb
        i < ref && j < ref ? V[i,j] = A[i,j] : NaN
        i > ref && j > ref ? V[i,j] = A[i-1,j-1] : NaN
        i > ref && j < ref ? V[i,j] = A[i-1,j] : NaN
        i < ref && j > ref ? V[i,j] = A[i,j-1] : NaN
    end
    return V
end

function build_B̆(B̂inv,refbus)
    Nb = size(B̂inv,1)+1
    B̆ = zeros(Nb,Nb)
    for i in 1:Nb, j in 1:Nb
        if i < refbus && j < refbus
            B̆[i,j] = B̂inv[i,j]
        end
        if i > refbus && j > refbus
            B̆[i,j] = B̂inv[i-1,j-1]
        end
        if i > refbus && j < refbus
            B̆[i,j] = B̂inv[i-1,j]
        end
        if i < refbus && j > refbus
            B̆[i,j] = B̂inv[i,j-1]
        end
    end
    return B̆
end

function W(i)
    W = zeros(length(ξ[:bus]),set[:k_t][i])
    for t in 1:i
        W[:,1] .= ξ[:w_f][1] * set[:w̅]
        if t > 1
            for i in 1:length(ξ[:bus])
                if i <= 4
                    W[i,set[:k_t][t]-2] = (ξ[:w_f][t] - ξ[:w_f][t-1]) * set[:w̅]
                end
                if i > 4 && i <= 9
                    W[i,set[:k_t][t]-1] = (ξ[:w_f][t] - ξ[:w_f][t-1]) * set[:w̅]
                end
                if i > 9
                    W[i,set[:k_t][t]] = (ξ[:w_f][t] - ξ[:w_f][t-1]) * set[:w̅]
                end
            end
        end
    end
    return W
end

function D(i)
    D = zeros(net[:D],set[:k_t][i])
    for t in 1:i
        D[:,1] .= net[:d]
    end
    return D
end

function S(t)
    S = zeros(set[:k_t][t],set[:k_t][set[:T]])
    S[diagind(S)] .= 1
    return S
end

function Δ(i)
    L = zeros(length(net[:N]),set[:k_t][i])
    L[:,1] = net[:δ]

    L_opf = CSV.File("data/gas_uncertainty.csv",header=0) |> Tables.matrix
    # L_opf = [zeros(i)';L_opf].*10
    L_opf = L_opf[:,1:set[:k_t][i]] * 3
    L[findall(x->x>0,L_opf[:,1]),:] = L_opf[findall(x->x>0,L_opf[:,1]),:]
    return L
end
