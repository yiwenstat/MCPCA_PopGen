################################################
## kyfan norm
################################################
function kyfan(X)
	lambda = (svd(X).S).^2;
	return cumsum(lambda)
end

################################################
## eigenq ##
## return first q eigenvalue and eigenvector
## in decreasing order
################################################
function eigenq(phi_mat, q)
	n, p = size(phi_mat)
    Z = phi_mat' * phi_mat ./ n;
    Dz, Vz = eigen(Z); 
    od = sortperm(Dz, rev=true)[1:q]; 
    Dz = Dz[od];
    Vz = Vz[:,od];
    return(Dz, Vz)
end

################################################
## normalize matrix
## columsn of x mean zero and norm sqrt(n)
################################################
function normalize_matrix(x)
     (n,p) = size(x);
     for i = 1:p
         x[:,i] = x[:,i] .- mean(x[:,i]);
         if norm(x[:,i]) != 0
             x[:,i] = x[:,i] ./ norm(x[:,i]) .* sqrt(n);
         else
             println("one of the columns is all zero")
         end
     end
     return x
end

################################################
# updating column index in the phi_mat
# V is a p by q matrix, each column contains an eigenvector
################################################
function update_phi_disc(index, X, phi_mat, V, fun_cell)
    p,q = size(V);
    n = size(phi_mat)[1];

    w = zeros(n,1);
    for r = 1:q, i = 1:p
        if i!=index
            w += V[index,r] * V[i,r] * phi_mat[:,i];
        end
    end

    ################################################
    y_alpha = sort(unique(X[:,index])); # alphabets
    M = length(y_alpha); # number of alphabets

    w_alpha = zeros(M,1);
    phi_new = zeros(n,1);
    for i = 1:M
        ind_t = collect(1:n)[X[:,index] .== y_alpha[i]];
        w_alpha[i] = mean(w[ind_t,1]);
        phi_new[ind_t] .= w_alpha[i];
    end

    ################################################
    phi_new = phi_new .- mean(phi_new);
    if norm(phi_new) !=0
        phi_new = phi_new ./norm(phi_new) .* sqrt(n);
        phi_mat[:,index] = phi_new;
        fun_cell[index,1] = hcat(y_alpha,w_alpha);
    else()
     # disp("norm phi_new is zero, not updating")
    end
    return phi_mat, fun_cell
end

