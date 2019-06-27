function MCPCA_sample_disc(X, q, num_iter, control_init)
    ## initialization
    n,p = size(X);

    if control_init == 0
        phi_mat = copy(X);
    else
        phi_mat = X .* 0.0;
        for i = 1:p
            alphabet = sort(unique(X[:,i]));
            alphabet2 = randperm(length(alphabet));
            for j = 1:length(alphabet)
                ind_t = collect(1:n)[X[:,i] .== alphabet[j]];
                phi_mat[ind_t,i] .= alphabet2[j];
            end
        end
    end

    phi_mat = normalize_matrix(phi_mat);
    D, V= eigenq(phi_mat, q);

    obj_vec = zeros(1, num_iter+1);
    obj_vec[1,1] = sum(D);
    fun_cell = Array{Any, 2}(undef, (p, 1));

    for iter=1:num_iter
        ## updating phi_i
        for i=1:p
            phi_mat,fun_cell = update_phi_disc(i, X, phi_mat, V, fun_cell);
        end
        ## updating V
        D, V= eigenq(phi_mat, q);
        obj_vec[1,iter+1] = sum(D);
    end
    return obj_vec, phi_mat
end
