function MCPCA_sample_disc_wrapper(X_input, q, num_iter, num_init)

 ######################################
 ## Labeling alphabets with {1,2,...,B_i}
 
 n,p = size(X_input);
 X0 = zeros(n, p);
 map_alphabet0 = Array{Any, 2}(undef, (1, p));
 for i = 1:p
     alpha0 = sort(unique(X_input[:,i]))
     alpha1 = collect(1:length(alpha0))
     map_alphabet0[1,i] = hcat(alpha0,alpha1)
     for j = 1:length(alpha0)
         ind_t = collect(1:n)[X_input[:,i] .== alpha0[j]]
         X0[ind_t,i] .= j
     end
 end

 ######################################
 ## Applying MCPCA with identity initialization

 println("MCPCA: iteration 0")
 obj_vec_test, phi_mat_test = MCPCA_sample_disc(X0, q, num_iter, 0)

 ######################################
 ## extracting functions()
 fun_cell = Array{Any, 2}(undef, (p, 1));
 for ii = 1:p
     alpha0 = sort(unique(X_input[:,ii]))
     alpha1 = alpha0*0.0
     for jj = 1:length(alpha0)
         ind_t = collect(1:n)[X_input[:,ii] .== alpha0[jj]]
         if length(unique(phi_mat_test[ind_t,ii])) == 1
             alpha1[jj] = unique(phi_mat_test[ind_t,ii])[1]
         end
     end
     fun_cell[ii,1] = hcat(alpha0,alpha1)
 end

 ######################################
 ## saving results for init 0
 total_res = Array{Any, 2}(undef, (3,num_init+1));
 kyfan_vec = zeros(1, num_init+1)
 kyfan_vec[1] = obj_vec_test[end];
 total_res[1,1] = obj_vec_test;
 total_res[2,1] = phi_mat_test;
 total_res[3,1] = fun_cell;
 ######################################
 ## random initialization of MCPCA

 for i = 1:num_init
     println("MCPCA: iteration $(i)")
     obj_vec_test,phi_mat_test = MCPCA_sample_disc(X0, q, num_iter, 1)

     ## extracting functions()
     fun_cell = Array{Any, 2}(undef, (p, 1));
     for ii = 1:p
         alpha0 = sort(unique(X_input[:,ii]))
         alpha1 = alpha0*0.0
         for jj = 1:length(alpha0)
             ind_t = collect(1:n)[X_input[:,ii] .== alpha0[jj]]
             if length(unique(phi_mat_test[ind_t,ii])) == 1
                 alpha1[jj] = unique(phi_mat_test[ind_t,ii])[1]
             end
         end
         fun_cell[ii,1] = hcat(alpha0,alpha1)
     end

     ## saving results for init i
     kyfan_vec[i+1] = obj_vec_test[end]
     total_res[1,i+1] = obj_vec_test
     total_res[2,i+1] = phi_mat_test
     total_res[3,i+1] = fun_cell
 end

 ######################################
 ## selecting best initilization
 ind_m = argmax(kyfan_vec[:])
 phi_mat = total_res[2,ind_m];
 fun_cell = total_res[3,ind_m];

 ######################################
 ## making positive orientation
 n,p = size(X_input)
 for i = 1:p
     temp = cor(X_input[:,i],phi_mat[:,i])
     if temp < 0
         phi_mat[:,i] = -phi_mat[:,i]
         fun_temp = fun_cell[i,1]
         fun_temp[:,2] = -fun_temp[:,2]
         fun_cell[i,1] = fun_temp
     end
 end
 return phi_mat,fun_cell
end
