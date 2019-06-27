function MCPCA_PopGen(DS, kyfan_q; discretize_method="interval", num_iter=10,num_init=0)
	## input data matrix
	n, p = size(DS)

	## clean the data
	DS = DS[:,[length(unique(DS[:,i])) for i in 1:p] .> 1];
	DS = DS[:,[iqr(DS[:,i]) for i in 1:p] .!= 0]
	n, p = size(DS)

	## Discretize
	if discretize_method != "none" 
		DS_dis = Discretize(DS, dis_method=discretize_method);
	elseif discretize_method == "none"
		DS_dis = DS;
	end

	phimat, funcell = MCPCA_sample_disc_wrapper(DS_dis, kyfan_q, num_iter, num_init)

	return phimat

end
