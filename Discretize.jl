######################
## jenks discretizer 
######################
function Jenks_label(inputvec, num_breaks)
    R"breaks=getJenksBreaks($inputvec, k=$num_breaks)"
    breaks = rcopy(R"breaks")
    label = zeros(length(inputvec))
    for i=1:(length(breaks)-1)
        id = [(inputvec[j] >= breaks[i] && inputvec[j] <= breaks[i+1]) for j in 1:length(inputvec)]
        label[id] .= i
    end
    return label
end


## discretize
function Discretize(inputdata; dis_method="interval")
    (n, p) = size(inputdata);
    data_dis = Array{Int64,2}(undef, n, p);        
    
    if dis_method == "interval"
        for i=1:p
            num_bins=min( max(Int(ceil(n^(1/3)/iqr(inputdata[:, i]))),2), n-1);
            dis_map=LinearDiscretizer(binedges(DiscretizeUniformWidth(num_bins), inputdata[:, i]))
            data_dis[:,i] = encode(dis_map, inputdata[:,i])
        end    
    elseif dis_method == "counts"
        inputdata = inputdata + rand(Uniform(-1,1), n, p) .* 1e-13;
        for i=1:p
            um_bins=min( max(Int(ceil(n^(1/3)/iqr(inputdata[:, i]))),2), n-1);
            dis_map=LinearDiscretizer(binedges(DiscretizeUniformCount(num_bins), inputdata[:, i]))
            data_dis[:,i] = encode(dis_map, inputdata[:,i])
        end
    elseif dis_method == "Jenks"
    	R"""
            dyn.load('./jenksBrks.so')
            source('./getJenksBreaks.R')
        """
        for i=1:p
            inputvec=[round(inputdata[j,i] , digits=7) for j in 1:n];
            num_bins=min( max(Int(ceil(n^(1/3)/iqr(inputvec))),2), n-1 );
            data_dis[:,i] = Jenks_label(inputvec, num_bins+1)
        end 
    end
    return data_dis
end
