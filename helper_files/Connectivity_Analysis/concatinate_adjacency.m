function comb_adj =  concatinate_adjacency(adj,nl,index,binary,fractional)
    

    if binary == 1
        neuron_deg = adj
        neuron_deg(neuron_deg>0) = 1
    elseif binary == 0
        neuron_deg = adj;
    else
        error('Incorrect input for binary. 1: yes 0: no')
    end
    size(adj)
    size(neuron_deg)
    comb_adj = zeros(max(index));
    for i = 1:max(index);
        group_ind_i = find(index == i);
        for j = 1:max(index);
            group_ind_j = find(index == j);
            if fractional == 0 
                comb_adj(i,j) = sum(reshape(neuron_deg(group_ind_i,group_ind_j),numel(neuron_deg(group_ind_i,group_ind_j)),1));
            elseif fractional == 1 && binary == 0
                neurons_group_j = nl(group_ind_j);
                input_counts = arrayfun(@(x) vertcat(neurons_group_j(x).Inputs.conid), 1:length(neurons_group_j), 'UniformOutput',false) ;        
                input_total = sum(arrayfun(@(x) length(input_counts{x}), 1:length(input_counts)));
                
                comb_adj(i,j) = sum(reshape(neuron_deg(group_ind_i,group_ind_j),numel(neuron_deg(group_ind_i,group_ind_j)),1))/input_total;
            elseif fractional ==1 && binary == 1
                comb_adj(i,j) = sum(reshape(neuron_deg(group_ind_i,group_ind_j),numel(neuron_deg(group_ind_i,group_ind_j)),1))/numel(neuron_deg(group_ind_i,group_ind_j));

            else
                error('Incorrect input for fractional. 1: yes 0: no')
            end
            clear group_ind_j
        end
        clear group_ind_i
    end