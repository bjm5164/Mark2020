function [sim_mat,dists] = synapse_similarity_io_overlap(nl,sigma,omega)
% calculate synapse similarity between inputs and outputs of a population
% of neurons such that sim_mat(i,j) is the similarity between the
% presynaptic sites of neuron i and the postsynaptic sites of neuron j.  
if exist('omega','var') == 0
    omega = sigma
end

for i = 1:length(nl)
    if isempty(nl(i).Inputs.treenodeID) == 0
        Input_Coords{i} = nl(i).Inputs.xyz;
    else
    end
    if isempty(nl(i).Outputs.treenodeID) == 0
        outputs = arrayfun(@(x) repmat(nl(i).Outputs.xyz(x,:),nl(i).Outputs.polyadics(x),1),[1:length(nl(i).Outputs.polyadics)],'UniformOutput',false);
        Output_Coords{i} = cat(1,outputs{:});
    else
    end
end


for i = 1:length(nl)
        D_isis = pdist2(Output_Coords{i},Output_Coords{i},'euclidean');   
    if isempty(D_isis)
        D_isis = inf;
    else
    end
    for j = 1:length(nl)
        clear D_jkjk and D_isjk and k_vals and k_idx 
        D_jkjk = pdist2(Input_Coords{j},Input_Coords{j},'euclidean'); % intraneuron synapse distances for neuron j
   
        if isempty(D_jkjk) == 1
            D_jkjk = inf;     
        else 
        end
        
        D_isjk = pdist2(Output_Coords{i},Input_Coords{j},'euclidean'); % interneuron synapse distances for neurons i,j
  
        if isempty(D_isjk) == 1
            D_isjk = inf;
        else 
        end
            
        [k_vals,k_idx] = min(D_isjk,[],2);  % for every synapse of neuron i, smallest distance to a synapse of neuron j and corresponding index.
        clear sim and n_is and n_jk and Dsk
        for s = 1:size(D_isjk,1);
            Dsk(s) = k_vals(s);
            n_is(s) = numel(find(D_isis(s,:)<omega))/size(D_isis,2);
            n_jk(s) = numel(find(D_jkjk(k_idx(s),:)<omega))/size(D_jkjk,2);
            sim(s) = exp((-1*(Dsk(s)^2))/(2*(sigma^2)))*exp(-1*abs(n_is(s)-n_jk(s))/(n_is(s)+n_jk(s)));
            if i == 1 & j == 1 & s == 1
                dists = Dsk;
            else
                dists = [dists,Dsk];
            end
        end
        similarity_matrix(i,j) = mean(sim);
        
    end
end
similarity_matrix(find(isnan(similarity_matrix))) = 0
sim_mat = (similarity_matrix);
%sim_mat(boolean(eye(length(sim_mat)))) = 1

end

    