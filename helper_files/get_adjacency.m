function neuron_adj = get_adjacency(Neuron_List,fractions)
    % Generate an adjacency matrix from a neuron list.  It does this by
    % looking at the partner skids associated with each connector.
    % Inputs:      Neuron_List:  a standard neuron list imported using
    %              CATMAID_IMPORT.
    %              fractions:  1 or 0.  1: divide the edge weight by the
    %              total number of inputs onto that neuron. 0: raw edge
    %              weights.
    % Outputs:  adjacency matrix
    neuron_adj = zeros(length(Neuron_List),length(Neuron_List));
    skids = [Neuron_List(:).SkIDs];

    for i = 1:length(Neuron_List)
        partners = arrayfun(@(x) str2num(Neuron_List(i).Outputs.Partner_skids{x}),1:length(Neuron_List(i).Outputs.conid),'UniformOutput',false);
        ap = transpose(cat(2,partners{:}));
        [C,ia,ib] = unique(ap);
        a_counts = accumarray(ib,1); % Count the number of synapses at each unique set of coordinates
        value_counts = [C, a_counts] ;
        connected = ismember(C,skids);
        con_table = value_counts(connected,:); 
        for j = 1:size(con_table,1)
            if fractions == 1
            neuron_adj(i,find(skids==con_table(j,1))) = ...
            con_table(j,2)/length(Neuron_List(find(skids==con_table(j,1))).Inputs.treenodeID);
            else
            neuron_adj(i,find(skids==con_table(j,1))) = con_table(j,2);
        end
    end
end