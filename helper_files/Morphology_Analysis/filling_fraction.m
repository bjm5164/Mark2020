function [filling_fraction,n_real,n_possible] = filling_fraction(pre_neuron,post_neuron,density_col)

% Compute a filling fraction for two neurons.  Default settings are to
% compare the filling fraction between the presynaptic axon and
% postsynaptic dendrite. 

if exist('density_col','var') == 0
    density_col = 1
end

min_strahler = 1; 
sn = strahler_number( post_neuron );
%%
% Find the arbors that fit the strahler threshold that are part of the
% dendrite 
input_neurite = find(sn(:)<=min_strahler & post_neuron.input_density(:,density_col) > 5);
post_xyz = post_neuron.Skeleton_Coords(input_neurite,:);

% Next get the index of the presynaptic sites on the axon of the pre-neuron
presyn_node_ids = arrayfun(@(x) pre_neuron.skeleton_data.Tree2Ind(pre_neuron.Outputs.treenodeID(x)), 1:length(pre_neuron.Outputs.treenodeID));
pre_xyz = pre_neuron.Skeleton_Coords(presyn_node_ids,:);


distance_range = 2000;
clear nprox and nfill and ff
for i = 1:length(distance_range)
    D = min( pdist2(post_xyz,pre_xyz), [], 1 );
    D_thresh = D < (distance_range(i));

    n_possible(i) = sum( D_thresh );
    n_real(i) = length(intersect( post_neuron.Inputs.conid, pre_neuron.Outputs.conid( D_thresh ) ) );
    filling_fraction(i) = n_real(i)/n_possible(i);
    if filling_fraction(i) == inf
        filling_fraction(i) = -1;
    elseif isnan(filling_fraction(i)) == 1
        filling_fraction(i) = -1;
    end
end
% figure; hold on
% plot(distance_range(1:end)./1000,ff(1:end),'Marker','o','Color','b','MarkerSize',15,'LineWidth',3)
% set(gcf,'Color','w')
% set(gca,'FontSize',18)
% xlabel('Potential Synapse Radius (µm)')
% ylabel('Filling Fraction')
end
