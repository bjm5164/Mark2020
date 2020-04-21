function [ratiometric_distance] = normalized_neurite_distance(nl,use_mean)
% Calculate a cortex neurite length for each neuron that is normalized to
% the width of the neuropil at the location of that neuron. This can be
% done using the neurite of individual neurons, or all of the neurites of a
% particular neuron.  If use_mean == 1, it will use all neurites of the
% neuron list.  I have only tested this on A1 neurons.  It generates a
% tract vector by using two points: the mean of the tract coordinates and
% the mean +/- a standard deviation. The left and the right side standard
% deviations need to be calculated slightly differently given their axis,
% so this may not be super generalizable. 

% If using the mean, set the neurite coordinates as all of the coordinates
% for the neuron list. 
if use_mean == 1
    for i = 1:length(nl)
        sc{i} = nl(i).skeleton_data.Neurite_xyz
    end
    n_coords = cat(1,sc{:})
    neurite_coords = arrayfun(@(x) n_coords,1:length(nl),'UniformOutput',false)
else
    for i = 1:length(nl)
        neurite_coords{i} = nl(i).skeleton_data.Neurite_xyz
    end
end

% Generate neuropil and CNS triangulations
load Neuropil_Mesh_Object.mat
load CNS_mesh.mat
NP = triangulation(NPM.faces,NPM.vertices);
CP = triangulation(CNS_mesh.faces,CNS_mesh.vertices);

figure; hold on
for i = 1:length(nl)       
    % Generate a vector that is the average of the neurite tract.
    if contains(nl(i).Names{1}(end),'l')
        p1 = mean(neurite_coords{i},1);
        p2 = p1 + std(neurite_coords{i},1);
    elseif contains(nl(i).Names{1}(end),'r')
        p1 = mean(neurite_coords{i},1);
        p2 = p1 + std(neurite_coords{i},1);
        p2(:,2) = p1(:,2)-std(neurite_coords{i}(:,2),1);
        p2(:,3) = p1(:,3)-std(neurite_coords{i}(:,3),1);
    else error('Name not formatted correctly for this analysis')
    end    
    u = (p2-p1)/norm(p2-p1);   % unit vector, p1 to p2
    d = (-50000:50:50000)';       % displacement from p1, along u
    neurite_vector = p1 + d*u;


    
    % Generate distance measures for each point along the neurite vector to
    % either the edge of the neuropil or the edge of the CNS
    clear CP_id and CP_d and NP_id and NP_d
    
    % Find points on neurite vector in the neuropil and in the cns
    in_neuropil = find(inpolyhedron(NPM,neurite_vector)==1);
    out_cns = find(inpolyhedron(CNS_mesh,neurite_vector)==0);
    
     d_neuropil = sqrt((neurite_vector(in_neuropil,1)-p1(1)).^2 + ...
            (neurite_vector(in_neuropil,2)-p1(2)).^2 + ...
             (neurite_vector(in_neuropil,3)-p1(3)).^2);
     neuropil_point = neurite_vector(in_neuropil(d_neuropil == min(d_neuropil)),:);
    
    d_CNS = sqrt((neurite_vector(out_cns,1)-p1(1)).^2 + ...
            (neurite_vector(out_cns,2)-p1(2)).^2 + ...
             (neurite_vector(out_cns,3)-p1(3)).^2);
     cortex_point = neurite_vector(out_cns(d_CNS == min(d_CNS)),:);
    
    cortex_width(i) = pdist2(cortex_point,neuropil_point); 
    npd(i) = nl(i).skeleton_data.Distance_To_Neuropil;
    
    if use_mean == 0
        plot_neurons(nl(i),'k',.1,3,1,0);
        scatter3(neuropil_point(1),neuropil_point(2),neuropil_point(3),100,'c')
        scatter3(cortex_point(1),cortex_point(2),cortex_point(3),100,'g')
        plot3([neuropil_point(1),cortex_point(1)],[neuropil_point(2),cortex_point(2)],[neuropil_point(3),cortex_point(3)]);
    else
    end
        

    
end
if use_mean == 1
    surfaces(NPM,'k',.05,3);
    surfaces(CNS_mesh,'k',.02,3);
    plot_neurons(nl,'k',.1,3,1,0);
    scatter3(neuropil_point(1),neuropil_point(2),neuropil_point(3),100,'c')
    scatter3(cortex_point(1),cortex_point(2),cortex_point(3),100,'g')
    plot3([neuropil_point(1),cortex_point(1)],[neuropil_point(2),cortex_point(2)],[neuropil_point(3),cortex_point(3)]);
    
else
end
ratiometric_distance = npd./cortex_width;
end


