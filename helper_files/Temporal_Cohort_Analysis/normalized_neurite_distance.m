function [ratiometric_distance] = normalized_neurite_distance(nl,use_mean)
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
load Neuropil_Mesh_Object.mat
load CNS_mesh.mat

NP = triangulation(NPM.faces,NPM.vertices);
CP = triangulation(CNS_mesh.faces,CNS_mesh.vertices);

for i = 1:length(nl)       
    % Generate a vector that is the average of the neurite tract.
    if contains(nl(i).Names{1}(end),'l')
        p1 = mean(neurite_coords{i},1);
        p2 = p1 + std(neurite_coords{i},1);
    elseif contains(nl(i).Names{1}(end),'r')
        p1 = mean(neurite_coords{i},1);
        p2 = p1 + std(neurite_coords{i},1);
        p2(:,2) = p1(:,2)-std(neurite_coords{i}(:,2),1);
    else error('Name not formatted correctly for this analysis')
    end
        
        
    u = (p2-p1)/norm(p2-p1);   % unit vector, p1 to p2
    d = (-25000:500:15000)';       % displacement from p1, along u
    neurite_vector = p2 + d*u;

    % Generate neuropil and CNS triangulations
    
    % Generate distance measures for each point along the neurite vector to
    % either the edge of the neuropil or the edge of the CNS
    clear CP_id and CP_d and NP_id and NP_d
    
    % Find points on neurite vector in the neuropil and in the cns
    in_neuropil = inpolyhedron(NPM,neurite_vector)
    in_cns = inpolyhedron(CNS_mesh,neurite_vector)
    
    % Now find the distance between the center of the neurite and the
    % points inside the the neuropil.  The point closest to the center is
    % the entry point of the neurite vector.
    d_neuropil = sqrt((neurite_vector(in_neuropil==1,1)-p1(1)).^2 + ...
                 (neurite_vector(in_neuropil==1,2)-p1(2)).^2 + ...
                 (neurite_vector(in_neuropil==1,3)-p1(3)).^2)
    
    % Now find the distance between the center of the neurite and the
    % points outside the CNS.  The point closest to center is the CNS
    % boundary.
    d_CNS = sqrt((neurite_vector(in_cns==0,1)-p1(1)).^2 + ...
            (neurite_vector(in_cns==0,2)-p1(2)).^2 + ...
            (neurite_vector(in_cns == 0),3-p1(3)).^2)
    
    for ii = 1:length(neurite_vector)
        [CNS_id,CNS_d] = nearestNeighbor(CP,neurite_vector(ii,:));
        [n_ID,neuropil_d] = nearestNeighbor(NP,neurite_vector(ii,:));
        
        CP_id(ii) = CNS_id;
        CP_d(ii) = CNS_d;
        
        NP_id(ii) = n_ID;
        NP_d(ii) = neuropil_d;
    end
    
  cortex_point = neurite_vector(find(CP_d == min(CP_d)),:);
  neuropil_point = neurite_vector(find(NP_d == min(NP_d)),:);
  
  cortex_width(i) = pdist2(cortex_point,neuropil_point);
    
%     if isfield(nl,'Mean_Neuropil_Distance')
%         npd(i) = nl(i).Mean_Neuropil_Distance;
%     else
        npd(i) = nl(i).skeleton_data.Distance_To_Neuropil;
%    end
    
end
figure; hold on
surfaces(NPM,'k',.05,3);
surfaces(CNS_mesh,'k',.02,3);
plot_neurons(nl,'k',.1,3,1,0);
scatter3(neuropil_point(1),neuropil_point(2),neuropil_point(3),100,'c')
scatter3(cortex_point(1),cortex_point(2),cortex_point(3),100,'g')
plot3([neuropil_point(1),cortex_point(1)],[neuropil_point(2),cortex_point(2)],[neuropil_point(3),cortex_point(3)]);
ratiometric_distance = npd./cortex_width;
end


