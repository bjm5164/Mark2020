function [cohorts,t] = get_temporal_cohorts_by_cell_body(lineage,cutoff,bilateral,normalization)
    if  nargin>2 & bilateral == 1
        an_l = parse_annotations(lineage,1);
        index_l = find(an_l.Side_Index == 0);
        index_r = find(an_l.Side_Index == 1);
        for i = 1:length(index_l)
            soma_coords_l(i,:) = lineage(index_l(i)).Soma_Coords;
            soma_coords_r(i,:) = lineage(index_r(i)).Soma_Coords;
        end
        
        if normalization == 1
            npd_l = normalized_neurite_distance(lineage(index_l),1)
            npd_r = normalized_neurite_distance(lineage(index_r),1)
        else
             npd_l = arrayfun(@(x) lineage(index_l(x)).skeleton_data.Distance_To_Neuropil,1:length(index_l)); %Get cortex neurite lengths
             npd_r = arrayfun(@(x) lineage(index_r(x)).skeleton_data.Distance_To_Neuropil,1:length(index_r));
        end
           
        npd = .5*(npd_l+npd_r); % Get mean cortex neurite length
        npd_bi = zeros(length(npd)*2,1);
        npd_bi(index_l) = npd;
        npd_bi(index_r) = npd;
        npd = .5*(npd_l+npd_r);
        d_l = pdist(soma_coords_l,'euclidean'); %Get distances between somas    
        d_r = pdist(soma_coords_r,'euclidean');
        d = .5*(d_l + d_r);
    else
        % Get soma coordinates
        for i = 1:length(lineage)
            soma_coords(i,:) = lineage(i).Soma_Coords;
        end
        if normalization == 1
            npd = normalized_neurite_distance(lineage,1)
        else
            npd = arrayfun(@(x) lineage(x).skeleton_data.Distance_To_Neuropil,1:length(lineage)); %Get cortex neurite lengths
        end
        d = pdist(soma_coords,'euclidean'); %Get distances between somas
    end
    
    Z = linkage(d,'ward'); % Cluster based on soma positions
    t = cluster(Z,'cutoff',cutoff,'criterion','distance'); %Set cutoff for clusters
%      figure;
%      if bilateral == 1
%         dendrogram(Z,'Labels',[lineage(1:2:end).Names],'ColorThreshold',cutoff)
%      else
%         dendrogram(Z,'Labels',[lineage.Names],'ColorThreshold',cutoff)
%      end
    for i = 1:max(t)
        cluster_distance(i) = mean(npd(t==i)); % Get mean cluster neurite length
    end
    if normalization == 1
        c_edge = sort(discretize(cluster_distance,[0,.39,.64,.85,1.2]));
    else
        c_edge = sort(discretize(cluster_distance,[0,8863,16088,21544,27000]));
    end
    [~,index] = sort(cluster_distance); %Sort clusters by neurite length
    if bilateral == 1
        t = kron(t,ones(2,1));
    else
    end

    for i = 1:length(lineage)
        cohorts(i) = c_edge(ismember(index,t(i))); %get temporal cohorts
    end
    map = cool(4);
    figure; hold on
    for i = 1:max(c_edge)
        plot_neurons(lineage(cohorts==i),map(i,:),1,3,1,0)
    end
    view([190 90]);
    camlight ;
    axis off;
    set(gcf,'Color','w');
    xlim([.8e4,11e4])
    ylim([5.75e4,12e4])  
end