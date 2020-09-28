function [h2 T perm1 Z] = Synapse_Distance_Clustering_v2(sim_mat,Names,clusters,dis_sim)
if nargin >= 3
    clusters = clusters;
else 
    clusters = 0;
end
cmap =cbrewer('seq', 'Blues',64)


if nargin == 3
   dis_sim = input('Distance (1) / Similarity (0)')
else
end

if dis_sim == 1
        d_mat = 1-sim_mat;
        dist = squareform(d_mat);  
        mat = d_mat;
else
        dist = pdist(sim_mat,'cosine')
        mat = sim_mat;
end
% eva_gap = evalclusters(sim_mat,'linkage','gap','KList',[1:20])
% eva_sil = evalclusters(sim_mat,'linkage','silhouette','KList',[1:20])
display('doing cluster')
tic, Z = linkage(dist,'average'); toc
leafOrder = optimalleaforder(Z,dist,'Transformation','inverse')

% figure;
[h2 T perm2] = dendrogram(Z,0,'Orientation','Left','ColorThreshold' ,.8,'Labels',Names,'Reorder',leafOrder)

%figure 
%[h1 t1 perm_c] = dendrogram(Z,clusters,'Orientation','Left','ColorThreshold' ,1,'Labels',Names,'Reorder',leafOrder);
figure('rend','painters','pos',[1 1 1940 952]); 
subplot(3,4,[1 5 9]);
display('doing dendrogram');
[h1 t1 perm1] = dendrogram(Z,0,'Orientation','Left','ColorThreshold' ,.8,'Labels',Names,'Reorder',leafOrder);
set(h1,'LineWidth',3)
subplot(3,4,[2 3 4 6 7 8 10 11 12]);
colormap(cmap)
imagesc(mat(perm1,perm1)); axis xy; set(gcf,'Color','w'); axis off; c = colorbar; c.Label.String = 'Synapse Similarity Score'; 
%[h1 t2 perm1] = dendrogram(Z,clusters,'Orientation','Left','ColorThreshold' ,2);
set(findall(gcf,'-property','FontSize'),'FontSize',18)

end