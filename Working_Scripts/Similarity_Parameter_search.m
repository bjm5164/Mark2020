d_vals = reshape(dists,numel(dists),1);
d_m = []
for i = 1:length(d_vals)
    d_1 = cat(2,d_vals{i}{:});
    d_m = [d_m,d_1];
end

%%
D = 0:10:6000;
sigma = 1000:1000:10000;
map = plasma(length(sigma));
figure; 
for i = 1:length(sigma)
  kernel = exp((-1*(D.^2))./(2*(sigma(i)^2)));
  subplot(2,4,[1,2,3,5,6,7]); hold on
  plot(D,kernel,'Color',map(i,:),'LineWidth',5);
  subplot(2,4,[4,8]); hold on
  histogram(kernel,0:.01:1,'FaceColor',map(i,:),'Normalization','probability')
  dev(i) = std(kernel)
  ax = gca;
  ax.XAxisLocation = 'Bottom'; xticks([]);
  view([-270, -90]); 
  
end

colormap(map)

subplot(2,4,[1,2,3,5,6,7])
ylabel('Convolved Distance')
xlabel('Distance')
set(gca,'FontSize',18)
subplot(2,4,[4,8])
ylabel('Frequency')
set(gca,'FontSize',18)
c = colorbar
c.Label.String = 'Sigma'
c.Ticks=[0:.1:1]
c.TickLabels=[1000:1000:10000]

%%
sigma = 1000:1000:5000
map = plasma(length(sigma))
figure; hold on
for i = 1:length(sigma)
    sim_mat(:,:,i) = synapse_similarity_io_overlap(nl,5000,sigma(i))
    histogram(sim_mat(:,:,i),0:.01:1,'FaceColor',map(i,:),'FaceAlpha',.15)
end