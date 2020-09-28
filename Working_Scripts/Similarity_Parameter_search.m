load mark2020_neurons_temporal_cohorts.mat
[sim_mat,dists] = synapse_similarity_io_overlap(nl,2000,5000);
%%
figure; 
histogram(dists,10,'Normalization','probability')
xlabel('Distance between pre/post synapse (µm)')
ylabel('Frequency')

%%
for i = 1:length(nl)
    tl = get_twigs(nl(i),2);
    if i == 1
        twig_lengths = tl;
    else
        twig_lengths = [twig_lengths,tl];
    end
    clear tl
    
end
%% 
twig_lengths(isnan(twig_lengths)) = []

figure; hold on
bar(1,mean(dists),'FaceColor',[.2 .2 .2])
errorbar(1,mean(dists),std(dists),'k')
bar(2,mean(twig_lengths),'FaceColor',[.5 .5 .5])
errorbar(2,mean(twig_lengths),std(twig_lengths),'k')

figure; 
histogram(twig_lengths,20)
set(gca,'YScale','log')


%%
D = 0:1:3000;
sigma = 0:500:5000;

D = D./1000;
sigma = sigma./1000
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
ylabel('Similarity')
xlabel('Distance')
set(gca,'FontSize',18)
subplot(2,4,[4,8])
ylabel('Frequency')
set(gca,'FontSize',18)
c = colorbar
c.Label.String = 'Sigma'
c.Ticks=[.1:.1:1]
c.TickLabels=[.1:.1:10]

%%
sigma = 1000:1000:5000
map = plasma(length(sigma))
figure; hold on
for i = 1:length(sigma)
    sim_mat(:,:,i) = synapse_similarity_io_overlap(nl,5000,sigma(i))
    histogram(sim_mat(:,:,i),0:.01:1,'FaceColor',map(i,:),'FaceAlpha',.15)
end