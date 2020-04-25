figure;hold on; 

scatter3(neurite_coords{1}(:,1),neurite_coords{1}(:,2),neurite_coords{1}(:,3),'k')


scatter3(p1(:,1),p1(:,2),p1(:,3),50,'c')
scatter3(p2(:,1),p2(:,2),p2(:,3),'c')

plot3(neurite_vector(:,1),neurite_vector(:,2),neurite_vector(:,3),'r','LineWidth',3)

%scatter3(neurite_vector(out_cns,1),neurite_vector(out_cns,2),neurite_vector(out_cns,3))
surfaces(CNS_mesh,'k',.02,3)
surfaces(NPM,'k',.02,3)
%%
figure;
histogram(neurite_coords{1}(:,1))
