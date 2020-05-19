function Plot_Strahler(nl)
for u = 1:length(nl)
    figure
    pg = plot_neurons(nl(u),'k',1,3,1,0,0)
    map = hsv(max(nl(u).Strahler_Number));
    for i = 1:length(map)
        s = find(nl(u).Strahler_Number == i)
        arrayfun(@(x) highlight(pg,x,successors(nl(u).skeleton_data.Graph,x),'EdgeColor',map(i,:)),s)
    end

end
