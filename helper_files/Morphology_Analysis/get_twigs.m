function [twig_length,twigs] = get_twigs(neuron,min_strahler)

    if exist('min_strahler','var')
    else
        min_strahler = 1
    end


    sn = strahler_number(neuron)
    
    twig_g = subgraph(graph(neuron.skeleton_data.Adj_dir+neuron.skeleton_data.Adj_dir'),find(sn<=min_strahler));
    twigs = conncomp(twig_g);

    for i = 1:max(twigs)
        twig_temp = subgraph(twig_g,find(twigs==i));
        if size(twig_temp.Edges,1) > 1
            twig_dmat = distances(twig_temp);
            twig_length(i) = max(twig_dmat(:));
        else
            twig_length(i) = nan;
        end
    end
     
%      neuron.Strahler_Number = sn
%      Plot_Strahler(neuron)
%     

end
