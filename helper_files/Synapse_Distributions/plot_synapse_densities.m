function plot_synapse_densities(nl_cell,clr,direction,threshold_cutoff)
    
if ischar(clr) == 1
        color_table = containers.Map({'b','k','r','g','y','c','m','w'},{[0,0,1],[0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,1,1],[1,0,1],[1,1,1]});
        clr = color_table(clr);
    elseif iscell(clr)
        color_table = containers.Map({'b','k','r','g','y','c','m','w'},{[0,0,1],[0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,1,1],[1,0,1],[1,1,1]});
        for i = 1:length(clr)
            if ischar(clr{i}) == 1
                colors(i,:) = color_table(clr{i});
            else
                colors(i,:) = clr{i};
            end
        end
        clr = colors;
    else
end





if iscell(nl_cell)
    n = length(nl_cell)

    figure('pos',[1975 153 1375 715],'rend','painters'); hold on
    if n > 1
        if iscell(clr)
            synapse_density_distributions(nl_cell{1},clr{1},direction,threshold_cutoff,1)
            for i = 2:n
                synapse_density_distributions(nl_cell{i},clr{i},direction,threshold_cutoff,0)
            end
        else
            synapse_density_distributions(nl_cell{1},clr(1,:),direction,threshold_cutoff,1)
            for i = 2:n
                synapse_density_distributions(nl_cell{i},clr(i,:),direction,threshold_cutoff,0)
            end
        end
    else
        synapse_density_distributions(nl_cell,clr,direction,threshold_cutoff,1)
    end
else
    synapse_density_distributions(nl_cell,clr,direction,threshold_cutoff,1)
end
    
    ax_1 = subplot(3,7,[16 17 18]);
    pos = ax_1.Position;
    posnew = pos; posnew(4) = posnew(4) +.08; set(ax_1, 'Position', posnew); clear pos and posnew;
    ax_1.YLabel.Position(2) = ax_1.YLabel.Position(2) + .00005
    
    
    ax_2 = subplot(3,7,[12,13,14,19,20,21]); hold on;
    pos = ax_2.Position;
    posnew = pos; posnew(1) = posnew(1) + 0.08; posnew(4) = posnew(4) +.08; set(ax_2, 'Position', posnew); clear pos and posnew;
    
    ax_3 = subplot(3,7,[5,6,7]); hold on
    pos = ax_3.Position
    posnew = pos; posnew(1) = posnew(1) + 0.08; set(ax_3, 'Position', posnew); clear pos and posnew
end
        
        