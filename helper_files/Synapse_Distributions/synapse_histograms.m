function synapse_histograms(nl,clr,direction,alpha,draw_surface,synapse_index)
%This script makes the synapse distribution plots for figures 2 and 3.  
% Synapse_Distributions(nl,clr,direction,alpha,draw_surface)
% Inputs:   nl:  A neuron list 
%           clr: color
%           direction:  Pre (1) or post (2) or both (3)
%           alpha: synapse alpha
%           draw_surface:  Draw an outline of the neuropil for the a/p.
%           This is an option because the neuropil mesh object can be large
%           when exported as an svg.

% A few details about how this works:  The brain is slightly offset on the
% M/L axis, so to make the left and right hemisegments aligned for the
% cross sectional plots, we rotate all of those points 12 degrees.  The
% synapses are scaled by how many inputs they get.  Basically a pre-synapse
% with one output gets a size of 100, and a pre_synapse with n partners is
% size n*100.  The hisogram plots are 1D kernel density estimates with a 1µm
% bandwith.  


    
load Neuropil_Mesh_Object.mat %Load neuropil mesh for skeletons

for i = 1:length(nl)
    if isempty(nl(i).Inputs.treenodeID) == 0
        if exist('synapse_index') & direction == 2
            Input{i} = nl(i).Inputs.xyz(synapse_index{i});
        else
            Input{i} = nl(i).Inputs.xyz;
        end
    else
        Input{i} = [nan,nan,nan];
    end
    if isempty(nl(i).Outputs.treenodeID) == 0
        if exist('synapse_index') & direction == 1
            output_index = find(synapse_index{i} == 1);
            outputs = arrayfun(@(x) repmat(nl(i).Outputs.xyz(output_index(x),:),nl(i).Outputs.polyadics(output_index(x)),1),[1:length(output_index)],'UniformOutput',false);
        else
            outputs = arrayfun(@(x) repmat(nl(i).Outputs.xyz(x,:),nl(i).Outputs.polyadics(x),1),[1:length(nl(i).Outputs.polyadics)],'UniformOutput',false);
        end    
         Output{i} = cat(1,outputs{:});
    else
        Output{i} = [nan,nan,nan];
    end
end




%% Correct for offset brain
% n_lin = length(nl);
xrot = 0;
yrot = -1;
zrot = -12;


NPM.vertices = rotate_points(NPM.vertices,zrot,3)
NPM.vertices = rotate_points(NPM.vertices,xrot,1)
NPM.vertices = rotate_points(NPM.vertices,yrot,2)

%Find limits each set of axes
axis_lims_MLDV = [min(NPM.vertices(NPM.vertices(:,3)> 105550 & NPM.vertices(:,3)< 144000,1)), max(NPM.vertices(NPM.vertices(:,3)> 105550 & NPM.vertices(:,3)< 144000,1)) , min(NPM.vertices(NPM.vertices(:,3)> 105550 & NPM.vertices(:,3)< 144000,2)), max(NPM.vertices(NPM.vertices(:,3)> 105550 & NPM.vertices(:,3)< 144000,2))];
axis_lims_APDV = [min(NPM.vertices(:,3)), max(NPM.vertices(:,3)),min(NPM.vertices(:,1)), max(NPM.vertices(:,1))];



%% Organize synapse coords
if nargin == 2
direction = input('Presynapses:1 Postsynapses:2 Both:3');
else 
end

    if direction == 1
        % Rotate Coords
        if isempty(nl(i).Outputs.treenodeID) == 0
            Output_Coords = cat(1,Output{:});
            Output_Coords = rotate_points(Output_Coords,zrot,3);
            Output_Coords = rotate_points(Output_Coords,xrot,1);
            Output_Coords = rotate_points(Output_Coords,yrot,2);


            d = 'Presynaptic';
            % Get polyadic information
            [C,ia,ic] = unique(Output_Coords,'rows'); % Search for all unique synapse coordinates
            a_counts = accumarray(ic,1); % Count the number of synapses at each unique set of coordinates
            value_counts = [C, a_counts] ;
            polyadics_ap = value_counts(:,4); % Get the size of each polyadic for plotting.

            apdv_points_counts = Output_Coords(:,[1,3]);
            apdv_points_plot = Output_Coords(ia,[1,3]);     
            % Remove points beyond a2/t3 for cross-sectional analysis
            mldv_points_counts = Output_Coords(apdv_points_counts(:,2) >= 105550 & apdv_points_counts(:,2) <= 144000,1:2);
            mldv_points_plot = Output_Coords(ia,1:2);
            mldv_points_plot = mldv_points_plot(apdv_points_plot(:,2) > 105550 & apdv_points_plot(:,2) < 144000,1:2);
            polyadics_ml = value_counts(apdv_points_plot(:,2) > 105550 & apdv_points_plot(:,2) < 144000,4);
        else
        end
        
    elseif direction == 2
        if isempty(nl(i).Inputs.treenodeID) == 0
            Input_Coords = cat(1,Input{:});
            Input_Coords = rotate_points(Input_Coords,zrot,3);
            Input_Coords = rotate_points(Input_Coords,xrot,1);
            Input_Coords = rotate_points(Input_Coords,yrot,2);

            apdv_points_counts = Input_Coords(:,[1,3]);
            apdv_points_plot = apdv_points_counts;

            mldv_points_counts = Input_Coords(apdv_points_counts(:,2) > 105550 & apdv_points_counts(:,2) < 144000,1:2);
            mldv_points_plot = mldv_points_counts;

            polyadics_ap = ones(length(apdv_points_plot),1)*1;
            polyadics_ml = ones(length(mldv_points_plot),1)*1;
            d = 'Postsynaptic';
        else
        end
    
    elseif direction == 3
        if isempty(nl(i).Outputs.treenodeID) == 0
            Output_Coords = cat(1,Output{:});
            Output_Coords = rotate_points(Output_Coords,zrot,3);
            Output_Coords = rotate_points(Output_Coords,xrot,1);
            Output_Coords = rotate_points(Output_Coords,yrot,2);

            [C,ia,ic] = unique(Output_Coords,'rows'); % Search for all unique synapse coordinates
            a_counts = accumarray(ic,1); % Count the number of synapses at each unique set of coordinates
            value_counts = [C, a_counts] ;
            polyadics_outputs_ap = value_counts(:,4); % Get the size of each polyadic for plotting.

            apdv_outputs_counts = Output_Coords(:,[1,3]);
            apdv_outputs_plot = Output_Coords(ia,[1,3]);     

            mldv_outputs_counts = Output_Coords(apdv_outputs_counts(:,2) > 105550 & apdv_outputs_counts(:,2) < 144000,1:2);
            mldv_outputs_plot = Output_Coords(ia,1:2);
            mldv_outputs_plot = mldv_outputs_plot(apdv_outputs_plot(:,2) > 105550 & apdv_outputs_plot(:,2) < 144000,1:2);
            polyadics_outputs_ml = value_counts(apdv_outputs_plot(:,2) > 105550 & apdv_outputs_plot(:,2) < 144000,4);
        else
        end
        if isempty(nl(i).Inputs.treenodeID) == 0
            Input_Coords = cat(1,Input{:})
            Input_Coords = rotate_points(Input_Coords,zrot,3);
            Input_Coords = rotate_points(Input_Coords,xrot,1);
            Input_Coords = rotate_points(Input_Coords,yrot,2);

            apdv_inputs_counts = Input_Coords(:,[1,3]);
            apdv_inputs_plot = apdv_inputs_counts;

            mldv_inputs_counts = Input_Coords(apdv_inputs_counts(:,2) > 105550 & apdv_inputs_counts(:,2) < 144000,1:2);
            mldv_inputs_plot = mldv_inputs_counts;

            polyadics_inputs_ap = ones(length(apdv_inputs_plot),1)*3;
            polyadics_inputs_ml = ones(length(mldv_inputs_plot),1)*3;
        else
        end
        
        mldv_points_counts = vertcat(mldv_inputs_counts,mldv_outputs_counts);
        mldv_points_plot = vertcat(mldv_inputs_plot,mldv_outputs_plot);
        polyadics_ml = vertcat(polyadics_inputs_ml,polyadics_outputs_ml);
        
        apdv_points_counts = vertcat(apdv_inputs_counts,apdv_outputs_counts);
        apdv_points_plot = vertcat(apdv_inputs_plot,apdv_outputs_plot);
        
        polyadics_ap = vertcat(polyadics_inputs_ap,polyadics_outputs_ap);
         
        d = 'Combined Synaptic ';
    
    else error('Incorrect Direction')
    end
  
    
    
    % Plot the DV histogram
    subplot(3,1,1); hold on
    pd_MLDV_2 = fitdist(mldv_points_counts(:,2)*-1,'kernel','Kernel','normal','BandWidth',1000);
    x2_MLDV = axis_lims_MLDV(4)*-1:200:axis_lims_MLDV(3)*-1;
    y2_MLDV = pdf(pd_MLDV_2,x2_MLDV);
    area(x2_MLDV,y2_MLDV,'EdgeColor',clr,'LineWidth',2,'FaceColor',clr,'FaceAlpha',.1);
    xlim([axis_lims_MLDV(4)*-1,axis_lims_MLDV(3)*-1]);
    ax_dv = gca;
    xticks([])
    yl = ylim;
    xlabel('D/V Density')
    yticks([yl(2)/2,yl(2)]);
    ytickformat('%.0f')
   
    % Plot the M/L histogram
    subplot(3,1,2); hold on
    pd_MLDV_1 = fitdist(mldv_points_counts(:,1),'kernel','kernel','normal','BandWidth',1000);
    x1_MLDV = axis_lims_MLDV(1):200:axis_lims_MLDV(2);
    y1_MLDV = pdf(pd_MLDV_1,x1_MLDV);
    area(x1_MLDV,y1_MLDV,'EdgeColor',clr,'LineWidth',2,'FaceColor',clr,'FaceAlpha',.1)
    xlabel('M/L Density')

  
    ax_ml = gca; ; 
    xticks([]); ax_ml.YAxisLocation = 'Left';
    xlim([axis_lims_MLDV(1),axis_lims_MLDV(2)]);
    yl = ylim;
    yticks([yl(2)/2,yl(2)]);
    ytickformat('%.0f')
    

    % Plot the AP histogram
    ax_ap2 = subplot(3,1,3); hold on
    pd_apdv = fitdist(apdv_points_counts(:,2),'kernel','kernel','normal','BandWidth',2000);
    x1_apdv = axis_lims_APDV(1):800:axis_lims_APDV(2);
    y1_apdv = pdf(pd_apdv,x1_apdv);
    area(x1_apdv,y1_apdv,'EdgeColor',clr,'LineWidth',2,'FaceColor',clr,'FaceAlpha',.1)
    ylabel(strcat(d,' Density'));
    xlim([axis_lims_APDV(1),axis_lims_APDV(2)]); 
    xticks([]);
    xlabel('A/P Density')
    yl = ylim;
    yticks([(yl(2)/2),(yl(2))]);
    ytickformat('%.0f')
    hold off
    
    
    
    % legend(Lins,'Location','SouthWest') ; %legendmarkeradjust(25)
    set(gcf,'Color','w');
    set(findall(gcf,'-property','FontSize'),'FontSize',24);
    
    %linkaxes([ax_ap1,ax_ap2],'x')


end
