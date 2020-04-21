function [contours] = synapse_density_contours(nl,clr,direction,contour_number,axes_to_use,threshold_cutoff)
load Neuropil_Mesh_Object.mat %Load neuropil mesh for skeletons
%load All_Synapse_Coords; % Load synase coords for A1 for comparison to total.

for i = 1:length(nl)
    if isempty(nl(i).Inputs.treenodeID) == 0
        Input{i} = nl(i).Inputs.xyz;
    else
        Input{i} = [nan,nan,nan]
    end
    if isempty(nl(i).Outputs.treenodeID) == 0
        outputs = arrayfun(@(x) repmat(nl(i).Outputs.xyz(x,:),nl(i).Outputs.polyadics(x),1),[1:length(nl(i).Outputs.polyadics)],'UniformOutput',false);
        Output{i} = cat(1,outputs{:});
    else
        Output{i} = [nan,nan,nan]
    end
end


%% Correct for offset brain
n_lin = length(nl)
xrot = 0;
yrot = -1;
zrot = -12;

NPM.vertices = rotate_points(NPM.vertices,zrot,3);
NPM.vertices = rotate_points(NPM.vertices,xrot,1);
NPM.vertices = rotate_points(NPM.vertices,yrot,2);

%% Make Plots

 if direction == 1
        Output_Coords = cat(1,Output{:})
        Output_Coords = rotate_points(Output_Coords,zrot,3);
        Output_Coords = rotate_points(Output_Coords,xrot,1);
        Output_Coords = rotate_points(Output_Coords,yrot,2);
        
                    
        d = 'Presynaptic';
        [C,ia,ic] = unique(Output_Coords,'rows'); % Search for all unique synapse coordinates
        a_counts = accumarray(ic,1); % Count the number of synapses at each unique set of coordinates
        value_counts = [C, a_counts] ;
        polyadics_ap = value_counts(:,4) % Get the size of each polyadic for plotting.
        
        
        apdv_points_counts = Output_Coords(:,[1,3]);
        apdv_points_plot = Output_Coords(ia,[1,3]);     
        
        mldv_points_counts = Output_Coords(apdv_points_counts(:,2) >= 105550 & apdv_points_counts(:,2) <= 144000,1:2)
        mldv_points_plot = Output_Coords(ia,1:2);
        mldv_points_plot = mldv_points_plot(apdv_points_plot(:,2) > 105550 & apdv_points_plot(:,2) < 144000,1:2)
        polyadics_ml = value_counts(apdv_points_plot(:,2) > 105550 & apdv_points_plot(:,2) < 144000,4)
        
        
    
    elseif direction == 2
        Input_Coords = cat(1,Input{:})
        Input_Coords = rotate_points(Input_Coords,zrot,3);
        Input_Coords = rotate_points(Input_Coords,xrot,1);
        Input_Coords = rotate_points(Input_Coords,yrot,2);
        
        apdv_points_counts = Input_Coords(:,[1,3]);
        apdv_points_plot = apdv_points_counts;
        
        mldv_points_counts = Input_Coords(apdv_points_counts(:,2) > 105550 & apdv_points_counts(:,2) < 144000,1:2);
        mldv_points_plot = mldv_points_counts;
        
        polyadics_ap = ones(length(apdv_points_plot),1);
        polyadics_ml = ones(length(mldv_points_plot),1);
        d = 'Postsynaptic';
    
    elseif direction == 3
        Output_Coords = cat(1,Output{:})
        Output_Coords = rotate_points(Output_Coords,zrot,3);
        Output_Coords = rotate_points(Output_Coords,xrot,1);
        Output_Coords = rotate_points(Output_Coords,yrot,2);
        
        [C,ia,ic] = unique(Output_Coords,'rows'); % Search for all unique synapse coordinates
        a_counts = accumarray(ic,1); % Count the number of synapses at each unique set of coordinates
        value_counts = [C, a_counts] ;
        polyadics_outputs_ap = value_counts(:,4) % Get the size of each polyadic for plotting.
        
        
        apdv_outputs_counts = Output_Coords(:,[1,3]);
        apdv_outputs_plot = Output_Coords(ia,[1,3]);     
        
        mldv_outputs_counts = Output_Coords(apdv_outputs_counts(:,2) > 105550 & apdv_outputs_counts(:,2) < 144000,1:2)
        mldv_outputs_plot = Output_Coords(ia,1:2);
        mldv_outputs_plot = mldv_outputs_plot(apdv_outputs_plot(:,2) > 105550 & apdv_outputs_plot(:,2) < 144000,1:2)
        polyadics_outputs_ml = value_counts(apdv_outputs_plot(:,2) > 105550 & apdv_outputs_plot(:,2) < 144000,4)
        
        
        Input_Coords = cat(1,Input{:})
        Input_Coords = rotate_points(Input_Coords,zrot,3);
        Input_Coords = rotate_points(Input_Coords,xrot,1);
        Input_Coords = rotate_points(Input_Coords,yrot,2);
        
        apdv_inputs_counts = Input_Coords(:,[1,3]);
        apdv_inputs_plot = apdv_inputs_counts;
        
        mldv_inputs_counts = Input_Coords(apdv_inputs_counts(:,2) > 105550 & apdv_inputs_counts(:,2) < 144000,1:2);
        mldv_inputs_plot = mldv_inputs_counts;
        
        polyadics_inputs_ap = ones(length(apdv_inputs_plot),1);
        polyadics_inputs_ml = ones(length(mldv_inputs_plot),1);
        

        mldv_points_counts = vertcat(mldv_inputs_counts,mldv_outputs_counts);
        mldv_points_plot = vertcat(mldv_inputs_plot,mldv_outputs_plot);
        polyadics_ml = vertcat(polyadics_inputs_ml,polyadics_outputs_ml);
        
        
        apdv_points_counts = vertcat(apdv_inputs_counts,apdv_outputs_counts);
        apdv_points_plot = vertcat(apdv_inputs_plot,apdv_outputs_plot);
        
        polyadics_ap = vertcat(polyadics_inputs_ap,polyadics_outputs_ap);
        
        
        d = 'Combined Synaptic ';
    
    else error('Incorrect Direction')
    end




if axes_to_use == 1
    axis_lims= [min(NPM.vertices(NPM.vertices(:,3)>.9e5,1)), max(NPM.vertices(NPM.vertices(:,3)>.9e5,1)) , min(NPM.vertices(NPM.vertices(:,3)>.9e5,2)), max(NPM.vertices(NPM.vertices(:,3)>.9e5,2))];    axis_index = [1,2];
    Synapses = mldv_points_counts;
elseif axes_to_use == 2
    axis_lims = [min(NPM.vertices(:,3)), max(NPM.vertices(:,3)),min(NPM.vertices(:,1)), max(NPM.vertices(:,1))];
    axis_index = [1,3]
    P = NPM.vertices(:,[1,3]);
    k = boundary(P,.4);
    Synapses = apdv_points_counts;
else error('Wrong Axes')
end

%%

hold on
if axes_to_use == 1
    %surfaces([NPM.v(NPM.v(:,3)>.9e5,1),NPM.v(NPM.v(:,3)>.9e5,2)],'k',.05,'-');
    mesh_x = axis_lims(1):100:axis_lims(2);
    mesh_y = axis_lims(3):100:axis_lims(4);
else
    surfaces(NPM.v(:,[1,3]),'k',.05,'-'); hold on
    mesh_x = min(NPM.v(:,1)):100:max(NPM.v(:,1));
    mesh_y = min(NPM.v(:,3)):100:max(NPM.v(:,3));
end


[x1,x2] = meshgrid(mesh_x,mesh_y);

x_1 = x1(:);
x_2 = x2(:);
xi = [x_1 x_2];

[f,xi,bw] = ksdensity(Synapses,xi,'Bandwidth',1000);

f_cut = prctile(f(f>0),threshold_cutoff);
f(f<f_cut) = 0;


if contour_number == 1
    [C,h] = contour(x1,x2,reshape(f,size(x1)),[f_cut f_cut]);
    h.LineColor = clr;
else
    [C,h] = contour(x1,x2,reshape(f,size(x1)),min(f(f>0)):(max(f))/contour_number:max(f));
    h.Fill = 'On';
end
    


% figure; subplot(2,1,1)
% [C,h] = contour(x1,x2,reshape(f,size(x1)),10); h.Fill = 'On'
% subplot(2,1,2)
% [C,h] = contour(x1,x2,reshape(f,size(x1)),f_cutt:(max(f)-f_cutt)/contour_number:max(f)); h.Fill = 'On';





view([180,90])
set(gcf,'Color','w')
set(gca,'FontSize',14)


end
