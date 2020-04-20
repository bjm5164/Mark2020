function annotations = parse_annotations(Neuron_List,names)

% p = inputParser;
% 
% defaultAnnotations = 'Hemilineage_Side';
% validAnnotation = {'Hemilineage','Lineage','Hemilineage_Side','Lineage_Side'};
% checkAnnotation = @(x) any(validatestring(x,validAnnotation));
% 
% addRequired(p,'Neuron_List');
% addOptional(p,'Annotations',defaultAnnotations,checkAnnotation)
% 
% parse(p,Neuron_List,varargin{:})
% 
% 
% disp(['Annotations: ', p.Results.Annotations])
% 
% if ~isempty(fieldnames(p.Unmatched))
%    disp('Extra inputs:')
%    disp(p.Unmatched)
% end
% if ~isempty(p.UsingDefaults)
%    disp('Using defaults: ')
%    disp(p.UsingDefaults)
% end


for i = 1:length(Neuron_List)
    annots = arrayfun(@(x) ischar(Neuron_List(i).Annotations{x}),1:length(Neuron_List(i).Annotations));
    Neuron_List(i).Annotations(annots == 0) = [];
    % If annotations contain NB, then they are either a MN or an
    % Interneuron
    inter_idx = arrayfun(@(x) find(contains(Neuron_List(i).Annotations(x),'NB')),1:length(Neuron_List(i).Annotations),'UniformOutput',false);
    IN_MN(i) = sum(cat(1,inter_idx{:}))>0;

    % If annotations contain MN, they're an MN
    if ismember(Neuron_List(i).Names{1},'"MNs')
        mns(i) = 1
    else
        mns(i) = contains(Neuron_List(i).Names{1},'MN');
    end


    % If annotations contain Sensory, they're sensory
    sensory_idx = arrayfun(@(x) find(contains(Neuron_List(i).Annotations(x),'Sensory')),1:length(Neuron_List(i).Annotations),'UniformOutput',false);
    sensoryneurons(i) = sum(cat(1,sensory_idx{:}))>0;
end

[interneurons,~] = setdiff(IN_MN,mns,'rows');

for i = 1:length(Neuron_List)
    if interneurons(i) == 1
        dv_index(i) = contains(cat(1,[Neuron_List(i).Annotations{:}]),'Dorsal');
    elseif mns(i) == 1
        dv_index(i) = 2;
    elseif sensoryneurons(i) == 1
        dv_index(i) = 3;
    else
    end
end

interneuron_index = find(interneurons == 1);
for i = 1:length(Neuron_List)
    if interneurons(i) == 1
        f = strfind([Neuron_List(i).Annotations{:}],'NB');
        a =  [Neuron_List(i).Annotations{:}];
        lineages{i} = a(f:f+4);
    elseif mns(i) == 1
        lineages{i} = 'MN';
    elseif sensoryneurons(i) == 1
        lineages{i} = 'SN';
    else
    end
end
unique_lin = unique(lineages);
[~,lin_index] = ismember(lineages,unique_lin);


for i = 1:length(Neuron_List)
    if contains([Neuron_List(i).Annotations{:}],'A1L') == 1 | contains([Neuron_List(i).Annotations{:}],'T3L') == 1
        side_index(i) = 0;
    elseif contains([Neuron_List(i).Annotations{:}],'A1R') == 1 | contains([Neuron_List(i).Annotations{:}],'T3R') == 1
        side_index(i) = 1;
    else
        side_index(i) = nan;
    end
end

if names == 1
    for i = 1:length(Neuron_List)
        if side_index(i) == 0
            side = 'l';
        elseif side_index(i) == 1
            side = 'r';
        else 
        end
   
        if dv_index(i) == 0
            dv = 'Ventral';
        elseif dv_index(i) == 1
            dv = 'Dorsal';
        else 
            dv = '';
        end
         Names{i} = strcat(Neuron_List(i).Annotations{1}(2:end),',',dv,',',side);

    end
    annotations = table(Names(:),dv_index',lineages(:),lin_index',side_index','VariableNames',{'Names','DV_Index','Lineage','Lineage_Index','Side_Index'});
else
    annotations = table(dv_index',lineages(:),lin_index',side_index','VariableNames',{'DV_Index','Lineage','Lineage_Index','Side_Index'});

end
end




    


