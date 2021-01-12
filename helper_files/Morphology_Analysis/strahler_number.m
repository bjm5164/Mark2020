
function sn = strahler_number( Neuron )
% Code adapted from: 10.7554/eLife.12059
sn = zeros( size( Neuron.skeleton_data.Adj_dir,1 ),1);

Adj_dir = spones(Neuron.skeleton_data.Adj_dir);
is_branch = sum(Adj_dir,1)>1;
is_branch = full(is_branch);
work_inds = find( sum(Adj_dir,1) == 0);

while ~isempty(work_inds)
    
    rel_ind = work_inds(1);
    work_inds(1) = [];
    
    ch_inds = find(Adj_dir( :, rel_ind ));
    if isempty(ch_inds)
        sn(rel_ind) = 1;
    elseif length(ch_inds)==1
        sn(rel_ind) = sn(ch_inds);
    elseif any( sn(ch_inds)==0 )
        work_inds(end+1) = rel_ind;
        continue
    else
        sn(rel_ind) = max( sn(ch_inds) ) + ( sum( sn(ch_inds)==max(sn(ch_inds)) )>1 );
    end
    
    while true
        rel_ind = find( Adj_dir(rel_ind,:) ); % Look for the parent
        
        if isempty(rel_ind)
            break
        else
            if is_branch( rel_ind )
                work_inds(end+1) = rel_ind;
                break
            else
                sn(rel_ind) = sn(Adir_uw(:,rel_ind)>0);
            end
        end
     
    end
    
    
end
end
