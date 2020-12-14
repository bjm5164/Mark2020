load mark2020_neurons_temporal_cohorts.mat
%%
supp_table = an_in(1:2:end,[1,3,6])

supp_table.Hemilineage = ones(height(supp_table),1)
dv_index = table2array(an_in(1:2:end,2))
for i = 1:height(supp_table)
    if dv_index(i) == 0
        hemi{i} = 'Ventral'
    else 
        hemi{i} = 'Dorsal'
    end
end
    
       
supp_table.Hemilineage = hemi'



writetable(supp_table,'/Users/brandon/Desktop/Supp_Table_1.csv')