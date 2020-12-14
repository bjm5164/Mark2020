cohort_sim(similarity(:,:,1),tc_index)
%%

for i = 1:1000
    
    rand_index = tc_index(randi(length(tc_index),length(tc_index),1));
    cs_rand(i) = cohort_sim(similarity(:,:,1),rand_index);
end
%%
function sd = cohort_sim(sim,index)

   
    for i = 1:max(index)
        clear temporal_cohort_pairs temporal_cohort_pairs
        if length(find(index==i)) > 1
            temporal_cohort_pairs = nchoosek(find(index==i),2);
            for j = 1:size(temporal_cohort_pairs,1);

                n_1_sim = sim(temporal_cohort_pairs(j,1),:);
                n_2_sim = sim(temporal_cohort_pairs(j,1),:);
                
                sim_dot{i}(j) = dot(n_1_sim,n_2_sim);
            end
        else
                   
        end
        
      
    

    end

 sd = mean(cat(2,sim_dot{:}));   

end