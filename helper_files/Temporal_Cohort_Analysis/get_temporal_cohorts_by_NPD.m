 function [cohorts] = get_temporal_cohorts_by_NPD(lineage,bandwidth,bilateral,normalization)
        if bilateral == 1
            npd = arrayfun(@(x) lineage(x).Mean_Neuropil_Distance,1:length(lineage)); %Get cortex neurite lengths
        else
            npd = arrayfun(@(x) lineage(x).skeleton_data.Distance_To_Neuropil,1:length(lineage)); %Get cortex neurite lengths
        end
        figure; hold on
        neurite_kernel = fitdist(npd(:),'Kernel','BandWidth',bandwidth); 
        x = 0:.1:max(npd);
        yKernel = pdf(neurite_kernel,x);
        TF = islocalmin(yKernel,'MinSeparation',10000,'MinProminence',1e-5);
        
        plot(x,yKernel,x(TF),yKernel(TF),'r*');
        
        cut_points = x(TF);
        if normalization == 1
            edges = [0,cut_points,1.2]
        else
            edges = [0,cut_points,28000];
        end
        c = discretize(npd,edges);
        for ii = 1:max(c)
            mean_c(ii) = mean(npd(c==ii));
        end
        
        if normalization == 1
            temporal_cohort = discretize(mean_c,[0,.39,.64,.85,1.2]);
        else      
            temporal_cohort = discretize(mean_c,[0,8863,16088,21544,27000]);
        end
        cohorts = zeros(length(c),1);
        for ii = 1:max(c)
            cohorts(c == ii) = temporal_cohort(ii);
        end
 end