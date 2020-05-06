function pg = plot_neurons(Neurons,color,alpha,dim,soma,synapses,strahler_sizing,soma_color)

for i = 1:length(Neurons)
    
        
    hold on
    if dim == 3
       pg = plot(Neurons(i).skeleton_data.Graph,'XData',Neurons(i).Skeleton_Coords(:,1),'YData',Neurons(i).Skeleton_Coords(:,2),'ZData',Neurons(i).Skeleton_Coords(:,3),'EdgeColor',color,'NodeColor',color,'MarkerSize',.001,'LineWidth',1,'EdgeAlpha',alpha,'ShowArrows','Off');
       if nargin > 6 
           if strahler_sizing == 1
               sn = strahler_number(Neurons(i));
                for kk = 1:length(sn)
                    highlight(pg,kk,successors(Neurons(i).skeleton_data.Graph,kk),'LineWidth',sn(kk));
                    highlight(pg,kk,'MarkerSize',sn(kk));
                end
           else
           end
       else
       end
       
        if exist('soma_color','var')
        else
            soma_color = color
        end
        if soma == 1
            [X, Y, Z] = sphere(300);
            X = 2000*X + Neurons(i).Soma_Coords(1);
            Y= 2000*Y + Neurons(i).Soma_Coords(2);
            Z = 2000*Z + Neurons(i).Soma_Coords(3);
            surf(X,Y,Z,'FaceColor',soma_color,'LineStyle','none','FaceAlpha',1);
            
        else 
        end
        
    elseif dim == 2
        pg = plot(Neurons(i).skeleton_data.Graph,'XData',Neurons(i).Skeleton_Coords(:,1),'YData',Neurons(i).Skeleton_Coords(:,2),'EdgeColor',color,'NodeColor',color,'MarkerSize',.001,'LineWidth',2,'EdgeAlpha',alpha,'ShowArrows','Off');
    
        if soma == 1
%             scatter(Neurons(i).Soma_Coords(1),Neurons(i).Soma_Coords(2),5000,'.','MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha);
%             camlight 
               [X, Y, Z] = sphere(300);
               X = 2000*X + Neurons(i).Soma_Coords(1);
               Y= 2000*Y + Neurons(i).Soma_Coords(2);
               Z = 2000*Z + Neurons(i).Soma_Coords(3);
               surf(X,Y,Z,'FaceColor',color,'LineStyle','none','FaceAlpha',alpha);
               
        else
        end
    elseif contains(dim,'d')
        pg = plot(Neurons(i).skeleton_data.Graph,'XData',Neurons(i).Skeleton_Coords(:,3),'YData',Neurons(i).Skeleton_Coords(:,1),'EdgeColor',color,'NodeColor',color,'MarkerSize',.001,'LineWidth',2,'EdgeAlpha',alpha,'ShowArrows','Off');
        if soma == 1
            scatter(Neurons(i).Soma_Coords(3),Neurons(i).Soma_Coords(1),1000,'.','MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha);
             
        else
        end
    else
    end
    end
    if nargin >5
        if synapses == 1
             
%                 for u = 1:length(Neurons(i).Inputs.treenodeID)
%                     highlight(pg,Neurons(i).Tree2Ind(Neurons(i).Inputs.treenodeID(u)),'NodeColor','c','MarkerSize',5)
%                 end
%                 for k = 1:length(Neurons(i).Outputs.treenodeID)
%                     highlight(pg,Neurons(i).Tree2Ind(Neurons(i).Outputs.treenodeID(k)),'NodeColor','r','MarkerSize',5)
%                 end
                if isempty(Neurons(i).Inputs.treenodeID) == 0
                    if dim == 3
                        scatter3(Neurons(i).Inputs.xyz(:,1),Neurons(i).Inputs.xyz(:,2),Neurons(i).Inputs.xyz(:,3),200,'.','MarkerFaceColor','c','MarkerEdgeColor','c')
                    else
                        scatter(Neurons(i).Inputs.xyz(:,1),Neurons(i).Inputs.xyz(:,2),200,'.','MarkerFaceColor','c','MarkerEdgeColor','c')
                    end
                else
                end
                if isempty(Neurons(i).Outputs.treenodeID) == 0
                    if dim == 3
                        scatter3(Neurons(i).Outputs.xyz(:,1),Neurons(i).Outputs.xyz(:,2),Neurons(i).Outputs.xyz(:,3),200,'.','MarkerFaceColor','r','MarkerEdgeColor','r')
                    else
                        scatter(Neurons(i).Outputs.xyz(:,1),Neurons(i).Outputs.xyz(:,2),200,'.','MarkerFaceColor','r','MarkerEdgeColor','r')
                    end
                else
                end

             
        else
        end
    else
    end    

set(gca,'Color','w')
end
