%在所有参数中选出AUPR最大的参数
%同时给出该参数下的AUPR，AUC
%同时给出最大的AUC
neighbor_alpha=0.1:0.1:0.9;
LP_alpha=0.1:0.1:0.9;
lncRNA_alpha=0:0.05:1.0;
max_aupr = 0;
max_auc = 0;
best_neighbor_alpha = 0;
best_LP_alpha = 0;
best_lncRNA_alpha = 0;
count = 0;
all_aupr = zeros(length(neighbor_alpha)*length(LP_alpha)*length(lncRNA_alpha), 4);

for i =1:length(neighbor_alpha)
    for j =1:length(LP_alpha)
        for k =1:length(lncRNA_alpha)
            count = count + 1;
            load(strcat('LNSLP_', num2str(neighbor_alpha(i)), '_', num2str(LP_alpha(j)), '_', num2str(lncRNA_alpha(k)), '.mat'));
            all_aupr(count,:) = [final_result(1), neighbor_alpha(i), LP_alpha(j), lncRNA_alpha(k)];
            if (final_result(1)>max_aupr)
                best_neighbor_alpha = neighbor_alpha(i);
                best_LP_alpha = LP_alpha(j);
                best_lncRNA_alpha = lncRNA_alpha(k);
                max_aupr=final_result(1);
            end
            if (final_result(2)>max_auc)
                max_auc=final_result(2);
            end
        end
    end
end
best_neighbor_alpha
best_LP_alpha
best_lncRNA_alpha
max_aupr
load(strcat('LNSLP_', num2str(best_neighbor_alpha), '_', num2str(best_LP_alpha), '_', num2str(best_lncRNA_alpha), '.mat'));
final_result(2)
max_auc