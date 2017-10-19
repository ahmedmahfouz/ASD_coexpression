
clear all;

%%%------------------------------------------------------------------------
%%% Calculate the total sum of reads per sample (brain structure) - from
%%% the original excel file including 580 samples for the HK genes
%%%------------------------------------------------------------------------
% load('files\origExpMat');
% load('files\hkGeneIDs.mat');
% 
% [hkSL_IDs hkSL_gNames] = xlsread('files\DataFiles\HKshortlist.xls');
% 
% readsPerSample_ALL = sum(origExpMat(:, 2:end));
% readsPerSample_HK = sum(origExpMat(hkIDs, 2:end));
% readsPerSample_HKSL = sum(origExpMat(hkSL_IDs, 2:end));
% 
% figure, hold on
% % bar(readsPerSample_ALL, 'blue'), grid on; 
% % bar(readsPerSample_HK, 'r'), grid on;
% bar(readsPerSample_HKSL, 'g'), grid on;
% xlabel('samples'); ylabel('RPKM Sum');
% title('RPKM Sum per Sample'); %legend('All Genes (22327)', 'Housekeeping Genes (1830)', 'Housekeeping Genes (25 Genes)');
% hold off
%%%------------------------------------------------------------------------

%%%------------------------------------------------------------------------
%%% Plot the expression profile of the long list of HK genes (1830 genes)
%%% provided by Chang et. al
%%%------------------------------------------------------------------------
% load('files\origExpMat');
% load('files\hkGeneIDs.mat');
% 
% expMat_HK = origExpMat(hkIDs, 2:end);
% 
% colors = hsv(length(hkIDs));
% figure, hold on
% for i = 1 : length(hkIDs)
%     
%     plot(expMat_HK(i,:), 'Color', colors(i,:), 'LineWidth', 2), grid on
%     xlabel('Samples', 'fontweight', 'bold'); ylabel('RPKM', 'fontweight', 'bold');
%     title('Expression Profile of 1830 HK Genes', 'fontweight', 'bold'); 
%     
% end
% hold off
%%%-----------------------------------------------------------------------

%%%------------------------------------------------------------------------
%%% Plot the expression profile of the short list of HK genes (25 genes)
%%% provided by ark Ziats
%%%------------------------------------------------------------------------
% load('files\origExpMat');
% 
% [hkSL_IDs hkSL_gNames] = xlsread('files\DataFiles\HKshortlist.xls');
% 
% hkSL_IDs(3:7) = [];
% 
% expMat_HKSL = origExpMat(hkSL_IDs, 2:end);
% 
% colors = hsv(length(hkSL_IDs));
% figure, hold on
% for i = 1 : length(hkSL_IDs)
%     
%     plot(expMat_HKSL(i,:), 'Color', colors(i,:), 'LineWidth', 2), grid on
%     xlabel('Samples', 'fontweight', 'bold'); ylabel('RPKM', 'fontweight', 'bold');
%     title('Expression Profile of 25 HK Genes', 'fontweight', 'bold'); 
%     
% end
% % legend('PSMB2', 'VPS72', 'RAC1', 'EIF1', 'ATP5A1', 'GPS1', 'H3F3A', ...
% %     'AAMP', 'HSP90AB1', 'CENPB', 'PSMB3', 'PTDSS1', 'ESD', 'LAMB1', ...
% %     'POLR2A', 'CANX', 'HNRNPC', 'TCEA1P2', 'COL6A1', 'HMGB1')
% legend('PSMB2', 'VPS72', ...
%     'AAMP', 'HSP90AB1', 'CENPB', 'PSMB3', 'PTDSS1', 'ESD', 'LAMB1', ...
%     'POLR2A', 'CANX', 'HNRNPC', 'TCEA1P2', 'COL6A1', 'HMGB1')
% hold off
%%%-----------------------------------------------------------------------



%%%------------------------------------------------------------------------
%%% Plot the mean and std of the expression profile of the short list of HK 
%%% genes (25 genes) provided by Mark Ziats
%%%------------------------------------------------------------------------
load('files\origExpMat');

[hkSL_IDs hkSL_gNames] = xlsread('files\DataFiles\HKshortlist.xls');

% hkSL_IDs(3:7) = [];

expMat_HKSL = origExpMat(hkSL_IDs, 2:end);

hkM = mean(expMat_HKSL, 1);
hkSTD = std(expMat_HKSL, 1);

xlswrite('hk_M&STD.xls', hkM', 1, 'B2');
xlswrite('hk_M&STD.xls', hkSTD', 1, 'C2');

figure,
boxplot(expMat_HKSL, 'plotstyle', 'compact'), grid on
% xlabel('Samples', 'fontweight', 'bold'); ylabel('RPKM', 'fontweight', 'bold');
title('Expression Profile of 25 HK Genes', 'fontweight', 'bold');
%%%-----------------------------------------------------------------------



