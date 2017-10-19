%%% Analyze the molecular Interaction between genes (gene correlation
%%% networks) according to the new design of the experiments (11 July 2012)

% clear all;
filesFolder = 'E:\Ahmed\HP\work\Data\ASD_Coexpression_Network_Analysis\';
resultsFolder = 'E:\Ahmed\HP\work\Results\ASD_Coexpression_Network_Analysis\';

%%% Gene Names
load([filesFolder 'files/gNames.mat']);
load([filesFolder 'files/genesStatus_5RPKM.mat']);
gNames(find(genesStatus_5RPKM == 0)) = [];

%%% Structures Names
load([filesFolder 'files/strucs.mat']);
strucIndInc = [2, 3, 5, 6, 7, 8, 9, 15, 18, 19, 20, 21, 22, 23, 24, 26];
for i = 1 : length(strucIndInc)
    S{i} = strucs{strucIndInc(i)};
end

%%% Array of reduced structures (11 NCX structures combined)
nonCs = [5,7,10,12,16];
Cs = [1,2,3,4,6,8,9,11,13,14,15];
newS = S(nonCs);
newS(length(newS)+1) = {'NCX'};

%%% Order structureNames
orderedStrucInd = [nonCs, Cs];
ordS = S(orderedStrucInd);

group{1} = [1, 2, 3 ,4];
group{2} = [5, 6, 7, 8];
group{3} = [9, 10, 11, 12];
group{4} = [13, 14, 15, 16];
group{5} = [17, 18, 19, 20];
group{6} = [21, 22, 23, 24];
group{7} = [25, 26, 27, 28, 29, 30];


%%% -----------------------------------------------------------------------
%%% Experiment#1: Cluster genes based on the average correlation of all
%%% donors together  
%%% -----------------------------------------------------------------------
%%% (1) load the averaged correlation matrix (Spearman Correlation)
% load(['results/geneCorr/Avg(log2+epsilon)/All_corrMat.mat']);
% 
% %%% (2) Adjust the correlation matrix to look like the putput of the
%%% "pdist" function
% t = ones(size(corrMat)); t = triu(t); t = t * 10000;
% p = find(t == 10000); clear t;
% corrDist = abs(corrMat); 
% corrDist = corrMat;  clear corrMat;
% corrDist(p) = []; clear p;
% corrDist = 1 - corrDist;
% % 
% % %%% (3) Define Clustering Parameters
% T = 1.51; % cutoff threshold
% linkType = 'complete'; % linkage type
% clusterTree = linkage(corrDist, linkType);
% clusters = cluster(clusterTree, 'cutoff', T, 'criterion', 'distance');
% save('results\hClustResults_NEW\exp1\allDonorsAvgCorr_clusters.mat', 'clusters');
% f1 = figure;
% [H, Leafs, P] = dendrogram(clusterTree, 13563, 'colorthreshold', T, 'orientation', 'left');
% title('All Donors Avg Correlation', 'fontweight', 'bold');
% xlabel([linkType ' Linkage, Threshold: ' num2str(T)], 'fontweight', 'bold');
% saveas(f1, 'results\hClustResults_NEW\exp1\allDonorsAvgCorr_dendogram.fig');
% saveas(f1, 'results\hClustResults_NEW\exp1\allDonorsAvgCorr_dendogram.jpg');
% save('results\hClustResults_NEW\exp1\allDonorsAvgCorr_genePerm.mat', 'P');
% clear clusters; clear clusterTree; clear corrDist; clear P; clear Leafs; 

% %%% (4) plot heatmaps of the gene expression of each donor. Genes ordered
% %%% according to the dendogram
% load('files/donorsExpMat_5RPKM.mat');
% donorsExpMat_5RPKM = normalizeExpMat(donorsExpMat_5RPKM);
% donorsExpMat_5RPKM = log2(donorsExpMat_5RPKM + (rand(size(donorsExpMat_5RPKM))*(10^-5)));
% %%% Normalize the exp Mt for visualization
% for d = 1 : size(donorsExpMat_5RPKM,3)
%     NexpMat(:,:,d) = bsxfun(@minus, donorsExpMat_5RPKM(:,:,d), mean(donorsExpMat_5RPKM(:,:,d),2));
%     NexpMat(:,:,d) = bsxfun(@rdivide, NexpMat(:,:,d), std(donorsExpMat_5RPKM(:,:,d),[],2));
% end
% load('results\hClustResults_NEW\exp1\allDonorsAvgCorr_clusters.mat');
% load('results\hClustResults_NEW\exp1\allDonorsAvgCorr_genePerm.mat');
% ordP = P(1,end:-1:1);
% for k = 1 : 30
%     f = figure;%('Visible', 'Off');  
%     imagesc(NexpMat(ordP,orderedStrucInd,k)), colormap(redbluecmap); colorbar;
%     title(['Donor ' num2str(k)], 'fontweight', 'bold');
%     set(gca, 'XTick', 1:length(ordS));
%     set(gca, 'XTickLabel', ordS);
%     xlabel('Brain Structures', 'fontweight', 'bold');
%     set(gca, 'YTick', 1:length(gNames));
%     set(gca, 'YTickLabel', gNames(ordP));
%     ylabel('Genes', 'fontweight', 'bold');
%     saveas(f, ['results\hClustResults_NEW\exp1\donorsHM\donor' num2str(k) '_HM.jpg']);
%     saveas(f, ['results\hClustResults_NEW\exp1\donorsHM\donor' num2str(k) '_HM.fig']);
%     close(f);
% end


% %%% (44) plot heatmaps of the gene expression of each stage. Genes ordered
% %%% according to the dendogram
% load('files/donorsExpMat_5RPKM.mat');
% donorsExpMat_5RPKM = normalizeExpMat(donorsExpMat_5RPKM);
% donorsExpMat_5RPKM = log2(donorsExpMat_5RPKM + (rand(size(donorsExpMat_5RPKM))*(10^-5)));
% %%% Normalize the exp Mt for visualization
% for d = 1 : size(donorsExpMat_5RPKM,3)
%     NexpMat(:,:,d) = bsxfun(@minus, donorsExpMat_5RPKM(:,:,d), mean(donorsExpMat_5RPKM(:,:,d),2));
%     NexpMat(:,:,d) = bsxfun(@rdivide, NexpMat(:,:,d), std(donorsExpMat_5RPKM(:,:,d),[],2));
% end
% load('results\hClustResults_NEW\exp1\allDonorsAvgCorr_clusters.mat');
% load('results\hClustResults_NEW\exp1\allDonorsAvgCorr_genePerm.mat');
% ordP = P(1,end:-1:1);
% for k = 1 : length(group)
%     ds = group{k};
%     f = figure;%('Visible', 'Off');  
%     imagesc(mean(NexpMat(ordP,orderedStrucInd,ds), 3)), 
%     colormap(redbluecmap); colorbar;
%     title(['Stage ' num2str(k)], 'fontweight', 'bold');
%     set(gca, 'XTick', 1:length(ordS));
%     set(gca, 'XTickLabel', ordS);
%     xlabel('Brain Structures', 'fontweight', 'bold');
%     set(gca, 'YTick', 1:length(gNames));
%     set(gca, 'YTickLabel', gNames(ordP));
%     ylabel('Genes', 'fontweight', 'bold');
%     saveas(f, ['results\hClustResults_NEW\exp1\stagesHM\stage' num2str(k) '_HM.jpg']);
%     saveas(f, ['results\hClustResults_NEW\exp1\stagesHM\stage' num2str(k) '_HM.fig']);
%     close(f);
% end


%%% (5) plot the 1st PC of each cluster for all donors
% load('files/donorsExpMat_5RPKM.mat');
% donorsExpMat_5RPKM = normalizeExpMat(donorsExpMat_5RPKM);
% donorsExpMat_5RPKM = log2(donorsExpMat_5RPKM + (rand(size(donorsExpMat_5RPKM))*(10^-5)));
% load('results\hClustResults_NEW\exp1\allDonorsAvgCorr_clusters.mat');
% for k = 1 : length(unique(clusters))
%     clustInd = find(clusters == k);
%     for i = 1 : 30
%         [pc, zscores, pcvars] = princomp(donorsExpMat_5RPKM(clustInd,orderedStrucInd,i));
%         pcMat(i,:) = pc(:,1)';
%     end
%     f1 = figure('Visible', 'Off');
%     imagesc(pcMat), colormap(redbluecmap), colorbar
%     title(['Cluster#' num2str(k) ' - ' num2str(length(clustInd))], 'fontweight', 'bold');
%     set(gca, 'XTick', 1:length(ordS));
%     set(gca, 'XTickLabel', ordS);
%     xlabel('Brain Structures', 'fontweight', 'bold');
%     set(gca, 'YTick', 1:30);
%     ylabel('Donors', 'fontweight', 'bold');
%     saveas(f1, ['results\hClustResults_NEW\exp1\clustersPC1\cluster' num2str(k) '_pc1.jpg']);
%     saveas(f1, ['results\hClustResults_NEW\exp1\clustersPC1\cluster' num2str(k) '_pc1.fig']);
%     close(f1);
% end
    
%%% (6) plot the 1st PC of each cluster for all donors
% load('files/donorsExpMat_5RPKM.mat');
% donorsExpMat_5RPKM = normalizeExpMat(donorsExpMat_5RPKM);
% donorsExpMat_5RPKM = log2(donorsExpMat_5RPKM + (rand(size(donorsExpMat_5RPKM))*(10^-5)));
% newExpMat = donorsExpMat_5RPKM(:,nonCs,:);
% newExpMat(:,size(newExpMat,2)+1,:) = mean(donorsExpMat_5RPKM(:,Cs,:),2);
% clear donorsExpMat_5RPKM;
% load('results\hClustResults_NEW\exp1\allDonorsAvgCorr_clusters.mat');
% for k = 1 : length(unique(clusters))
%     clustInd = find(clusters == k);
%     for i = 1 : 30
%         [pc, zscores, pcvars] = princomp(newExpMat(clustInd,:,i));
%         pcMat(i,:) = pc(:,1)';
%     end
%     newPC(1,:) = mean(pcMat(1:4,:));
%     newPC(2,:) = mean(pcMat(5:8,:));
%     newPC(3,:) = mean(pcMat(9:12,:));
%     newPC(4,:) = mean(pcMat(13:16,:));
%     newPC(5,:) = mean(pcMat(17:20,:));
%     newPC(6,:) = mean(pcMat(21:24,:));
%     newPC(7,:) = mean(pcMat(25:30,:));
%     f1 = figure('Visible', 'Off');
%     plot(newPC, 'linewidth', 2); grid on
%     title(['Cluster#' num2str(k) ' - ' num2str(length(clustInd))], 'fontweight', 'bold');
%     set(gca, 'XTick', 1:30);
%     xlabel('Donors', 'fontweight', 'bold');
%     ylabel('PC1', 'fontweight', 'bold');
%     legend(newS);
%     axis([0 8 0.25 0.6])
%     saveas(f1, ['results\hClustResults_NEW\exp1\clustersPC1_str\cluster' num2str(k) '_pc1.jpg']);
%     saveas(f1, ['results\hClustResults_NEW\exp1\clustersPC1_str\cluster' num2str(k) '_pc1.fig']);
%     close(f1);
% end


%%% (7) test for clusters enrichment of ASDs/SZ genes
% load('results\hClustResults_NEW\exp1\allDonorsAvgCorr_clusters.mat');
% dNameList = {'ASDs', 'SZ', 'ND'};
% outFolder = 'results\hClustResults_NEW\exp1\clustersEnrich\';
% enrichAnalyze(clusters, dNameList, outFolder)
%%% -----------------------------------------------------------------------


%%% -----------------------------------------------------------------------
%%% Experiment #2: Cluster genes based on the correlation of all smaples (
%%% Structures & donors together
%%% -----------------------------------------------------------------------
%%% (1) load and rearrange the expression matrix and cluster
load([filesFolder 'files/donorsExpMat_5RPKM.mat']);
newMat_2 = donorsExpMat_5RPKM(:,Cs,:);
newMat_2 = reshape(newMat_2, size(newMat_2,1), size(newMat_2,2)*size(newMat_2,3));
donorsExpMat_5RPKM = permute(donorsExpMat_5RPKM, [1 3 2]);
newMat_1 = donorsExpMat_5RPKM(:,:,nonCs);
newMat_1 = reshape(newMat_1, size(newMat_1,1), size(newMat_1,2)*size(newMat_1,3));
newExpMat = [newMat_1 newMat_2];
newExpMat = normalizeExpMat(newExpMat);
newExpMat = log2(newExpMat + (rand(size(newExpMat))*(10^-5)));
clear newMat_1; clear newMat_2; clear donorsExpMat_5RPKM;
% save('files/finalWGMINmat.mat', 'newExpMat');
T = 1.4; % cutoff threshold
corrType = 'spearman';
linkType = 'complete'; % linkage type
load([resultsFolder '\hClustResults_NEW\exp2\allSamples_corrDist.mat']);
% corrDist = pdist(newExpMat, corrType);
clusterTree = linkage(corrDist, linkType);
clusters = cluster(clusterTree, 'cutoff', T, 'criterion', 'distance');
% clusters = cluster(clusterTree, 'maxclust', 32, 'criterion', 'distance');
% save('results\hClustResults_NEW\exp2\allSamples_clusters.mat', 'clusters');
% save('results\hClustResults_NEW\exp2\allSamples_corrDist.mat', 'corrDist');
f1 = figure;
[H, Leafs, P] = dendrogram(clusterTree, 13563, 'colorthreshold', T, 'orientation', 'left');
title('All Samples', 'fontweight', 'bold');
xlabel([linkType ' Linkage, Threshold: ' num2str(T)], 'fontweight', 'bold');
saveas(f1, 'results\hClustResults_NEW\exp2\allSamples_dendogram.fig');
saveas(f1, 'results\hClustResults_NEW\exp2\allSamples_dendogram.jpg');
save('results\hClustResults_NEW\exp2\allSamples_genePerm.mat', 'P');
% clear clusters; clear clusterTree; clear corrDist; clear P; clear Leafs;


%%% (2) Plot the HM of the reordered matrix
% % load('files/donorsExpMat_5RPKM.mat');
% % newMat_2 = donorsExpMat_5RPKM(:,Cs,:);
% % newMat_2 = reshape(newMat_2, size(newMat_2,1), size(newMat_2,2)*size(newMat_2,3));
% % donorsExpMat_5RPKM = permute(donorsExpMat_5RPKM, [1 3 2]);
% % newMat_1 = donorsExpMat_5RPKM(:,:,nonCs);
% % newMat_1 = reshape(newMat_1, size(newMat_1,1), size(newMat_1,2)*size(newMat_1,3));
% % newExpMat = [newMat_1 newMat_2];
% % newExpMat = normalizeExpMat(newExpMat);
% % newExpMat = log2(newExpMat + (rand(size(newExpMat))*(10^-5)));
% % clear newMat_1; clear newMat_2; clear donorsExpMat_5RPKM;
% NexpMat = bsxfun(@minus, newExpMat, mean(newExpMat,2));
% NexpMat = bsxfun(@rdivide, NexpMat, std(newExpMat,[],2));
% load('results\hClustResults_NEW\exp2\allSamples_genePerm.mat');
% ordP = P(1,end:-1:1);
% f = figure;%('Visible', 'Off');  
% HeatMap(NexpMat(ordP,:), 'ColorMap', 'redbluecmap', 'DisplayRange', 3);
% title('All Samples', 'fontweight', 'bold');
% set(gca, 'XTick', 1:length(ordS));
% set(gca, 'XTickLabel', ordS);
% xlabel('Brain Structures', 'fontweight', 'bold');
% set(gca, 'YTick', 1:length(gNames));
% set(gca, 'YTickLabel', gNames(ordP));
% ylabel('Genes', 'fontweight', 'bold');
% saveas(f, ['results\hClustResults_NEW\exp2\allSamplesHM.jpg']);
% saveas(f, ['results\hClustResults_NEW\exp2\allSamplesHM.fig']);
% close(f);


%%% (3) Clustering using clustergram to get the heatmap
% load('files/donorsExpMat_5RPKM.mat');
% newMat_2 = donorsExpMat_5RPKM(:,Cs,:);
% newMat_2 = reshape(newMat_2, size(newMat_2,1), size(newMat_2,2)*size(newMat_2,3));
% donorsExpMat_5RPKM = permute(donorsExpMat_5RPKM, [1 3 2]);
% newMat_1 = donorsExpMat_5RPKM(:,:,nonCs);
% newMat_1 = reshape(newMat_1, size(newMat_1,1), size(newMat_1,2)*size(newMat_1,3));
% newExpMat = [newMat_1 newMat_2];
% newExpMat = normalizeExpMat(newExpMat);
% newExpMat = log2(newExpMat + (rand(size(newExpMat))*(10^-5)));
% clear newMat_1; clear newMat_2; clear donorsExpMat_5RPKM;
% T = 1.4; % cutoff threshold
% corrType = 'spearman';
% linkType = 'complete'; % linkage type
% C = clustergram(newExpMat, 'RowLabels', gNames, ...
%         'Standardize', 'none', 'Cluster', 1, 'RowPDist', corrType, ...
%         'Linkage', linkType, 'Dendrogram', T, 'Colormap', redbluecmap);
% h = plot(C);
% saveas(h, 'results\hClustResults_NEW\exp2\allSamples_Clustergram.fig');
% saveas(h, 'results\hClustResults_NEW\exp2\allSamples_Clustergram.jpg');
% close all; clear C;


%%% (4) Plot the 1st PC of each cluster
% load('files/donorsExpMat_5RPKM.mat');
% donorsExpMat_5RPKM = normalizeExpMat(donorsExpMat_5RPKM);
% donorsExpMat_5RPKM = log2(donorsExpMat_5RPKM + (rand(size(donorsExpMat_5RPKM))*(10^-5)));
% load('results\hClustResults_NEW\exp2\allSamples_clusters.mat');
% for k = 1 : length(unique(clusters))
%     clustInd = find(clusters == k);
%     for i = 1 : 30
%         [pc, zscores, pcvars] = princomp(donorsExpMat_5RPKM(clustInd,orderedStrucInd,i));
%         pcMat(i,:) = pc(:,1)';
%     end
%     f1 = figure('Visible', 'Off');
%     imagesc(pcMat), colormap(redbluecmap), colorbar
%     title(['Cluster#' num2str(k) ' - ' num2str(length(clustInd))], 'fontweight', 'bold');
%     set(gca, 'XTick', 1:length(ordS));
%     set(gca, 'XTickLabel', ordS);
%     xlabel('Brain Structures', 'fontweight', 'bold');
%     set(gca, 'YTick', 1:30);
%     ylabel('Donors', 'fontweight', 'bold');
%     saveas(f1, ['results\hClustResults_NEW\exp2\clustersPC1\cluster' num2str(k) '_pc1.jpg']);
%     saveas(f1, ['results\hClustResults_NEW\exp2\clustersPC1\cluster' num2str(k) '_pc1.fig']);
%     close(f1);
% end


%%% (5) Plot the 1st PC of each cluster
% load('files/donorsExpMat_5RPKM.mat');
% donorsExpMat_5RPKM = normalizeExpMat(donorsExpMat_5RPKM);
% donorsExpMat_5RPKM = log2(donorsExpMat_5RPKM + (rand(size(donorsExpMat_5RPKM))*(10^-5)));
% newExpMat = donorsExpMat_5RPKM(:,nonCs,:);
% newExpMat(:,size(newExpMat,2)+1,:) = mean(donorsExpMat_5RPKM(:,Cs,:),2);
% clear donorsExpMat_5RPKM;
% % % % newExpMat = bsxfun(@minus, newExpMat, mean(newExpMat,2));
% % % % newExpMat = bsxfun(@rdivide, newExpMat, std(newExpMat,[],2));
% load('results\hClustResults_NEW\exp2\allSamples_clusters.mat');
% for k = 1 : length(unique(clusters))
%     clustInd = find(clusters == k);
%     for i = 1 : 30
%         [pc, zscores, pcvars] = princomp(newExpMat(clustInd,:,i));
%         pcMat(i,:) = pc(:,1)';
%     end
%     newPC(1,:) = mean(pcMat(1:4,:));
%     newPC(2,:) = mean(pcMat(5:8,:));
%     newPC(3,:) = mean(pcMat(9:12,:));
%     newPC(4,:) = mean(pcMat(13:16,:));
%     newPC(5,:) = mean(pcMat(17:20,:));
%     newPC(6,:) = mean(pcMat(21:24,:));
%     newPC(7,:) = mean(pcMat(25:30,:));
%     f1 = figure('Visible', 'Off');
%     plot(newPC, 'linewidth', 2); grid on
%     title(['Cluster#' num2str(k) ' - ' num2str(length(clustInd))], 'fontweight', 'bold');
%     set(gca, 'XTick', 1:30);
%     xlabel('Donors', 'fontweight', 'bold');
%     ylabel('PC1', 'fontweight', 'bold');
%     legend(newS);
% %     axis([0 8 0.15 0.6])
%     saveas(f1, ['results\hClustResults_NEW\exp2\clustersPC1_str\cluster' num2str(k) '_pc1.jpg']);
%     saveas(f1, ['results\hClustResults_NEW\exp2\clustersPC1_str\cluster' num2str(k) '_pc1.fig']);
%     close(f1);
% end


%%% (6) Plot the heat map of each cluster
% load('files/donorsExpMat_5RPKM.mat');
% newMat_2 = donorsExpMat_5RPKM(:,Cs,:);
% newMat_2 = reshape(newMat_2, size(newMat_2,1), size(newMat_2,2)*size(newMat_2,3));
% donorsExpMat_5RPKM = permute(donorsExpMat_5RPKM, [1 3 2]);
% newMat_1 = donorsExpMat_5RPKM(:,:,nonCs);
% newMat_1 = reshape(newMat_1, size(newMat_1,1), size(newMat_1,2)*size(newMat_1,3));
% newExpMat = [newMat_1 newMat_2];
% newExpMat = normalizeExpMat(newExpMat);
% newExpMat = log2(newExpMat + (rand(size(newExpMat))*(10^-5)));
% clear newMat_1; clear newMat_2; clear donorsExpMat_5RPKM;
% corrType = 'spearman';
% linkType = 'complete'; % linkage type
% load('results\hClustResults_NEW\exp2\allSamples_clusters.mat');
% for k = 1 : length(unique(clusters))
%     clustInd = find(clusters == k);
%     C = clustergram(newExpMat(clustInd,:), 'RowLabels', gNames(clustInd), ...
%         'Standardize', 'none', 'Cluster', 1, 'RowPDist', corrType, ...
%         'Linkage', linkType, 'Colormap', redbluecmap);
%     title = ['Cluster#' num2str(k) ' - ' num2str(length(clustInd)) ' genes'];
%     addTitle(C, title);
%     h = plot(C); 
%     saveas(h, ['results\hClustResults_NEW\exp2\clustersHM\cluster' num2str(k) '_HM.fig']);
%     saveas(h, ['results\hClustResults_NEW\exp2\clustersHM\cluster' num2str(k) '_HM.jpg']);
%     close all; clear C; clear clustInd;
% end


%%% (7) Plot the avg expression pattern of each cluster
% load('files/donorsExpMat_5RPKM.mat');
% donorsExpMat_5RPKM = normalizeExpMat(donorsExpMat_5RPKM);
% donorsExpMat_5RPKM = log2(donorsExpMat_5RPKM + (rand(size(donorsExpMat_5RPKM))*(10^-5)));
% newExpMat = donorsExpMat_5RPKM(:,nonCs,:);
% newExpMat(:,size(newExpMat,2)+1,:) = mean(donorsExpMat_5RPKM(:,Cs,:),2);
% clear donorsExpMat_5RPKM;
% load('results\hClustResults_NEW\exp2\allSamples_clusters.mat');
% for k = 1 : length(unique(clusters))
%     clear cp
%     clustInd = find(clusters == k);
%     cp(:,:) = median(newExpMat(clustInd,:,:), 1);
%     newCP(:,1) = mean(cp(:,1:4),2);
%     newCP(:,2) = mean(cp(:,5:8),2);
%     newCP(:,3) = mean(cp(:,9:12),2);
%     newCP(:,4) = mean(cp(:,13:16),2);
%     newCP(:,5) = mean(cp(:,17:20),2);
%     newCP(:,6) = mean(cp(:,21:24),2);
%     newCP(:,7) = mean(cp(:,25:30),2);
%     f1 = figure('Visible', 'Off');
%     plot(newCP', 'linewidth', 3); grid on
%     title({['Module ' num2str(k)], ['Number of Genes = ' num2str(length(clustInd))]}, 'fontweight', 'bold');
%     set(gca, 'XTick', 0:8);
%     xlabel('Developmental Stages', 'fontweight', 'bold');
%     ylabel('Avgerage Expression', 'fontweight', 'bold');
%     legend(newS);
%     axis([0 8 min(min(newCP))-0.2 max(max(newCP))+0.2])
%     saveas(f1, ['results\hClustResults_NEW\exp2\clustersPattern\cluster' num2str(k) '_pc1.jpg']);
%     saveas(f1, ['results\hClustResults_NEW\exp2\clustersPattern\cluster' num2str(k) '_pc1.fig']);
%     close(f1);
% end


%%% (8) test for enrichment
% load('results\hClustResults_NEW\exp2\allSamples_clusters.mat');
% dNameList = {'Parkinsons'};
% outFolder = 'results\hClustResults_NEW\exp2\clustersEnrich\';
% % enrichAnalyze(clusters, dNameList, outFolder)
% % enrichAnalyze_MOD(clusters, dNameList, outFolder)
% % [num txt] = xlsread(['Data\newASDgenes.xls']);
% [num txt] = xlsread(['Data\Donors\Parkinsons.xlsx']);
% list = txt(2:end, 1);
% listName = txt{1};
% clear num; clear txt;
% enrichAnalyze2(clusters, list, listName, outFolder);

%%% (9) test for cell-type enrichment
% load('results\hClustResults_NEW\exp2\allSamples_clusters.mat');
% ctNameList = {'Neurons', 'Oligodendrytes', 'Astrocytes'};
% outFolder = 'results\hClustResults_NEW\exp2\clustersEnrich\';
% enrichAnalyze_cellType(clusters, ctNameList, outFolder);

%%% (4Dec2013) test for enrichment
load([resultsFolder 'hClustResults_NEW\exp2\allSamples_clusters.mat']);
outFolder = [resultsFolder 'hClustResults_NEW\'];
markersFile = [filesFolder 'files\DataFiles\cellTypeMarkers.xlsx'];
enrichAnalyze_4Dec2013(clusters, markersFile, outFolder, filesFolder, resultsFolder);

%%% (10) export to cytoscape for visualization
% load('files/donorsExpMat_5RPKM.mat');
% newMat_2 = donorsExpMat_5RPKM(:,Cs,:);
% newMat_2 = reshape(newMat_2, size(newMat_2,1), size(newMat_2,2)*size(newMat_2,3));
% donorsExpMat_5RPKM = permute(donorsExpMat_5RPKM, [1 3 2]);
% newMat_1 = donorsExpMat_5RPKM(:,:,nonCs);
% newMat_1 = reshape(newMat_1, size(newMat_1,1), size(newMat_1,2)*size(newMat_1,3));
% newExpMat = [newMat_1 newMat_2];
% newExpMat = normalizeExpMat(newExpMat);
% newExpMat = log2(newExpMat + (rand(size(newExpMat))*(10^-5)));
% clear donorsExpMat_5RPKM; clear newMat_2; clear newMat_1;
% corrType = 'spearman';
% corrDist = pdist(newExpMat, corrType);
% clear newExpMat;
% fcm = squareform(corrDist);
% save('results\hClustResults_NEW\exp2\cytoscape\finalCorrMat.mat', 'fcm');

% load('results\hClustResults_NEW\exp2\allSamples_clusters.mat');
% load('results\hClustResults_NEW\exp2\cytoscape\finalCorrMat.mat');
% for j = 1 : max(clusters)
%     clustInd = find(clusters == j);
%     clustCorrMat = fcm(clustInd, clustInd);
% %     [num txt] = xlsread('results\hClustResults_NEW\exp2\clustersEnrich\clustInfo.xls', j);
% %     clustGeneNames = txt(2:end,1);
% %     clear num; clear txt;
%     outF = ['results\hClustResults_NEW\exp2\cytoscape\clust' num2str(j) '.csv'];
%     csvwrite(outF, clustCorrMat, 1, 1);
% end


%%% (11) export to cytoscape for visualization II
% load('files/donorsExpMat_5RPKM.mat');
% newMat_2 = donorsExpMat_5RPKM(:,Cs,:);
% newMat_2 = reshape(newMat_2, size(newMat_2,1), size(newMat_2,2)*size(newMat_2,3));
% donorsExpMat_5RPKM = permute(donorsExpMat_5RPKM, [1 3 2]);
% newMat_1 = donorsExpMat_5RPKM(:,:,nonCs);
% newMat_1 = reshape(newMat_1, size(newMat_1,1), size(newMat_1,2)*size(newMat_1,3));
% newExpMat = [newMat_1 newMat_2];
% newExpMat = normalizeExpMat(newExpMat);
% newExpMat = log2(newExpMat + (rand(size(newExpMat))*(10^-5)));
% clear donorsExpMat_5RPKM; clear newMat_2; clear newMat_1;
% corrType = 'spearman';

% load('C:\Users\amahfouz\Documents\MATLAB\ASD_study_new_BACKUP\results\hClustResults_NEW\exp2\allSamples_corrDist.mat');
% corrDist2 = squareform(corrDist);
% clear corrDist;
% 
% load('results\hClustResults_NEW\exp2\allSamples_clusters.mat');
% for j = 1 : max(clusters)
%     clear clustInd; clear clustPDist; clear xN1; clear xN2; clear clustGeneNames;
%     clustInd = find(clusters == j);
%     clustPDist = corrDist2(clustInd, clustInd);
% %     clustPDist = pdist(newExpMat(clustInd,:), corrType);
%     [num txt] = xlsread('results\hClustResults_NEW\exp2\clustersEnrich\clustInfo.xls', j);
%     clustGeneNames = txt(2:end,1);
%     clear num; clear txt;
%     count = 0;
%     for g1 = 1 : length(clustGeneNames)
%         for g2 = (g1+1) : length(clustGeneNames)
%             count = count + 1;
%             xN1{count} = clustGeneNames{g1};
%             xN2{count} = clustGeneNames{g2};
%         end
%     end
%     outF = ['results\hClustResults_NEW\exp2\cytoscape\clust' num2str(j) '.csv'];
%     outArr = cell(length(clustPDist), 3);
%     outArr(:,1) = xN1';
%     outArr(:,2) = xN2';
%     outArr(:,3) = num2cell(clustPDist)';
% 
%     TEMP(j) = length(find(clustPDist <= 0.1));
    
%     fid = fopen(outF, 'w');
%     for r = 1:size(outArr, 1)
%         fprintf(fid, '%s', outArr{r, 1});
%         fprintf(fid, ',');
%         fprintf(fid, '%s', outArr{r, 2});
%         fprintf(fid, ',');
%         fprintf(fid, '%f', outArr{r, 3});
%         fprintf(fid, '\n');
%     end
%     fclose(fid);
% end
%%% -----------------------------------------------------------------------






%%% -----------------------------------------------------------------------
%%% Experiment #3: Cluster genes based on the correlation of all smaples
%%% (within each age stage separately)
%%% -----------------------------------------------------------------------
%%% (1) load and rearrange the expression matrix and cluster
% load('files/donorsExpMat_5RPKM.mat');
% donorsExpMat_5RPKM = normalizeExpMat(donorsExpMat_5RPKM);
% donorsExpMat_5RPKM = log2(donorsExpMat_5RPKM + (rand(size(donorsExpMat_5RPKM))*(10^-5)));
% for i = 1: length(group)
%     load(['results/geneCorr/Avg(log2+epsilon)/ageGroup' num2str(i) '_corrMat.mat']);
%     t = ones(size(corrMat)); t = triu(t); t = t * 10000;
%     p = find(t == 10000); clear t;
%     corrDist = abs(corrMat); 
%     corrDist = corrMat;  clear corrMat;
%     corrDist(p) = []; clear p;
%     corrDist = 1 - corrDist;
%     T = 1.5; % cutoff threshold
%     corrType = 'spearman';
%     linkType = 'complete'; % linkage type
%     clusterTree = linkage(corrDist, linkType);
%     clusters = cluster(clusterTree, 'cutoff', T, 'criterion', 'distance');
%     save(['results\hClustResults_NEW\exp3\ageGroup ' num2str(i) '_clusters.mat'], 'clusters');
%     clear corrDist; clear clusterTree; clear clusters;
% end


%%% (2) test for enrichment
% for i = 1 : length(group)
%     load(['results\hClustResults_NEW\exp3\ageGroup ' num2str(i) '_clusters.mat']);
%     dNameList = {'ASDs', 'SZ', 'ND'};
%     outFolder = 'results\hClustResults_NEW\exp3\clustersEnrich\';
%     if ~exist(outFolder)
%         mkdir(outFolder);
%     end
%     enrichAnalyze(clusters, dNameList, outFolder)
% end


%%% (3) Plot the heat map of each cluster
% load('files/donorsExpMat_5RPKM.mat');
% donorsExpMat_5RPKM = normalizeExpMat(donorsExpMat_5RPKM);
% donorsExpMat_5RPKM = log2(donorsExpMat_5RPKM + (rand(size(donorsExpMat_5RPKM))*(10^-5)));
% newMat_2 = donorsExpMat_5RPKM(:,Cs,:);
% newMat_2 = reshape(newMat_2, size(newMat_2,1), size(newMat_2,2)*size(newMat_2,3));
% donorsExpMat_5RPKM = permute(donorsExpMat_5RPKM, [1 3 2]);
% newMat_1 = donorsExpMat_5RPKM(:,:,nonCs);
% newMat_1 = reshape(newMat_1, size(newMat_1,1), size(newMat_1,2)*size(newMat_1,3));
% newExpMat = [newMat_1 newMat_2];
% T = 1.5;
% corrType = 'spearman';
% linkType = 'complete'; % linkage type
% for i = 1 : length(group)
%     clear newMat_1; clear newMat_2;
%     load(['results\hClustResults_NEW\exp3\ageGroup ' num2str(i) '_clusters.mat']);
%     for k = 1 : length(unique(clusters))
%         clustInd = find(clusters == k);
%         C = clustergram(newExpMat(clustInd,:), 'RowLabels', gNames(clustInd), ...
%             'Standardize', 'none', 'Cluster', 1, 'RowPDist', corrType, ...
%             'Linkage', linkType, 'Colormap', redbluecmap);
%         title = ['Cluster#' num2str(k) ' - ' num2str(length(clustInd)) ' genes'];
%         addTitle(C, title);
%         h = plot(C); 
%         saveas(h, ['results\hClustResults_NEW\exp3\HM\ageGroup' num2str(i) '_HM.fig']);
%         saveas(h, ['results\hClustResults_NEW\exp3\HM\ageGroup' num2str(i) '_HM.jpg']);
%         close all; clear C; clear clustInd;
%     end
%     clear newExpMat;
% end




%%% -----------------------------------------------------------------------
%%% Experiment #4: clustering of ASDs and SZ genes
%%% -----------------------------------------------------------------------
% load('files/asdGeneIDs_5RPKM.mat');
% load('files/szGeneIDs_5RPKM.mat');
% load('files/donorsExpMat_5RPKM.mat');
% donorsExpMat_5RPKM = normalizeExpMat(donorsExpMat_5RPKM);
% donorsExpMat_5RPKM = log2(donorsExpMat_5RPKM + (rand(size(donorsExpMat_5RPKM))*(10^-5)));
% newMat_2 = donorsExpMat_5RPKM(:,Cs,:);
% newMat_2 = reshape(newMat_2, size(newMat_2,1), size(newMat_2,2)*size(newMat_2,3));
% donorsExpMat_5RPKM = permute(donorsExpMat_5RPKM, [1 3 2]);
% newMat_1 = donorsExpMat_5RPKM(:,:,nonCs);
% newMat_1 = reshape(newMat_1, size(newMat_1,1), size(newMat_1,2)*size(newMat_1,3));
% newExpMat = [newMat_1 newMat_2];
% 
% T = 1.3;
% corrType = 'spearman';
% linkType = 'complete'; % linkage type

%%%(1)ASDs
% corrDist = pdist(newExpMat(asdIDs_5RPKM,:), corrType);
% clusterTree = linkage(corrDist, linkType);
% clusters = cluster(clusterTree, 'cutoff', T, 'criterion', 'distance');
% save(['results\hClustResults_NEW\exp4\ASD_clusters.mat'], 'clusters');
% clear corrDist; clear clusterTree; clear clusters;

% C = clustergram(newExpMat(asdIDs_5RPKM,:), 'RowLabels', gNames(asdIDs_5RPKM), ...
%     'Standardize', 'none', 'Cluster', 1, 'RowPDist', corrType, ...
%     'Linkage', linkType, 'Dendrogram', T, 'Colormap', redbluecmap);
% title = 'ASDs-related Genes';
% addTitle(C, title);
% h = plot(C); 
% saveas(h, ['results\hClustResults_NEW\exp4\ASDs']);
% saveas(h, ['results\hClustResults_NEW\exp4\ASDs']);
% close all; clear C; clear clustInd;

% load('results\hClustResults_NEW\exp4\ASD_clusters.mat');
% for k = 1 : length(unique(clusters))
%     clustInd = find(clusters == k);
%     C = clustergram(newExpMat(clustInd,:), 'RowLabels', gNames(clustInd), ...
%         'Standardize', 'none', 'Cluster', 1, 'RowPDist', corrType, ...
%         'Linkage', linkType, 'Colormap', redbluecmap);
%     title = ['Cluster#' num2str(k) ' - ' num2str(length(clustInd)) ' genes'];
%     addTitle(C, title);
%     h = plot(C); 
%     saveas(h, ['results\hClustResults_NEW\exp4\ASDclustersHM\cluster' num2str(k) '_HM.fig']);
%     saveas(h, ['results\hClustResults_NEW\exp4\ASDclustersHM\cluster' num2str(k) '_HM.jpg']);
%     close all; clear C; clear clustInd;
% end

% load('results\hClustResults_NEW\exp4\ASD_clusters.mat');
% load('files/donorsExpMat_5RPKM.mat');
% newExpMat2 = donorsExpMat_5RPKM(:,nonCs,:);
% newExpMat2(:,size(newExpMat2,2)+1,:) = mean(donorsExpMat_5RPKM(:,Cs,:),2);
% for k = 1 : length(unique(clusters))
%     clear cp
%     clustInd = find(clusters == k);
%     cp(:,:) = median(newExpMat2(clustInd,:,:), 1);
%     newCP(:,1) = mean(cp(:,1:4),2);
%     newCP(:,2) = mean(cp(:,5:8),2);
%     newCP(:,3) = mean(cp(:,9:12),2);
%     newCP(:,4) = mean(cp(:,13:16),2);
%     newCP(:,5) = mean(cp(:,17:20),2);
%     newCP(:,6) = mean(cp(:,21:24),2);
%     newCP(:,7) = mean(cp(:,25:30),2);
%     f1 = figure;%('Visible', 'Off');
%     plot(newCP', 'linewidth', 2); grid on
%     title(['Cluster#' num2str(k) ' - ' num2str(length(clustInd))], 'fontweight', 'bold');
%     set(gca, 'XTick', 0:8);
%     xlabel('Age Stages', 'fontweight', 'bold');
%     ylabel('Avg. Exp.', 'fontweight', 'bold');
%     legend(newS);
%     axis([0 8 min(min(newCP))-0.2 max(max(newCP))+0.2])
%     saveas(f1, ['results\hClustResults_NEW\exp4\ASDclustersP\cluster' num2str(k) '_HM.fig']);
%     saveas(f1, ['results\hClustResults_NEW\exp4\ASDclustersP\cluster' num2str(k) '_HM.jpg']);    
%     close(f1);
% end

%%%(2)SZ
% corrDist = pdist(newExpMat(szIDs_5RPKM,:), corrType);
% clusterTree = linkage(corrDist, linkType);
% clusters = cluster(clusterTree, 'cutoff', T, 'criterion', 'distance');
% save(['results\hClustResults_NEW\exp4\SZ_clusters.mat'], 'clusters');
% clear corrDist; clear clusterTree; clear clusters;
% 
% C = clustergram(newExpMat(szIDs_5RPKM,:), 'RowLabels', gNames(szIDs_5RPKM), ...
%     'Standardize', 'none', 'Cluster', 1, 'RowPDist', corrType, ...
%     'Linkage', linkType, 'Dendrogram', T, 'Colormap', redbluecmap);
% title = 'ASDs-related Genes';
% addTitle(C, title);
% h = plot(C); 
% saveas(h, ['results\hClustResults_NEW\exp4\SZ.fig']);
% saveas(h, ['results\hClustResults_NEW\exp4\SZ.jpg']);
% close all; clear C; clear clustInd;
% 
% load('results\hClustResults_NEW\exp4\SZ_clusters.mat');
% for k = 1 : length(unique(clusters))
%     clustInd = find(clusters == k);
%     C = clustergram(newExpMat(clustInd,:), 'RowLabels', gNames(clustInd), ...
%         'Standardize', 'none', 'Cluster', 1, 'RowPDist', corrType, ...
%         'Linkage', linkType, 'Colormap', redbluecmap);
%     title = ['Cluster#' num2str(k) ' - ' num2str(length(clustInd)) ' genes'];
%     addTitle(C, title);
%     h = plot(C); 
%     saveas(h, ['results\hClustResults_NEW\exp4\SZclustersHM\cluster' num2str(k) '_HM.fig']);
%     saveas(h, ['results\hClustResults_NEW\exp4\SZclustersHM\cluster' num2str(k) '_HM.jpg']);
%     close all; clear C; clear clustInd;
% end

% load('results\hClustResults_NEW\exp4\SZ_clusters.mat');
% load('files/donorsExpMat_5RPKM.mat');
% newExpMat2 = donorsExpMat_5RPKM(:,nonCs,:);
% newExpMat2(:,size(newExpMat2,2)+1,:) = mean(donorsExpMat_5RPKM(:,Cs,:),2);
% for k = 1 : length(unique(clusters))
%     clear cp
%     clustInd = find(clusters == k);
%     cp(:,:) = median(newExpMat2(clustInd,:,:), 1);
%     newCP(:,1) = mean(cp(:,1:4),2);
%     newCP(:,2) = mean(cp(:,5:8),2);
%     newCP(:,3) = mean(cp(:,9:12),2);
%     newCP(:,4) = mean(cp(:,13:16),2);
%     newCP(:,5) = mean(cp(:,17:20),2);
%     newCP(:,6) = mean(cp(:,21:24),2);
%     newCP(:,7) = mean(cp(:,25:30),2);
%     f1 = figure;%('Visible', 'Off');
%     plot(newCP', 'linewidth', 2); grid on
%     title(['Cluster#' num2str(k) ' - ' num2str(length(clustInd))], 'fontweight', 'bold');
%     set(gca, 'XTick', 0:8);
%     xlabel('Age Stages', 'fontweight', 'bold');
%     ylabel('Avg. Exp.', 'fontweight', 'bold');
%     legend(newS);
%     axis([0 8 min(min(newCP))-0.2 max(max(newCP))+0.2])
%     saveas(f1, ['results\hClustResults_NEW\exp4\SZclustersP\cluster' num2str(k) '_HM.fig']);
%     saveas(f1, ['results\hClustResults_NEW\exp4\SZclustersP\cluster' num2str(k) '_HM.jpg']);    
%     close(f1);
% end



