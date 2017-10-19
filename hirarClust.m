%%% Hirarchial clustering of genes


clear all;

% filesDir = '/tudelft.net/staff-bulk/ewi/mm/DBL/amahfouz/MatlabFiles/';

%%% Calculate the gene correlation networks for each donor separately
load(['files/donorsExpMat_5RPKM.mat']);
load(['files/gNames.mat']);
load(['files/genesStatus_5RPKM.mat']);
load(['files/strucs.mat']);

gNames(find(genesStatus_5RPKM == 0)) = [];

strucIndInc = [2, 3, 5, 6, 7, 8, 9, 15, 18, 19, 20, 21, 22, 23, 24, 26];
for i = 1 : length(strucIndInc)
    S{i} = strucs{strucIndInc(i)};
end

donorsExpMat_5RPKM = normalizeExpMat(donorsExpMat_5RPKM);
donorsExpMat_5RPKM = log2(donorsExpMat_5RPKM + (rand(size(donorsExpMat_5RPKM))*(10^-5)));

% for i = 1 : size(donorsExpMat_5RPKM, 3)
%     
%     corrDist = pdist(donorsExpMat_5RPKM(:,:,i), 'spearman');
%     clusterTree = linkage(corrDist, 'complete');
% %     clusters = cluster(clusterTree, 'cutoff', 1.7, 'criterion', 'distance');
%     
% %     f = figure('Visible', 'off');
%     [H,T,perm] = dendrogram(clusterTree, 13563, 'colorthreshold', 1.7);  
%     reorderedMat = donorsExpMat_5RPKM(perm,:,i);
%     HeatMap(reorderedMat, 'ColumnLabels', S, 'Colormap', redbluecmap);
%     clear reorderedMat; clear perm; clear H; clear T;
%     
% %     title(['Donor' num2str(i)]);
% %     
% %     save(['results\hClustResults\donors_complete\donor' num2str(i) '_clusters.mat'], 'clusters');
% %     saveas(f, ['results\hClustResults\donors_complete\donor' num2str(i) '_dendogram.fig']);
% %     saveas(f, ['results\hClustResults\donors_complete\donor' num2str(i) '_dendogram.jpg']);
% %     clear corrDist; clear clusterTree; clear clusters;
%     
% %     save(['results\hClustResults\donor' num2str(i) '_clusters.mat'], 'clusters');
% %     clear corrDist; clear clusterTree; clear clusters;
% %     
%     C = clustergram(donorsExpMat_5RPKM(:,:,i), 'RowLabels', gNames, 'ColumnLabels', S, ...
%         'Standardize', 'none', 'Cluster', 1, 'RowPDist', 'spearman', ...
%         'Linkage', 'complete', 'Dendrogram', 1.7, 'Colormap', redbluecmap, ...
%         'DisplayRange', 3);
% % 
% %     h = plot(C);
% %     saveas(h, ['results\hClustResults\donor' num2str(i) '_dendogram.fig']);
% %     saveas(h, ['results\hClustResults\donor' num2str(i) '_dendogram.jpg']);
% %     save(['results/hClustResults/donor' num2str(i) '_C.mat'], 'C');
% %     close all; clear C;
%     
% % %     save([filesDir 'results/hClust/donor' num2str(i) '_clusters.mat'], 'clusters');
% % %     
% % %     h = figure;
% % %     set(h, 'visible', 'off');
% % %     clustergram(donorsExpMat_5RPKM(:,:,i), 'RowLabels', gNames, 'ColumnLabels', S,...
% % %         'Colormap',redbluecmap);
% % %     saveas(h, [filesDir 'results/hClust/donor' num2str(i) '_hClust.fig']);
% % 
% % %     load(['results\hClustResults\hClustTrees\donor' num2str(i) '_clusterTree_AVG.mat']);
% % %     dendrogram(clusterTree, 13563)
% % %     set(H,'LineWidth',2)
% % 
% % %     C = clustergram(donorsExpMat_5RPKM(:,:,i), 'RowLabels', gNames, 'ColumnLabels', S, ...
% % %         'Standardize', 'none', 'Cluster', 1, 'RowPDist', 'spearman', ...
% % %         'Linkage', 'average');%, 'RowGroupMarker', rm);
% %     
% end
% 
% 
% group{1} = [1, 2, 3 ,4];
% group{2} = [5, 6, 7, 8];
% group{3} = [9, 10, 11, 12];
% group{4} = [13, 14, 15, 16];
% group{5} = [17, 18, 19, 20];
% group{6} = [21, 22, 23, 24];
% group{7} = [25, 26, 27, 28, 29, 30];
% 
% % 
% for i = 1 : 7
%         
% %     corrDist = pdist(donorsExpMat_5RPKM(:,:,i), 'spearman');
% %     clusterTree = linkage(corrDist, 'average');    
% %     dendrogram(clusterTree, 13563)
%     
%     load(['results/geneCorr/Avg(log2+epsilon)/ageGroup' num2str(i) '_corrMat.mat']);
% %     
%     T = 1.7;
%     
%     t = ones(size(corrMat));
%     t = triu(t);
%     t = t * 10000;
%     p = find(t == 10000);
%     clear t;
%     
%     corrDist = abs(corrMat); 
%     corrDist = corrMat; 
%     clear corrMat;
%     corrDist(p) = [];
%     clear p;
% 
%     corrDist = 1 - corrDist;
%     
%     clusterTree = linkage(corrDist, 'complete');
%     clusters = cluster(clusterTree, 'cutoff', T, 'criterion', 'distance');
%     
% %     f1 = figure(1);
% %     [H,T,perm] = dendrogram(clusterTree, 13563, 'colorthreshold', T, ...
% %         'orientation', 'left');
% %     title(['Stage ' num2str(i)], 'fontweight', 'bold');
% %     xlabel('Complete Linkage', 'fontweight', 'bold');
% %     
% %     ds = group{i};
% %     for j = 1 : length(ds)
% %         clear reorderedMat; clear HM;
% %         reorderedMat = donorsExpMat_5RPKM(perm,:,ds(j));
% %         HM = HeatMap(reorderedMat, 'RowLabels', gNames(perm), 'ColumnLabels', S, ...
% %             'Colormap', redbluecmap);
% %         f2 = plot(HM);
% %         title(['Donor ' num2str(ds(j))], 'fontweight', 'bold');
% %         saveas(f2, ['results\hClustResults\Stages\stages_complete\donor' num2str(ds(j)) '_heatMap.fig']);
% %         saveas(f2, ['results\hClustResults\Stages\stages_complete\donor' num2str(ds(j)) '_heatMap.jpg']);
% %     end
% %    
% %     saveas(f1, ['results\hClustResults\Stages\stages_complete\ageGroup' num2str(i) '_dendogram.fig']);
% %     saveas(f1, ['results\hClustResults\Stages\stages_complete\ageGroup' num2str(i) '_dendogram.jpg']);
%     save(['results\hClustResults\Stages\stagesModules\stages_complete\ageGroup' num2str(i) '_clusters.mat'], 'clusters');
% %     close(f1); close all;
%     clear clusters; clear clusterTree; clear corrDist;
%     
% end

%%%-------------------------------------------------
% load('Data\Samples_str_age.mat');
% 
% % corrDist = pdist(donorsExpMat_5RPKM, 'spearman');
% % clusterTree = linkage(corrDist, 'average');
% % clusters = cluster(clusterTree, 'cutoff', 0.8, 'criterion', 'distance');
% 
% %     dendrogram(clusterTree, 13563, 'colorthreshold', 0.8);
% 
% % save(['results\hClustResults\donor' num2str(i) '_clusters.mat'], 'clusters');
% % clear corrDist; clear clusterTree; clear clusters;
% 
%     [num txt] = xlsread('temp.xls');
%     colNames = txt(:,1);
%     
%     C = clustergram(donorsExpMat_5RPKM, 'RowLabels', gNames, 'ColumnLabels', colNames, ...
%         'Standardize', 'none', 'Cluster', 1, 'RowPDist', 'spearman', ...
%         'Linkage', 'average', 'Dendrogram', 0.75, 'Colormap', redbluecmap);
% % 
% %     h = plot(C);
% %     saveas(h, ['results\hClustResults\donor' num2str(i) '_dendogram.fig']);
% %     saveas(h, ['results\hClustResults\donor' num2str(i) '_dendogram.jpg']);
% %     save(['results/hClustResults/donor' num2str(i) '_C.mat'], 'C');
% %     close all; clear C;
% 
% %     save([filesDir 'results/hClust/donor' num2str(i) '_clusters.mat'], 'clusters');

%%%-----------------------------------------------

nonCs = [5,7,10,12,16];
Cs = [1,2,3,4,6,8,9,11,13,14,15];

newS = S(nonCs);
% newS(length(newS)+1:16) = S(Cs);

labels(1:11) = S(Cs);
for i = 2 : 30    
    temp(1:11) = S(Cs);
    labels = [labels temp];
    clear temp;
end


% geneMat1 = donorsExpMat_5RPKM(:,nonCs,:);
geneMat1 = donorsExpMat_5RPKM(:,Cs,:);
geneMat1 = reshape(geneMat1, size(geneMat1,1), size(geneMat1,2)*size(geneMat1,3));
clear donorsExpMat_5RPKM;

for i = 1 : size(geneMat1, 2)
    
    T = 1.7;
    clear mat;
    mat(:,:) = geneMat1;
    newS{1} = 'NCX';
    corrDist = pdist(mat, 'spearman');
    clusterTree = linkage(corrDist, 'complete');
    clusters = cluster(clusterTree, 'cutoff', T, 'criterion', 'distance');
    
%     [H,T,perm] = dendrogram(clusterTree, 13563, 'colorthreshold', T);
%     
%     reorderedMat = mat(perm,:);
%     HM = HeatMap(reorderedMat, 'RowLabels', gNames(perm), 'ColumnLabels', labels, ...
%             'Colormap', redbluecmap, 'DisplayRange', 3);

    save(['results\hClustResults\Structures\structuresModules\structures_complete\' newS{i} '_clusters.mat'], 'clusters');
    clear corrDist; clear clusterTree; clear clusters;
    
    [num txt] = xlsread('files\DataFiles\donorSummary.xls');
    donorAge = txt(2:end,4);
%     
%     clear mat;
%     mat(:,:) = geneMat1(:,i,:);
% %  
    C = clustergram(mat, 'RowLabels', gNames, 'ColumnLabels', labels, ...
        'Standardize', 'none', 'Cluster', 1, 'RowPDist', 'spearman', ...
        'Linkage', 'complete', 'Dendrogram', T, 'Colormap', redbluecmap);
% 
    h = plot(C);
%     saveas(h, ['results\hClustResults\Structures\structures_complete\' newS{i} '.fig']);
    saveas(h, ['results\hClustResults\Structures\structures_complete\' newS{i} '.jpg']);
    close all; clear C;
    
%     save([filesDir 'results/hClust/donor' num2str(i) '_clusters.mat'], 'clusters');
%     
%     h = figure;
%     set(h, 'visible', 'off');
%     clustergram(donorsExpMat_5RPKM(:,:,i), 'RowLabels', gNames, 'ColumnLabels', S,...
%         'Colormap',redbluecmap);
%     saveas(h, [filesDir 'results/hClust/donor' num2str(i) '_hClust.fig']);

%     load(['results\hClustResults\hClustTrees\donor' num2str(i) '_clusterTree_AVG.mat']);
%     dendrogram(clusterTree, 13563)
%     set(H,'LineWidth',2)

%     C = clustergram(donorsExpMat_5RPKM(:,:,i), 'RowLabels', gNames, 'ColumnLabels', S, ...
%         'Standardize', 'none', 'Cluster', 1, 'RowPDist', 'spearman', ...
%         'Linkage', 'average');%, 'RowGroupMarker', rm);
    
end

