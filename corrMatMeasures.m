% Measurements on the correlation matrices

clear all;

% load('files\newGeneMat_5RPKM.mat');
load('files\asdGeneIDs_5RPKM.mat');
load('files\szGeneIDs_5RPKM.mat');
load('files\ndGeneIDs_5RPKM.mat');
load('files\hkGeneIDs_5RPKM.mat');


% for i = 1 : 4
%       
%     % ASDs
%     load(['results\diffNetworks\abs_abs\asdDiffNet_' num2str(i+1) num2str(i) '.mat']);
%     temp1 = triu(asdDiffMat,1); 
%     nonZ = find(temp1 ~= 0);
%     temp2 = temp1(nonZ);
%     N = isnan(temp2);
%     NN = find(N == 1);
%     temp2(NN) = [];
%     asdM(i) = mean(temp2);
%     clear asdDiffMat; clear temp1; clear temp2; clear nonZ;
%     
%     % SZ
%     load(['results\diffNetworks\abs_abs\szDiffNet_' num2str(i+1) num2str(i) '.mat']);
%     temp1 = triu(szDiffMat,1); 
%     nonZ = find(temp1 ~= 0);
%     temp2 = temp1(nonZ);
%     N = isnan(temp2);
%     NN = find(N == 1);
%     temp2(NN) = [];
%     szM(i) = mean(temp2);
%     clear szDiffMat; clear temp1; clear temp2; clear nonZ;
%     
%     % ND
%     load(['results\diffNetworks\abs_abs\ndDiffNet_' num2str(i+1) num2str(i) '.mat']);
%     temp1 = triu(ndDiffMat,1); 
%     nonZ = find(temp1 ~= 0);
%     temp2 = temp1(nonZ);
%     N = isnan(temp2);
%     NN = find(N == 1);
%     temp2(NN) = [];
%     ndM(i) = mean(temp2);
%     clear ndDiffMat; clear temp1; clear temp2; clear nonZ;
%     
% end
%%%------------------------------------------------------------------------

%%%------------------------------------------------------------------------
%%% Hirarchial Clustering
%%%------------------------------------------------------------------------
% load(['results\geneCorr\Avg(log2)\ageGroup1_corrMat.mat']);
% 
% Z = linkage(corrMat);
% save('results\geneCorr\Avg(log2)\ageGroup1_Z.mat', 'Z');

load('results\geneCorr\Avg(log2)\ageGroup1_Z.mat');
% set(0,'RecursionLimit',1000);

dendrogram(zNEW)








