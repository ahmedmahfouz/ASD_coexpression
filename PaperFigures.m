%%% PaperFigures

clear all;

%%% Gene Names
load(['files/gNames.mat']);
load(['files/genesStatus_5RPKM.mat']);
gNames(find(genesStatus_5RPKM == 0)) = [];

%%% Structures Names
load(['files/strucs.mat']);
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
%%%------------------------------------------------------------------------

% load('files/donorsExpMat_5RPKM.mat');
% newExpMat = donorsExpMat_5RPKM(:,nonCs,:);
% tempMat = donorsExpMat_5RPKM(:,Cs,:);
% tempMat = mean(tempMat, 2);
% newExpMat(:,end+1,:) = tempMat;
% clear tempMat;
% newExpMat = normalizeExpMat(newExpMat);
% newExpMat = log2(newExpMat + (rand(size(newExpMat))*(10^-5)));
% clear donorsExpMat_5RPKM;
% 
% for i = 1 : length(group)
%     
%     tempMat = newExpMat(:,:,group{i});
%     avgMat = mean(tempMat, 3);
%     C = clustergram(avgMat, 'Standardize', 'none', ...
%         'Cluster', 2, 'Colormap', redbluecmap);
%     h = plot(C);
%     saveas(h, ['figures\ageStage_' num2str(i) '_HM.fig']);
%     saveas(h, ['results\ageStage_' num2str(i) '_HM.jpg']);
%     close all; clear C;
%     
% end

%%%------------------------------------------------------------------------

D = 'ASDs';
%%% load the gene list correlation matrix
load(['results\corrChange\' D 'Corr.mat']);
t = ones(size(glCorr,1), size(glCorr,2)); 
t = triu(t); t = t * 10000;
p = find(t ~= 10000); clear t;
for i = 1 : 7    
    clear corrDist;
    corrDist(:,:) = glCorr(:,:,i);
    modCorr(:,i) = corrDist(p);
end

temp = modCorr;
T = 0.8;
for k = 1 : size(temp,1)
    x1 = find(temp(k,:) >= T);
    spC(k) = length(x1);
    x2 = find(temp(k,:) <= (-1*T));
    snC(k) = length(x2);
end
spcID = find(spC > 0);
sncID = find(snC > 0);
sConn = temp([spcID sncID], :);
newP = p([spcID sncID]);

load(['results\corrChange\clusters_' D '.mat']);
clustInd3 = find(clusters == 3);
cC3 = mean(sConn(clustInd3,:), 1);
clustInd4 = find(clusters == 4);
cC4 = mean(sConn(clustInd4,:), 1);

f = figure;
hold on;
% plot(sConn', 'linewidth', 0.25, 'Color', [0.6 1 0.6])
plot(sConn(clustInd3,:)', 'linewidth', 0.25, 'Color', [0.8 0.8 0.8]);
plot(sConn(clustInd4,:)', 'linewidth', 0.25, 'Color', [0.8 0.8 0.8]);
plot(cC3, 'linewidth', 4, 'Color', 'red'); 
plot(cC4, 'linewidth', 4, 'Color', 'blue'); 
axis([0.8 7.2 0 1])
grid on
hold off



