% Data Variance and threshold analysis

clear all;

%%%------------------------------------------------------------------------
%%% Overall log10(Expression) Distribution
%%%------------------------------------------------------------------------
% load('files\geneMat.mat')
% tempGeneMat = reshape(geneMat, 1, size(geneMat,1)*size(geneMat,2)*size(geneMat,3));
% clear geneMat;
% 
% temp2 = log10(tempGeneMat);
% 
% xA = -0.05 : 0.05 : max(temp2);
% n = histc(temp2, xA);
% figure, bar(xA, n), grid on;

% for i = 1 : length(n)
%     
%     nG = sum(n(1:i));
%     pG(i) = nG/sum(n);
%         
% end
%%%------------------------------------------------------------------------

%%%------------------------------------------------------------------------
%%% Overall log2(Expression) Distribution
%%%------------------------------------------------------------------------
load('files\donorsExpMat.mat')
load('files\genesStatus_5RPKM');
load('files\hkGeneIDs.mat');
[hkSL_IDs hkSL_gNames] = xlsread('files\DataFiles\HKshortlist.xls');

eG = find(genesStatus_5RPKM == 1);
egMat = expMat(eG,:,:);
tempGeneMat1 = reshape(egMat, 1, size(egMat,1)*size(egMat,2)*size(egMat,3));

neG = find(genesStatus_5RPKM == 0);
negMat = expMat(neG,:,:);
tempGeneMat2 = reshape(negMat, 1, size(negMat,1)*size(negMat,2)*size(negMat,3));

% hkgMat = expMat(hkSL_IDs,:,:);
hkgMat = expMat(hkIDs,:,:);
tempGeneMat3 = reshape(hkgMat, 1, size(hkgMat,1)*size(hkgMat,2)*size(hkgMat,3));

agMat = expMat;
tempGeneMat4 = reshape(agMat, 1, size(agMat,1)*size(agMat,2)*size(agMat,3));

clear egMat; clear negMat;
clear donorsExpMat; clear genesStatus_5RPKM;

temp1 = log2(tempGeneMat1); 
temp11 = log2(tempGeneMat1+(1*10^-5));
temp2 = log2(tempGeneMat2);
temp3 = log2(tempGeneMat3);
temp4 = log2(tempGeneMat4);

xA = -10 : 0.5 : 20;
n1 = histc(temp1, xA);
n11 = histc(temp11, xA);
n2 = histc(temp2, xA);
n3 = histc(temp3, xA);
n4 = histc(temp4, xA);
figure, hold on
plot(xA, n1, 'linewidth', 3, 'color', 'blue');
plot(xA, n11, 'linewidth', 3, 'linestyle', '--', 'color', 'yellow');
% plot(xA, n2, 'linewidth', 3, 'color', 'red');
% plot(xA, n3, 'linewidth', 3, 'color', 'green');
% plot(xA, n4, 'linewidth', 3, 'color', 'blue');
grid on, title('Distribution plot of RPKM')
xlabel('Log_2(RPKM)'); ylabel('Density');
% legend('expressed', 'non-expressed', 'housekeeping genes')
% legend('housekeeping genes (short list)')
hold off
%%%------------------------------------------------------------------------

%%%------------------------------------------------------------------------
%%% Overall log10(Expression) Distribution and Variance of Housekeeping
%%% Genes
%%%------------------------------------------------------------------------
% load('files\geneMat.mat')
% load('files\hkGeneIDs.mat')

%%%Overall Distribution
% tempGeneMat = geneMat(:,:,hkIDs);
% tempGeneMat = reshape(tempGeneMat, 1, size(tempGeneMat,1)*size(tempGeneMat,2)*size(tempGeneMat,3));
% 
% temp1 = log10(tempGeneMat);
% 
% xA = -0.05 : 0.05 : max(temp1);
% n = histc(temp1, xA);
% figure, bar(xA, n), grid on;
% title('Distribution of log10(Expression) of 1830 Housekeeping Genes');
% xlabel('log10(RPKM)'); ylabel('Density');
% 
% %%% max Expression per gene
% tempGeneMat = geneMat(:,:,hkIDs);
% temp2 = reshape(tempGeneMat, size(tempGeneMat,1)*size(tempGeneMat,2), size(tempGeneMat,3));
% 
% maxExp = log10(max(temp2, [], 1));
% xA = -0.05 : 0.25 : max(maxExp);
% n = histc(maxExp, xA);
% figure, bar(xA, n), grid on;
% title('Distribution of log10(Max. Expression) of 1830 Housekeeping Genes');
% xlabel('log10(RPKM)'); ylabel('# of Genes');

% %%% max Expression per gene
% tempGeneMat = geneMat(:,:,hkIDs);
% tempGeneMat = reshape(tempGeneMat, size(tempGeneMat,1)*size(tempGeneMat,2), size(tempGeneMat,3));
% 
% varExp = var(tempGeneMat, [], 1);
% xA = 0 : 1 : max(varExp);
% n = histc(varExp, xA);
% % figure, bar(xA, n), grid on;
%%%------------------------------------------------------------------------

%%%------------------------------------------------------------------------
%%% Number of genes removed by applying a threshold: 1RPKM & 5RPKM
%%%------------------------------------------------------------------------
% load('files\geneMat.mat')
% tempGeneMat = reshape(geneMat, size(geneMat,1)*size(geneMat,2), size(geneMat,3));
% maxExp = min(tempGeneMat, [], 1);
% 
% xA = 0 : 1 : max(maxExp);
% n = histc(maxExp, xA);
% figure, bar(xA, n), grid on;
% 
% RPKM1 = find(xA == 1);
% N_1RPKM = sum(n(1:RPKM1));
% 
% RPKM5 = find(xA == 5);
% N_5RPKM = sum(n(1:RPKM5));
%%%------------------------------------------------------------------------

%%%------------------------------------------------------------------------
%%% Overall variance per gene
%%%------------------------------------------------------------------------
% load('files\geneMat.mat')
% tempGeneMat = reshape(geneMat, size(geneMat,1)*size(geneMat,2), size(geneMat,3));
% varExp = var(log2(tempGeneMat));
% maxExp = max(tempGeneMat, [], 1);
% 
% 
% temp = find(isnan(varExp) == 1);
% varExp(temp) = [];
% maxExp(temp) = [];
% 
% figure, grid on, plot(varExp, 'r');
% figure, grid on, plot(maxExp, 'b');
% 
% xA = 0 : 0.1 : max(varExp);
% n = histc(varExp, xA);
% figure, bar(xA, n), grid on;
% title('Variance Histogram of log2(Expression)', 'fontweight', 'bold');
% xlabel('Variance', 'fontweight', 'bold');
% ylabel('Density', 'fontweight', 'bold')
%%%------------------------------------------------------------------------

%%%------------------------------------------------------------------------
%%% Spatial variance per gene
%%%------------------------------------------------------------------------
% tempGeneMat = mean(geneMat, 2);
% tempGeneMat = reshape(tempGeneMat, size(tempGeneMat,1), size(tempGeneMat,3));
% varExp = var(tempGeneMat);
% 
% xA = 0 : 0.1 : 2;
% n = histc(varExp, xA);
% figure, bar(xA, n), grid on;
%%%------------------------------------------------------------------------

%%%------------------------------------------------------------------------
%%% Temporal variance per gene
%%%------------------------------------------------------------------------
% tempGeneMat = mean(geneMat, 1);
% tempGeneMat = reshape(tempGeneMat, size(tempGeneMat,2), size(tempGeneMat,3));
% varExp = var(tempGeneMat);
% 
% xA = 0 : 0.1 : 2;
% n = histc(varExp, xA);
% figure, bar(xA, n), grid on;
%%%------------------------------------------------------------------------

%%%------------------------------------------------------------------------
%%% 2D Gene expression change over time (different plot per group)
%%%------------------------------------------------------------------------
% load('files\newGeneMat_5RPKM.mat');
% load('files\asdGeneIDs_5RPKM.mat');
% load('files\szGeneIDs_5RPKM.mat');
% load('files\ndGeneIDs_5RPKM.mat');
% load('files\hkGeneIDs_5RPKM.mat');
% 
% colors = jet(5);
% 
% group = {'ASD', 'SZ', 'ND', 'HK'};
% 
% for i = 1 : length(newGeneMat)
%     
%     geneMat = newGeneMat{i};
%     tempGeneMat = mean(geneMat, 3);
%     
%     figure(1), hold on
%     
%     % ASDs
%     asdGenes = tempGeneMat(asdIDs_5RPKM,:);
%     asdMean = mean(asdGenes, 1);
%     
%     subplot(2,2,1), hold on, grid on;
%     plot(asdMean, 'linewidth', 2, 'color', colors(i,:))
%     legend('2nd Trimester', '1st Year', 'Childhood', 'Adolescence', 'Adulthood');
%     axis([0 17 0 50])
%     title(group{1}); xlabel('Time Point'); ylabel('Exoression (RPKM)');
%     hold off
%     
%     % SZ
%     szGenes = tempGeneMat(szIDs_5RPKM,:);
%     szMean = mean(szGenes, 1);
%     
%     subplot(2,2,2), hold on, grid on;
%     plot(szMean, 'linewidth', 2, 'color', colors(i,:))
%     legend('2nd Trimester', '1st Year', 'Childhood', 'Adolescence', 'Adulthood');
%     axis([0 17 0 50])
%     title(group{2}); xlabel('Time Point'); ylabel('Exoression (RPKM)');
%     hold off
%     
%     % ND
%     ndGenes = tempGeneMat(ndIDs_5RPKM,:);
%     ndMean = mean(ndGenes, 1);
%     
%     subplot(2,2,3), hold on, grid on;
%     plot(ndMean, 'linewidth', 2, 'color', colors(i,:))
%     legend('2nd Trimester', '1st Year', 'Childhood', 'Adolescence', 'Adulthood');
%     axis([0 17 0 50])
%     title(group{3}); xlabel('Time Point'); ylabel('Exoression (RPKM)');
%     hold off
%     
%     % HK
%     hkGenes = tempGeneMat(hkIDs_5RPKM,:);
%     hkMean = mean(hkGenes, 1);
%     
%     subplot(2,2,4), hold on, grid on;
%     plot(hkMean, 'linewidth', 2, 'color', colors(i,:))
%     legend('2nd Trimester', '1st Year', 'Childhood', 'Adolescence', 'Adulthood');
%     axis([0 17 0 50])
%     title(group{4}); xlabel('Time Point'); ylabel('Exoression (RPKM)');
%     hold off
%     
% end
%%%------------------------------------------------------------------------

%%%------------------------------------------------------------------------
%%% 2D Gene expression change over time (different plot per group)
%%%------------------------------------------------------------------------
% load('files\newGeneMat_5RPKM.mat');
% load('files\asdGeneIDs_5RPKM.mat');
% load('files\szGeneIDs_5RPKM.mat');
% load('files\ndGeneIDs_5RPKM.mat');
% load('files\hkGeneIDs_5RPKM.mat');
% 
% colors = jet(4);
% age = {'2nd Trimester', '1st Year', 'Childhood', 'Adolescence', 'Adulthood'};
% 
% for i = 1 : length(newGeneMat)
%     
%     geneMat = newGeneMat{i};
%     tempGeneMat = mean(geneMat, 3);
%     
%     % ASDs
%     asdGenes = tempGeneMat(asdIDs_5RPKM,:);
%     asdMean = mean(asdGenes, 1);
%     
%     % SZ
%     szGenes = tempGeneMat(szIDs_5RPKM,:);
%     szMean = mean(szGenes, 1);
%     
%     % ND
%     ndGenes = tempGeneMat(ndIDs_5RPKM,:);
%     ndMean = mean(ndGenes, 1);
%     
%     % HK
%     hkGenes = tempGeneMat(hkIDs_5RPKM,:);
%     hkMean = mean(hkGenes, 1);
%     
%     figure(i), grid on, 
%     hold on
%     plot(asdMean, 'linewidth', 2, 'color', colors(1,:))
%     plot(szMean, 'linewidth', 2, 'color', colors(2,:))
%     plot(ndMean, 'linewidth', 2, 'color', colors(3,:))
%     plot(hkMean, 'linewidth', 2, 'color', colors(4,:))
%     hold off
%     legend('ASD', 'SZ', 'ND', 'HK');
%     axis([0 17 0 80])
%     title(age{i}); xlabel('Time Point'); ylabel('Exoression (RPKM)');
%     
% end
%%%------------------------------------------------------------------------

%%%------------------------------------------------------------------------
%%% 3D Gene expression change over time (different plot per group)
%%%------------------------------------------------------------------------
% load('files\newGeneMat_5RPKM.mat');
% load('files\asdGeneIDs_5RPKM.mat');
% load('files\szGeneIDs_5RPKM.mat');
% load('files\ndGeneIDs_5RPKM.mat');
% load('files\hkGeneIDs_5RPKM.mat');
% 
% for i = 1 : length(newGeneMat)
%     
%     geneMat = newGeneMat{i};
%     tempGeneMat = mean(geneMat, 3);
%     
%     newGeneMatAvg(:,:, i) = tempGeneMat;
%     
% end
% clear newGeneMat; clear geneMat; clear tempGeneMat;
% 
% % ASDs
% asdGenes = newGeneMatAvg(asdIDs_5RPKM, :, :);
% asdMean(:,:) = mean(asdGenes, 1);
% 
% f1 = figure(1), grid on, 
% hold on
% bar3(asdMean, 0.25, 'detached')
% legend('2nd Trimester', '1st Year', 'Childhood', 'Adolescence', 'Adulthood');
% axis([0 6 0 17 0 80])
% xlabel('Ages'); ylabel('Structures'); zlabel('Avg Expression');
% title('ASDs'), view(3), hold off
% saveas(f1, 'figures\3D Expression Profile (ASDs).jpg');
% saveas(f1, 'figures\3D Expression Profile (ASDs).fig');
% 
% % SZ
% szGenes = newGeneMatAvg(szIDs_5RPKM, :, :);
% szMean(:,:) = mean(szGenes, 1);
% 
% f2 = figure(2), grid on, 
% hold on
% bar3(szMean, 0.25, 'detached')
% legend('2nd Trimester', '1st Year', 'Childhood', 'Adolescence', 'Adulthood');
% axis([0 6 0 17 0 80])
% xlabel('Ages'); ylabel('Structures'); zlabel('Avg Expression');
% title('SZ'), view(3), hold off;
% saveas(f2, 'figures\3D Expression Profile (SZ).jpg');
% saveas(f2, 'figures\3D Expression Profile (SZ).fig');
% 
% % ND
% ndGenes = newGeneMatAvg(ndIDs_5RPKM, :, :);
% ndMean(:,:) = mean(ndGenes, 1);
% 
% f3 = figure(3), grid on, 
% hold on
% bar3(ndMean, 0.25, 'detached')
% legend('2nd Trimester', '1st Year', 'Childhood', 'Adolescence', 'Adulthood');
% axis([0 6 0 17 0 80])
% xlabel('Ages'); ylabel('Structures'); zlabel('Avg Expression');
% title('ND'), view(3), hold off;
% saveas(f3, 'figures\3D Expression Profile (ND).jpg');
% saveas(f3, 'figures\3D Expression Profile (ND).fig');
% 
% % HK
% hkGenes = newGeneMatAvg(hkIDs_5RPKM, :, :);
% hkMean(:,:) = mean(hkGenes, 1);
% 
% f4 = figure(4), grid on, 
% hold on
% bar3(hkMean, 0.25, 'detached')
% legend('2nd Trimester', '1st Year', 'Childhood', 'Adolescence', 'Adulthood');
% axis([0 6 0 17 0 80])
% xlabel('Ages'); ylabel('Structures'); zlabel('Avg Expression');
% title('HK'), view(3), hold off;
% saveas(f4, 'figures\3D Expression Profile (HK).jpg');
% saveas(f4, 'figures\3D Expression Profile (HK).fig');
%%%------------------------------------------------------------------------

%%%------------------------------------------------------------------------
%%% 3D Gene expression change over time (different plot per age)
%%%------------------------------------------------------------------------
% load('files\newGeneMat_5RPKM.mat');
% load('files\asdGeneIDs_5RPKM.mat');
% load('files\szGeneIDs_5RPKM.mat');
% load('files\ndGeneIDs_5RPKM.mat');
% load('files\hkGeneIDs_5RPKM.mat');
% 
% name = {'2nd Trimester', '1st Year', 'Childhood', 'Adolescence', 'Adulthood'};
% 
% for i = 1 : length(newGeneMat)
%     
%     geneMat = newGeneMat{i};
%     tempGeneMat = mean(geneMat, 3);
%     
%     % ASDs
%     asdGenes = tempGeneMat(asdIDs_5RPKM, :);
%     asdMean = mean(asdGenes, 1);
%     
%     % SZ
%     szGenes = tempGeneMat(szIDs_5RPKM, :);
%     szMean = mean(szGenes, 1);
%     
%     % ND
%     ndGenes = tempGeneMat(ndIDs_5RPKM, :);
%     ndMean = mean(ndGenes, 1);
%     
%     % HK
%     hkGenes = tempGeneMat(hkIDs_5RPKM, :);
%     hkMean = mean(hkGenes, 1);
%     
%     ageGroup = [asdMean; szMean; ndMean; hkMean];
%     
%     f(i) = figure(i), grid on, 
%     hold on
%     bar3(ageGroup', 0.125, 'detached')
%     hold off
%     legend('ASDs', 'SZ', 'ND', 'HK');
%     axis([0 4 0 17 0 80])
%     xlabel('Group'); ylabel('Structures'); zlabel('Avg Expression');
%     title(name{i});
%     saveas(f(i), ['figures\3D Expression Profile (' name{i} ').jpg']);
%     saveas(f(i), ['figures\3D Expression Profile (' name{i} ').fig']);
%     
% end
% clear newGeneMat; clear geneMat; clear tempGeneMat;
%%%------------------------------------------------------------------------



