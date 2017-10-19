%%% plot expression distribution of donors

clear all;

load('files\donorsExpMat_5RPKM.mat');

expMat = normalizeExpMat(donorsExpMat_5RPKM);

for i = 1 : size(donorsExpMat_5RPKM, 3)

%     clear geneMat1; clear geneMat2;
%     geneMat1 = log2(donorsExpMat_5RPKM(:,:,i) + (1*10^-5));
    geneMat2 = log2(expMat(:,:,i) + 5);
%     RHO = corr(geneMat, 'Type', 'Spearman');
    
    temp1 = reshape(geneMat2, 1, size(geneMat2,1)*size(geneMat2,2));
%     temp2 = log2(temp1);
%     temp3 = log2(temp1 + (1*10^-5));
    
%     temp2 = reshape(corr(log2(geneMat), 'Type', 'Spearman'), 1, 16*16);
%     temp3 = reshape(corr(geneMat2, 'Type', 'Spearman'), 1, 16*16);
    
%     temp1 = reshape(geneMat1, 1, size(geneMat1,1)*size(geneMat1,2));
%     temp2 = reshape(geneMat2, 1, size(geneMat2,1)*size(geneMat2,2));
    
    xA = 0 : 0.1 : 20;
%     xA = 0.5 : 0.05 : 1;
    n1 = histc(temp1, xA);
%     n2 = histc(temp2, xA);
%     n3 = histc(temp3, xA);
    
    figure(1), 
    subplot(6,5,i)
    hold on
%     bar(xA, n1, 'b')
%     bar(xA, n2, 0.4, 'r')
    
    plot(xA, n1, 'linewidth', 2, 'color', 'blue');
%     plot(xA, n3, 'linewidth', 3, 'linestyle', '--', 'color', 'yellow');
    grid on, 
%     axis([-20 20 0 4*10^4]),
%     axis([0 1 0 150]),
    title(['Donor' num2str(i)])
%     xlabel('Log_2(RPKM)'); ylabel('Density');
%     legend('Log_2(RPKM)', 'Log_2(RPKM+(1*10-5))')
    hold off
    
end


% load('files\donorsExpMat_5RPKM.mat');
% expMat = normalizeExpMat(donorsExpMat_5RPKM);
% 
% nonCs = [5,7,10,12,16];
% Cs = [1,2,3,4,6,8,9,11,13,14,15];
% 
% geneMat1 = donorsExpMat_5RPKM(:,nonCs,:);
% geneMat1(:,size(geneMat1,2)+1,:) = mean(donorsExpMat_5RPKM(:,Cs,:), 2);
% 
% geneMat2 = expMat(:,nonCs,:);
% geneMat2(:,size(geneMat2,2)+1,:) = mean(expMat(:,Cs,:), 2);
% 
% geneMat1 = log2(geneMat1 + (1*10^-5));
% geneMat2 = log2(geneMat2 + (1*10^-5));
% 
% S = {'AMY', 'HIP', 'STR', 'MD', 'CBC', 'NCx'};
% colors = hsv(6);
% f= figure;
% 
% for i = 1 : 30
%     
%     tempMat = geneMat1(:,:,i);
%     boxplot(tempMat);
%     
%     subplot(2,1,1), hold on
%     for j = 1 : 6
%     
%         plot(geneMat1(:,j,i), 'Color', colors(j,:), 'LineWidth', 2), grid on
%         
%     end
%     xlabel('Genes', 'fontweight', 'bold'); 
%     ylabel('log_2(Expression)', 'fontweight', 'bold');
%     title(['Donor' num2str(i)], 'fontweight', 'bold');
%     legend(S)
%     hold off
%     
% end

% 
% load('files\donorsExpMat.mat');
% load('files\genesStatus_5RPKM.mat');
% 
% expMat = reshape(expMat, size(expMat,1), size(expMat,2)*size(expMat,3));
% expMat = log2(expMat + (1*10^-5));
% gM = mean(expMat');
% 
% figure, 
% boxplot(gM, genesStatus_5RPKM')







