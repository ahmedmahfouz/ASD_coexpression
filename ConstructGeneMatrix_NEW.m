% %% --Load Expression File--------------------------------------------------
%     filename = 'genes_matrix_csv\genes_matrix.csv';
%     origExpMat = csvread(filename);
%     save('files\origExpMat.mat', 'origExpMat');
% %--------------------------------------------------------------------------

% %% --Load Ages File--------------------------------------------------------
%     ageFile = ('Ages.txt');
%     origAgeList = textread(ageFile, '%s', 'delimiter', '\n');
% 
%     ageCount = 1;
%     temp1 = origAgeList{1};
%     delim1 = strfind(temp1, ' ');
%     temp2= temp1(1:delim1-1);
%     age1 = str2num(temp2);
%     ages(ageCount) = age1;
%     sampleAge(1) = 1;
% 
%     for i = 2 : length(origAgeList)
% 
%         clear temp1; clear temp2; clear delim1; clear age1;
%         temp1 = origAgeList{i};
%         delim1 = strfind(temp1, ' ');
%         temp2= temp1(1:delim1-1);
%         age1 = str2num(temp2);
% 
%         if age1 ~= ages(length(ages))
%             ageCount = ageCount + 1;
%             ages(length(ages)+1) = age1;
%         end
%         sampleAge(i) = ageCount;
%     end
%     clear i; clear j;
%     save('files\sampleAge.mat', 'sampleAge');
%     save('files\ages.mat', 'ages')
% %--------------------------------------------------------------------------
% 
% %% --Load Structures File--------------------------------------------------
%     strucFile = ('Structures.txt');
%     origStrucList = textread(strucFile, '%s', 'delimiter', '\n');
% 
%     strucCount = 1;
%     temp1 = origStrucList{1};
%     strucs{strucCount} = temp1;
%     strucsID(strucCount) = 1;
%     sampleStruc(1) = strucCount;
% 
%     for i = 2 : length(origStrucList)
% 
%         clear temp1; clear temp2; clear delim1; clear age1;
%         temp1 = origStrucList{i};
% 
%         R = 0; j = length(strucs) + 1;
%         while R == 0
%             clear temp2;
%             j = j - 1;
%             temp2 = strucs{j};
%             R = strcmp(temp1, temp2); 
%             if j == 1
%                 break;
%             end        
%         end
% 
%         if R == 1
%             sampleStruc(i) = strucsID(j);
%         else
%             strucCount = strucCount + 1;
%             strucs{strucCount} = temp1;
%             strucsID(strucCount) = strucCount;
%             sampleStruc(i) = strucCount;
%         end
% 
%     end
%     clear i; clear j;
%     save('files\sampleStruc.mat', 'sampleStruc');
%     save('files\strucs.mat', 'strucs');
%     save('files\strucsID.mat', 'strucsID');
% %--------------------------------------------------------------------------
% 
% %% --Load entrezID File----------------------------------------------------
%     idFile = 'genes_matrix_csv\genes_rows_metadata.csv';
%     entrezIDs = csvread(idFile, 1, 4);
%     save('files\entrezIDs.mat', 'entrezIDs');
% %--------------------------------------------------------------------------

%%% --Load Donor IDs-------------------------------------------------------
%     donorIDs = csvread('DonorID.csv');
%     save('files\donorIDs.mat', 'donorIDs');
%--------------------------------------------------------------------------

%%% --Load Gene Names------------------------------------------------------
%     [num gNames] = xlsread('files\DataFiles\gNames.xls');
%     save('files\DataFiles\gNames.mat', 'gNames');
%--------------------------------------------------------------------------

%%% --Load ensemble_IDs----------------------------------------------------
%     [num enIDs] = xlsread('files\DataFiles\ensembl_IDs.xls');
%     save('files\enIDs.mat', 'enIDs');
%--------------------------------------------------------------------------

%%% --Construct Donors' Expreseeion Matrix (rows:Genes; columns:Structures; depth:Donors)

clear all;

load('files\origExpMat.mat');
load('files\sampleAge');
load('files\sampleStruc');
load('files\ages.mat');
load('files\strucs.mat');
load('files\DonorIDs.mat')
[num txt] = xlsread('files\DataFiles\donorsNames&IDs.xls');
txt(:,1) = [];
num(:,2:3) = [];

dAge = num(:,2);
[dAge_S sorter] = sort(dAge);
dIDs = num(:,1);
dIDs_S = dIDs(sorter);
dNames = txt(:,1);
dNames_S = dNames(sorter);
dAges = txt(:,2);
dAges_S = dAges(sorter);

% xlswrite('files\DataFiles\donorsNames&IDs.xls', dIDs_S, 1, 'F1');
% xlswrite('files\DataFiles\donorsNames&IDs.xls', dNames_S, 1, 'G1');
% xlswrite('files\DataFiles\donorsNames&IDs.xls', dAges_S, 1, 'H1');


strucIndInc = [2, 3, 5, 6, 7, 8, 9, 15, 18, 19, 20, 21, 22, 23, 24, 26];

expMat = zeros(size(origExpMat, 1), length(strucIndInc), length(dIDs_S));

for D = 1 : length(dIDs_S)

    clear currD; clear cuurD_samples; clear currD_strucs;
    
    currD = dIDs_S(D);
    cuurD_samples = find(donorIDs == currD);
    currD_strucs = sampleStruc(cuurD_samples);

    for j = 1 : length(strucIndInc)

        clear currStruc; clear test; clear currD_ind; clear tempExp1;

        currStruc = strucIndInc(j);
        test = find(currD_strucs == currStruc);

        if ~isempty(test);
            
            currD_ind = cuurD_samples(test);
            tempExp1 = origExpMat(:, currD_ind);
            
            expMat(:, j, D) = tempExp1;
            
        end

    end

end

save('files\donorsExpMat.mat', 'expMat');

%%% correct for missing samples
load('files\donorsExpMat.mat');
%%% (1) 16pcw [Donor: 9, Structure: 3,5,8,9,15,16]
mS = [8,9,15,16];
iD = 9; sD = 10;
expMat(:, mS, iD) = imv_v1(expMat, iD, sD, mS);
clear retMat; clear mS; clear iD; clear sD;

mS = [3, 5];
iD = 9; sD = [10,11];
for i = 1 : length(sD)
    retMat(:,:,i) = imv_v1(expMat, iD, sD(i), mS);
end
expMat(:, mS, iD) = mean(retMat,3);
clear retMat; clear mS; clear iD; clear sD;

%%% (2) 16pcw [Donor: 11, Structure: 8,9,15,16]
mS = [8,9,15,16];
iD = 11; sD = 10;
expMat(:, mS, iD) = imv_v1(expMat, iD, sD, mS);
clear retMat; clear mS; clear iD; clear sD;

%%% (3) 21pcw [Donor: 15, Structure: 13]
mS = [13];
iD = 15; sD = 14;
expMat(:, mS, iD) = imv_v1(expMat, iD, sD, mS);
clear retMat; clear mS; clear iD; clear sD;

%%% (4) 3yrs [Donor: 25, Structure: 2,6,8]
mS = [2,6,8];
iD = 25; sD = 26;
expMat(:, mS, iD) = imv_v1(expMat, iD, sD, mS);
clear retMat; clear mS; clear iD; clear sD;

%%% (5) 3yrs [Donor: 26, Structure: 10,12]
mS = [10,12];
iD = 26; sD = 25;
expMat(:, mS, iD) = imv_v1(expMat, iD, sD, mS);
clear retMat; clear mS; clear iD; clear sD;

%%% (6) 8yrs [Donor: 29, Structure: 9,10,12,15]
mS = [9,10,12,15];
iD = 29; sD = 28;
expMat(:, mS, iD) = imv_v1(expMat, iD, sD, mS);
clear retMat; clear mS; clear iD; clear sD;

%%% (7) 40yrs [Donor: 41, Structure: 6]
mS = [6];
iD = 41; sD = 40;
expMat(:, mS, iD) = imv_v1(expMat, iD, sD, mS);
clear retMat; clear mS; clear iD; clear sD;

%%% (8) 17pcw [Donor: 12; Structure: 3,9,15; w1: 9; w2: 16]
mS = [3,9,15];
iD = 12; sD = [9:16];
for i = 1 : length(sD)
    retMat(:,:,i) = imv_v1(expMat, iD, sD(i), mS);
end
expMat(:, mS, iD) = mean(retMat,3);
clear retMat; clear mS; clear iD; clear sD;

%%% (9) 19pcw [Donor: 13; Structure: 3,5,8,9,15,16; w1: 9; w2: 16]
mS = [3,5,8,9,15,16];
iD = 13; sD = [9:16];
for i = 1 : length(sD)
    retMat(:,:,i) = imv_v1(expMat, iD, sD(i), mS);
end
expMat(:, mS, iD) = mean(retMat,3);
clear retMat; clear mS; clear iD; clear sD;

%%% (10) 11yrs [Donor: 30; Structure: 10,12; w1: 28; w2: 35]
mS = [10,12];
iD = 30; sD = [28:35];
for i = 1 : length(sD)
    retMat(:,:,i) = imv_v1(expMat, iD, sD(i), mS);
end
expMat(:, mS, iD) = mean(retMat,3);
clear retMat; clear mS; clear iD; clear sD;

%%% (11) 15yrs [Donor: 32; Structure: 7; w1: 28; w2: 35]
mS = [7];
iD = 32; sD = [28:35];
for i = 1 : length(sD)
    retMat(:,:,i) = imv_v1(expMat, iD, sD(i), mS);
end
expMat(:, mS, iD) = mean(retMat,3);
clear retMat; clear mS; clear iD; clear sD;

%%% (12) 18yrs [Donor: 33; Structure: 5,10,12; w1: 28; w2: 35]
mS = [5,10,12];
iD = 33; sD = [28:35];
for i = 1 : length(sD)
    retMat(:,:,i) = imv_v1(expMat, iD, sD(i), mS);
end
expMat(:, mS, iD) = mean(retMat,3);
clear retMat; clear mS; clear iD; clear sD;

neg = find(expMat <= 0);
expMat(neg) = 0;

save('files\donorsExpMat.mat', 'expMat');
    
%%% Remove inconsistent Donors
load('files\donorsExpMat.mat');
donorsToRemove = [1, 2, 3, 4, 5, 6, 7, 8, 17, 18, 22];
expMat(:,:,donorsToRemove) = [];

save('files\donorsExpMat.mat', 'expMat');
%--------------------------------------------------------------------------

%%% Remove "non-expressing" Genes
clear all;
% load('files\origExpMat.mat');
load('files\donorsExpMat.mat');

for g = 1 : size(expMat, 1)
    
%     clear temp1; 
%     temp1(:,:) = expMat(g,:,:);
%     temp1 = reshape(temp1, 1, size(temp1,1)*size(temp1,2));
%     CV(g) = std(temp1) / mean(temp1);
    
    R = find(expMat(g,:,:) >= 5);
    gE(g) = length(R);
        
end

lowExp = find(gE <= 0);
highExp = find(gE > 0);

% lowCV = find(CV < 0.5);
% highCV = find(CV > 0.5);

nonExpGenes = lowExp;
expGenes = highExp;

genesStatus_5RPKM = zeros(length(gE),1);
genesStatus_5RPKM(expGenes) = 1;

donorsExpMat_5RPKM = expMat(expGenes, :, :);

save('files\donorsExpMat_5RPKM.mat', 'donorsExpMat_5RPKM');
save('files\genesStatus_5RPKM.mat', 'genesStatus_5RPKM');
%%%------------------------------------------------------------------------


%%% Save expression data into xls sheets
clear all;

load('files\donorsExpMat_5RPKM.mat');
load('files\gNames.mat');
load('files\genesStatus_5RPKM.mat');
load('files\strucs.mat');

% list of included structures + the entrezID header
strucIndInc = [2, 3, 5, 6, 7, 8, 9, 15, 18, 19, 20, 21, 22, 23, 24, 26];
for i = 1 : length(strucIndInc)
    S{i+1} = strucs{strucIndInc(i)};
end
S{1} = 'Gene Name';

gIND = find(genesStatus_5RPKM == 1);
gNames_5RPKM = gNames(gIND);

% fname = ['Data\Donors\All_log2.xls'];
fname = ['Data\Donors\sepFiles\'];

% Z = find(donorsExpMat_5RPKM == 0);
% donorsExpMat_5RPKM(Z) = 8.2924e-005;

for i = 1 : size(donorsExpMat_5RPKM, 3)
    
    clear dfname;
    dfname = ['Data\Donors\sepFiles\D' num2str(i) '.xls'];
    xlswrite(dfname, log2(donorsExpMat_5RPKM(:,:,i)+(1*(10^-5))), i, 'B2');
    xlswrite(dfname, gNames_5RPKM, i, 'A2');
    xlswrite(dfname, S, i, 'A1');
    
end
%%------------------------------------------------------------------------



