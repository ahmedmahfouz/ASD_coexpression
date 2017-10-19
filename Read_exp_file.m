% Read gene_matrix.csv

clear all;

%%%--Load Expression File--------------------------------------------------
% filename = 'genes_matrix_csv\genes_matrix.csv';
% origExpMat = csvread(filename);
% save('files\origExpMat.mat', 'origExpMat');
%--------------------------------------------------------------------------

%%%--Load Ages File--------------------------------------------------------
% ageFile = ('Ages.txt');
% origAgeList = textread(ageFile, '%s', 'delimiter', '\n');
% 
% ageCount = 1;
% temp1 = origAgeList{1};
% delim1 = strfind(temp1, ' ');
% temp2= temp1(1:delim1-1);
% age1 = str2num(temp2);
% ages(ageCount) = age1;
% sampleAge(1) = 1;
% 
% for i = 2 : length(origAgeList)
%     
%     clear temp1; clear temp2; clear delim1; clear age1;
%     temp1 = origAgeList{i};
%     delim1 = strfind(temp1, ' ');
%     temp2= temp1(1:delim1-1);
%     age1 = str2num(temp2);
%         
%     if age1 ~= ages(length(ages))
%         ageCount = ageCount + 1;
%         ages(length(ages)+1) = age1;
%     end
%     sampleAge(i) = ageCount;
% end
% clear i; clear j;
% save('files\sampleAge.mat', 'sampleAge');
% save('files\ages.mat', 'ages')
%--------------------------------------------------------------------------

%%%--Load Structures File--------------------------------------------------------
% strucFile = ('Structures.txt');
% origStrucList = textread(strucFile, '%s', 'delimiter', '\n');
% 
% strucCount = 1;
% temp1 = origStrucList{1};
% strucs{strucCount} = temp1;
% strucsID(strucCount) = 1;
% sampleStruc(1) = strucCount;
% 
% for i = 2 : length(origStrucList)
%     
%     clear temp1; clear temp2; clear delim1; clear age1;
%     temp1 = origStrucList{i};
%     
%     R = 0; j = length(strucs) + 1;
%     while R == 0
%         clear temp2;
%         j = j - 1;
%         temp2 = strucs{j};
%         R = strcmp(temp1, temp2); 
%         if j == 1
%             break;
%         end        
%     end
%     
%     if R == 1
%         sampleStruc(i) = strucsID(j);
%     else
%         strucCount = strucCount + 1;
%         strucs{strucCount} = temp1;
%         strucsID(strucCount) = strucCount;
%         sampleStruc(i) = strucCount;
%     end
%     
% end
% clear i; clear j;
% save('files\sampleStruc.mat', 'sampleStruc');
% save('files\strucs.mat', 'strucs');
% save('files\strucsID.mat', 'strucsID');
%--------------------------------------------------------------------------

%%%--Load entrezID File----------------------------------------------------
% idFile = 'genes_matrix_csv\genes_rows_metadata.csv';
% entrezIDs = csvread(idFile, 1, 4);
% save('files\entrezIDs.mat', 'entrezIDs');
%--------------------------------------------------------------------------

%%%--Construct Data Matrix (rows:Age; columns:Structures)------------------
% load('files\origExpMat.mat');
% load('files\sampleAge');
% load('files\sampleStruc');
% load('files\ages.mat');
% load('files\strucs.mat');
% 
% dataMat = zeros(30,26);
% 
% for i = 1 : size(ages, 2)
%     
%     clear ageInd; clear strucsPerAge;
%     ageInds = find(sampleAge == i);
%     strucsPerAge = sampleStruc(ageInds);
%     for j = 1 : length(strucsPerAge)    
%         dataMat(i, strucsPerAge(j)) = dataMat(i, strucsPerAge(j)) + 1;
%     end
%     
% end
% figure, bar(dataMat', 'stacked'),
%  csvwrite('dataMat.csv', dataMat');
%--------------------------------------------------------------------------

%%%--Construct Gene Matrix (rows:Age; columns:Structures; depth:Probes)----
% load('files\origExpMat.mat');
% load('files\sampleAge');
% load('files\sampleStruc');
% load('files\ages.mat');
% load('files\strucs.mat');
% 
% ageIndInc = [3, 4, 5, 8, 9, 12, 13, 15, 16, 17, 18, 19, 21, 24, 25, 26, ...
%     27, 28, 29, 30];
% strucIndInc = [2, 3, 5, 6, 7, 8, 9, 15, 18, 19, 20, 21, 22, 23, 24];
% 
% probesNo = size(origExpMat, 1);
% geneMat = zeros(length(strucIndInc), length(ageIndInc), probesNo);
% 
% for i = 1 : length(ageIndInc)   
%     
%     clear currAge; clear cuurAgeSamplesInd; clear currAgeStrucs;
%     clear j;
%     currAge = ageIndInc(i);
%     cuurAgeSamplesInd = find(sampleAge == currAge);
%     currAgeStrucs = sampleStruc(cuurAgeSamplesInd);
%     
%     for j = 1 : length(strucIndInc)
%         
%         clear currStruc; clear test; clear currAgeInd; clear tempExp1;
%         clear tempExp2; clear gneMat;
%         currStruc = strucIndInc(j);
%         test = find(currAgeStrucs == currStruc);
%         
%         if ~isempty(test);
%             currAgeInd = cuurAgeSamplesInd(test);
%             tempExp1 = origExpMat(:, currAgeInd)';
%             tempExp2 = sum(tempExp1) ./ size(tempExp1,1);
%             geneMat(j, i, :) = tempExp2;
%         end
%         
%     end
%     
% end
% save('files\geneMat.mat', 'geneMat');
%--------------------------------------------------------------------------









