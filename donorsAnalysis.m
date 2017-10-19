%%%-- Donors Analysis ------------------
load('files\origExpMat.mat');
load('files\sampleAge');
load('files\sampleStruc');
load('files\ages.mat');
load('files\strucs.mat');
load('files\donorIDs.mat');

% Get Donors Names
[num txt] = xlsread('files\DataFiles\genes_matrix_csv\genes_columns_metadata.xls', 'B2:C580');
[num2 ag] = xlsread('files\DataFiles\genes_matrix_csv\genes_columns_metadata.xls', 'D2:E580');
clear num2;

uniqueDonorIDs = unique(donorIDs);

for i = 1 : length(uniqueDonorIDs)
    
    clear dIND; clear strucsPerDonor; clear tempAge;
    dIND = find(donorIDs == uniqueDonorIDs(i));
    strucsPerDonor = sampleStruc(dIND);
    tempAge = sampleAge(dIND);
    donorAge(i) = tempAge(1);
    
    for j = 1 : length(strucsPerDonor)    
        donorsData(i, strucsPerDonor(j)) = 1;
    end
    
    dNameIND = txt(find(num == uniqueDonorIDs(i)), 2);
    dNames{i} = dNameIND{1};
    
    ageIND = ag(find(num == uniqueDonorIDs(i)), 1);
    ageNames{i} = ageIND{1};
    
    genderIND = ag(find(num == uniqueDonorIDs(i)), 2);
    genderNames{i} = genderIND{1};
    
end

% xlswrite('files\dataFiles\donorsData.xls', donorsData', 1, 'B4');
% xlswrite('files\dataFiles\donorsData.xls', ageNames, 1, 'B1');
% xlswrite('files\dataFiles\donorsData.xls', genderNames, 1, 'B2');
% xlswrite('files\dataFiles\donorsData.xls', dNames, 1, 'B3');

xlswrite('files\dataFiles\donorsNames&IDs.xls', uniqueDonorIDs, 1, 'A1');
xlswrite('files\dataFiles\donorsNames&IDs.xls', dNames', 1, 'B1');
xlswrite('files\dataFiles\donorsNames&IDs.xls', ageNames', 1, 'C1');
xlswrite('files\dataFiles\donorsNames&IDs.xls', donorAge', 1, 'D1');
%--------------------------------------------------------------------------