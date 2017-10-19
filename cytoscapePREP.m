% prepare ASD module genes for cytoscape

clear all;

C = [9 19 20 26];
for i = 1 : 4
    [num txt] = xlsread(['C:\Ahmed\Work\ASD work\Paper Drafts\hub genes results\original files\clust' num2str(C(i)) '.csv']);
    G1 = txt(:,1);
    G2 = txt(:,2);
    corr = num(:,2);
    clear txt; clear num;
    
    count = 0;
    gC = 1;
    count = count + 1;
    topConGenes(count) = G1(gC);
    count = count + 1;
    topConGenes(count) = G2(gC);
    gC = gC + 1;
    
    while count <= 50
        if isempty(find(ismember(topConGenes, G1(gC))) == 1);
            count = count + 1;
            topConGenes(count) = G1(gC);
        end
        
        if isempty(find(ismember(topConGenes, G2(gC))) == 1);
            count = count + 1;
            topConGenes(count) = G2(gC);
        end
        
        gC = gC + 1;
    end
    
    i1 = find(ismember(G1, topConGenes) == 1);
    i2 = find(ismember(G2, topConGenes) == 1);
    tcgInd = intersect(i1, i2);
    
    genes1 = G1(tcgInd);
    genes2 = G2(tcgInd);
    corrG = corr(tcgInd);
    xlswrite(['results\hClustResults_NEW\exp2\cytoscape\clust' num2str(C(i)) '.xls'], genes1, 1, 'A1');
    xlswrite(['results\hClustResults_NEW\exp2\cytoscape\clust' num2str(C(i)) '.xls'], genes2, 1, 'B1');
    xlswrite(['results\hClustResults_NEW\exp2\cytoscape\clust' num2str(C(i)) '.xls'], corrG, 1, 'C1');
    
    clear genes1; clear genes2; clear corrG; clear tcgInd; clear i1; clear i2;

end







