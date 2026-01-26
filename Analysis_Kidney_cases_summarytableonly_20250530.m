

%% Cacluate corr between titan and Orion

sum1 = tableCombine;

listT = 1:768;
listO = 769:848;

corrTitanOrion = zeros(length(listO),length(listT));

for i = 1:length(listO)
    for j=1:length(listT)
        corrTitanOrion(i,j)=corr(sum1{:,listO(i)},sum1{:,listT(j)});
    end
end
