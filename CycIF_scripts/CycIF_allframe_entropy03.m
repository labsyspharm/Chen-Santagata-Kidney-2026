%%  Calculation entropy for all frames
% Jerry lin 2017/1/2


%% Initilization

total_col = 17;
total_row = 13;
cell_limit = 1000;
myName = 'EGFR';


allentropy = zeros(total_row,total_col);
entropytable = NaN(total_col*total_row,6);

for co = 1:total_col;
  for ro = 1:total_row;
    i = co+(ro-1)*total_col;
    
    data1 = cleardata(cleardata.frame == i,:);
    cellno = size(data1,1);
    entropytable(i,1) = cellno;
    entropytable(i,2) = 0;
    entropytable(i,3) = 0;
    entropytable(i,4) = 0;
    entropytable(i,5) = co;
    entropytable(i,6) = ro;
        
    if(cellno > cell_limit)
        eval(strcat('data2 = data1.',myName,';'));
        ne1 = nentropy(data2,'shannon');
        sample1 = datasample(cleardata,cellno);
        eval(strcat('sample2 = sample1.',myName,';'));
        ne2 = nentropy(sample2,'shannon');
        entropytable(i,2) = ne1;
        entropytable(i,3) = ne2;
        entropytable(i,4) = ne1/ne2;
        %display(ne1);
        %display(ne2);
        allentropy(data1.ROW,data1.COL) = ne1 / ne2;
    end
  end
end

entropytable = array2table(entropytable,'VariableNames',{'cellno','ne1','ne2','ne1_ne2','COL','ROW'});
eval(strcat('entropytable_',myName,' = entropytable;'));
eval(strcat('allentropy_',myName,' = allentropy;'));
eval('clear allentropy entropytable;');





