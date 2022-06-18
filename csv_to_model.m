clc
clear all

ra1 = readmatrix('rain_asia1.csv');
ra1 = transpose(ra1);
ra2 = readmatrix('rain_asia2.csv');
ra2 = transpose(ra2);
ra3 = readmatrix('rain_asia3.csv');
ra3 = transpose(ra3);
ra4 = readmatrix('rain_asia4.csv');
ra4 = transpose(ra4);
raf1 = readmatrix('rain_africa1.csv');
raf1 = transpose(raf1);
raf2 = readmatrix('rain_africa2.csv');
raf2 = transpose(raf2);
re1 = readmatrix('rain_europe1.csv');
re1 = transpose(re1);
re2 = readmatrix('rain_europe2.csv');
re2 = transpose(re2);
rna1 = readmatrix('rain_na1.csv');
rna1 = transpose(rna1);
rna2 = readmatrix('rain_na2.csv');
rna2 = transpose(rna2);
rna3 = readmatrix('rain_na3.csv');
rna3 = transpose(rna3);
rna4 = readmatrix('rain_na4.csv');
rna4 = transpose(rna4);
rna5 = readmatrix('rain_na5.csv');
rna5 = transpose(rna5);
rna6 = readmatrix('rain_na6.csv');
rna6 = transpose(rna6);
rsa1 = readmatrix('rain_sa1.csv');
rsa1 = transpose(rsa1);
rsa2 = readmatrix('rain_sa2.csv');
rsa2 = transpose(rsa2);
rsa3 = readmatrix('rain_sa3.csv');
rsa3 = transpose(rsa3);
ro = readmatrix('rain_o.csv');
ro = transpose(ro);

rain_asia = zeros(12,18);
rain_africa = zeros(12,18);
rain_europe = zeros(12,18);
rain_oceania = zeros(12,18);
rain_NA = zeros(12,18);
rain_SA = zeros(12,18);
rain_global = zeros(12,18);

for i = 2:13
    for j = 1:18
        rain_asia(i-1,j) = (ra1(i,j) + ra2(i,j) + ra3(i,j) + ra4(i,j))/7000;
        rain_africa(i-1,j) = (raf1(i,j) + raf2(i,j))/3000;
        rain_europe(i-1,j) = (re1(i,j) + re2(i,j))/1700;
        rain_NA(i-1,j) = (rna1(i,j) + rna2(i,j) + rna3(i,j) + rna4(i,j) + rna5(i,j) + rna6(i,j))/5000;
        rain_SA(i-1,j) = (rsa1(i,j) + rsa2(i,j) + rsa3(i,j))/1400;
        rain_oceania(i-1,j) = ro(i,j)/875;
        rain_global(i-1,j) = (rain_asia(i-1,j)*45 + rain_africa(i-1,j)*30 + rain_europe(i-1,j)*10 + rain_NA(i-1,j)*25 + rain_SA(i-1,j)*18 + rain_oceania(i-1,j)*8)/136;
    end
end

writematrix(rain_asia,'rain_03.csv')
writematrix(rain_africa,'rain_02.csv')
writematrix(rain_oceania,'rain_06.csv')
writematrix(rain_europe,'rain_04.csv')
writematrix(rain_NA,'rain_05.csv')
writematrix(rain_SA,'rain_07.csv')
writematrix(rain_global,'rain_01.csv')