clc
close all

newcolors = [0.8500 0.3250 0.0980
             0.9290 0.6940 0.1250
             0.4940 0.1840 0.5560];
colororder(newcolors);

cr_exp = readmatrix("CR_exp.xlsx");
cr_nrg = readmatrix("CR_nrg.xlsx");
bwf = readmatrix("BWF.xlsx");
gwf = readmatrix("GWF.xlsx");
nrg = readmatrix("NRG.xlsx");
rec = readmatrix("REC.xlsx");

y = linspace(1957,2157,10400);

figure(1)

subplot(1,2,1)

for j = 2:4
    plot(y,bwf(2:10401,j),'.')
    hold on
    grid on
end
xlabel('Year')
ylabel('x10^9 m^3/week')
title('Bluewater footprint')
legend('Base','PopEx','ConsInc')
legend('Location','southeast')

subplot(1,2,2)

for j = 2:4
    plot(y,cr_exp(:,j),'.')
    hold on
    grid on
end
xlabel('Year')
ylabel('Criticality Ratio (%)')
title('Exploitable water stress')
legend('Base','PopEx','ConsInc')
legend('Location','southeast')

saveas(figure(1), 'CR_v_BWF.jpg')

figure(2)

subplot(1,2,1)

for j = 2:4
    plot(y,nrg(2:10401,j),'.')
    hold on
    grid on
end
xlabel('Year')
ylabel('x10^9 m^3/week')
title('Nonrenewable groundwater')
legend('Base','PopEx','ConsInc')

subplot(1,2,2)

for j = 2:4
    plot(y,cr_nrg(2:10401,j),'.')
    hold on
    grid on
end
xlabel('Year')
ylabel('x10^9 m^3')
title('Nonrenewable groundwater stress')
legend('Base','PopEx','ConsInc')

saveas(figure(2), 'CR_v_NRG.jpg')

figure(3)

for j = 2:4
    plot(y,gwf(2:10401,j),'.')
    hold on
    grid on
end
xlabel('Year')
ylabel('x10^9 m^3/week')
% title('Greenwater footprint')
legend('Base','PopEx','ConsInc')
legend('Location','southeast')

saveas(figure(3), 'GWF.jpg')

figure(4)

for j = 2:3
    plot(y,rec(2:10401,j),'.')
    hold on
    grid on
end
xlabel('Year')
ylabel('x10^9 m^3/week')
% title('Recycling analysis')
legend('Recycling disabled','Recycling enabled')

saveas(figure(4),'REC.jpg')