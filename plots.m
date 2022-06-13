clc
clear all

% =============================== IMPORT DATA =========================================================================

years_precip = 18;
years_ggsm = 200;
years = linspace(2000,2017,years_precip*12);
x = linspace(0,1.1,12);
t = linspace(2000,2199,200);
w = linspace(1957,2156,years_ggsm*52);
y = linspace(1957,2156,years_ggsm);
pap_copy = zeros(years_ggsm,7,2);

cc = 1;

amp = readmatrix(num2str(cc-1,"Climate Change_%01d/actual_monthly_precipitation.csv"));
pmp = readmatrix(num2str(cc-1,"Climate Change_%01d/predicted_monthly_precipitation.csv"));
dev = readmatrix(num2str(cc-1,"Climate Change_%01d/deviation.csv"));
    
    
% =============================== ACTUAL VS PREDICTED =========================================================================
    
figure(1)

for i = 1:7
    subplot(3,3,i)
    plot(years,amp(:,i),'b',years,pmp(1:216,i),'r')
    xlabel('Year')
    ylabel('mm/day')
%     legend('Actual','Predicted')
end
saveas(figure(1),['Plots/Climate Change_' num2str(cc-1) '/Actual v Predicted.jpg'])
    
% =============================== DEVIATION =========================================================================
    
figure(2)

for i = 1:7
    plot(x,dev(:,i),'*')
    hold on
end
xlabel('Fractional Error')
ylabel('Cumulative Probability Distribution')
legend('Global','Africa','Asia','Europe','North America','Oceania','South America')
hold off
saveas(figure(2),['Plots/Climate Change_' num2str(cc-1) '/Deviation.jpg'])


    
        
% =============================== GGSM - OLD VS NEW =========================================================================
    
for ggsm = 1:3

    ngo = readmatrix(['Climate Change_' num2str(cc-1) '/GGSM_' num2str(ggsm-1) '/old_nonrenewable_groundwater.csv']);
    ngn = readmatrix(['Climate Change_' num2str(cc-1) '/GGSM_' num2str(ggsm-1) '/new_nonrenewable_groundwater.csv']);
    meo = readmatrix(['Climate Change_' num2str(cc-1) '/GGSM_' num2str(ggsm-1) '/old_annual_exploitablewater.csv']);
    men = readmatrix(['Climate Change_' num2str(cc-1) '/GGSM_' num2str(ggsm-1) '/new_annual_exploitablewater.csv']);
    
    figure(3 + 2*(ggsm-1))
    for i = 1:7
        subplot(3,3,i)
        plot(w,ngo(:,i),'.b',w,ngn(:,i),'.r')
        xlabel('Year')
        ylabel('x 10^9 m^3')
%             legend('Without Precip. Model','With Precip. Model')
        hold on
    end
    saveas(figure(3 + 2*(ggsm-1)),['Plots/Climate Change_' num2str(cc-1) '/GGSM_' num2str(ggsm-1) '/Nonrenewable Groundwater.jpg'])

    figure(4 + 2*(ggsm-1))
    for i = 1:7
        subplot(3,3,i)
        plot(y,meo(:,i),'.b',y,men(:,i),'.r')
        xlabel('Year')
        ylabel('x 10^9 m^3/ week')
%             legend('Without Precip. Model','With Precip. Model')
        hold on
    end
    saveas(figure(4 + 2*(ggsm-1)),['Plots/Climate Change_' num2str(cc-1) '/GGSM_' num2str(ggsm-1) '/Exploitable Water.jpg'])

end

% =============================== GGSM - BASE VS POPEX VS CONSINC =========================================================================

figure(9)

for i = 1:7
    for ggsm = 1:3

        col = {'-b', '-r', '-g'};
    
        men = readmatrix(['Climate Change_' num2str(cc-1) '/GGSM_' num2str(ggsm-1) '/new_annual_exploitablewater.csv']);
    
        subplot(3,3,i)
        plot(y,men(:,i), col{ggsm} )
        xlabel('Year')
        ylabel('x 10^9 m^3/ week')
        hold on
    end
    
end

saveas(figure(9),['Plots/Climate Change_' num2str(cc-1) '/Base v PopEx v ConsInc/Exploitable Water.jpg'])

figure(10)

for i = 1:7

    for ggsm = 1:3
    
        col = {'-b', '-r', '-g'};

        ngn = readmatrix(['Climate Change_' num2str(cc-1) '/GGSM_' num2str(ggsm-1) '/new_nonrenewable_groundwater.csv']);
    
        subplot(3,3,i)
        plot(w,ngn(:,i), col{ggsm})
        xlabel('Year')
        ylabel('x 10^9 m^3/ week')
        hold on
    end
end
hold off
saveas(figure(10),['Plots/Climate Change_' num2str(cc-1) '/Base v PopEx v ConsInc/Nonrenewable Groundwater.jpg'])


% =============================== CLIMATE CHANGE =========================================================================

figure(11)

for cc = 1:2
    pap = readmatrix(num2str(cc-1,"Climate Change_%01d/predicted_annual_precipitation.csv"));
    plot(t,pap(:,1),'.')
    hold on
end

hold off
xlabel('Year')
ylabel('mm/day')
legend('Base','Climate Change')
saveas(figure(11),"Plots/Base v Climate Change.jpg")  








