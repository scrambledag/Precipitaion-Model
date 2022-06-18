clc
close all
clear all

% =============================== IMPORT DATA =========================================================================

t_precip = 18;
t_ggsm = 200;
years = linspace(2000,2017,t_precip*12);
x = linspace(0,1.1,12);
t = linspace(2000,2199,200);
w = linspace(1957,2156,t_ggsm*52);
y = linspace(1957,2156,t_ggsm);
% col = {'#7E2F8E', '.r', '.m'}; 
newcolors = [0.8500 0.3250 0.0980
             0.9290 0.6940 0.1250
             0.4940 0.1840 0.5560
             0.4660 0.6740 0.1880
             0.3010 0.7450 0.9330];
         
colororder(newcolors)

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
saveas(figure(1),['Plots/Actual v Predicted.jpg'])
    
% =============================== DEVIATION =========================================================================
    
figure(2)

for i = 1:7
    plot(x,dev(:,i),'*')
    hold on
end
grid on
grid minor
xlabel('Fractional Error')
ylabel('Cumulative Probability Distribution')
legend('Global','Africa','Asia','Europe','North America','Oceania','South America')
legend('Location','southeast')
hold off
saveas(figure(2),['Plots/Deviation.jpg'])


    
        
% =============================== GGSM - OLD VS NEW =========================================================================
    
for ggsm = 1:3

    figure(3 + 2*(ggsm-1))

    for i = 1:7
        for cc = 1:2
            ngo = readmatrix(['Climate Change_' num2str(cc-1) '/GGSM_' num2str(ggsm-1) '/old_nonrenewable_groundwater.csv']);
            ngn = readmatrix(['Climate Change_' num2str(cc-1) '/GGSM_' num2str(ggsm-1) '/new_nonrenewable_groundwater.csv']);
            
            subplot(3,3,i)

            plot(w,ngo(:,i),'.',w,ngn(:,i),'.')
            grid on
            hold on
        end
        xlabel('Year')
        ylabel('x 10^9 m^3')
    end
    
    saveas(figure(3 + 2*(ggsm-1)),['Plots/GGSM_' num2str(ggsm-1) '/Nonrenewable Groundwater.jpg'])
end

for ggsm = 1:3
    
    figure(4 + 2*(ggsm-1))

    for i = 1:7
        for cc = 1:2

            meo = readmatrix(['Climate Change_' num2str(cc-1) '/GGSM_' num2str(ggsm-1) '/old_annual_exploitablewater.csv']);
            men = readmatrix(['Climate Change_' num2str(cc-1) '/GGSM_' num2str(ggsm-1) '/new_annual_exploitablewater.csv']);
            subplot(3,3,i)
            plot(y,meo(:,i),'.',y,men(:,i),'.')
            grid on
            hold on
        end
        xlabel('Year')
        ylabel('x 10^9 m^3/ week')
    end

    saveas(figure(4 + 2*(ggsm-1)),['Plots/GGSM_' num2str(ggsm-1) '/Exploitable Water.jpg'])

end

% =============================== GGSM - BASE VS POPEX VS CONSINC =========================================================================

figure(9)

for i = 1:7
    for ggsm = 1:3
        men = readmatrix(['Climate Change_' num2str(cc-1) '/GGSM_' num2str(ggsm-1) '/new_annual_exploitablewater.csv']);
    
        subplot(3,3,i)
        plot(y,men(:,i), '.' )
        xlabel('Year')
        ylabel('x 10^9 m^3/ week')
        grid on
        grid minor
        hold on
    end 
end

hold off
saveas(figure(9),['Plots/Base v PopEx v ConsInc/Exploitable Water.jpg'])

figure(10)

for i = 1:7

    for ggsm = 1:3

        ngn = readmatrix(['Climate Change_' num2str(cc-1) '/GGSM_' num2str(ggsm-1) '/new_nonrenewable_groundwater.csv']);
    
        subplot(3,3,i)
        plot(w,ngn(:,i), '.')
        xlabel('Year')
        ylabel('x 10^9 m^3/ week')
        grid on
        grid minor
        hold on
    end
end
hold off
saveas(figure(10),['Plots/Base v PopEx v ConsInc/Nonrenewable Groundwater.jpg'])


% =============================== CLIMATE CHANGE =========================================================================

figure(11)

for cc = 1:2
    pap = readmatrix(num2str(cc-1,"Climate Change_%01d/predicted_annual_precipitation.csv"));
    plot(t,pap(:,1),'.')
    grid on
    hold on
end

hold off
xlabel('Year')
ylabel('mm/day')
legend('Base','Climate Change')
legend('Location','northwest')
saveas(figure(11),"Plots/Base v Climate Change.jpg")  


