nc = 7; % Number of compartments

% 1 - Atmosphere
% 2 - Nonexploitable freshwater
% 3 - Exploitable renewable water
% 4 - Non-renewable groundwater
% 5 - Human sectors
% 6 - Water recycling
% 7 - Non-freshwater reservoir

% Flow Compartments: w_flow(X,Y,i) = Flow from X to Y in timestep i
% Stock Compartments: w_stock(X,i) = Stock of X at the start of timestep i

year = zeros(Tfinal/52,1);
for i = 1:Tfinal/52
    year(i) = i;
end
for i = 1:nc
    w_stock(:,i) = zeros(Tfinal,1);
    net_flow(:,i) = zeros(Tfinal,1);
    for j = 1:nc
        w_flow(i,j,:) = zeros(Tfinal,1);
    end
end

% Stock Initialisation

w_stock(1,1) = 1E2; % cu.km
w_stock(1,2) = 1E5; % cu.km
w_stock(1,3) = 0; % cu.km
w_stock(1,4) = 1E5; % cu.km
w_stock(1,5) = 0; % cu.km
w_stock(1,6) = 0;
w_stock(1,7) = 1E9; % cu.km
w_total = 0;

for i = 1:nc
    w_total = w_total + w_stock(1,i);
end

BWF = zeros(Tfinal,1); GWF = zeros(Tfinal,1); REC = zeros(Tfinal,1); EXP = zeros(Tfinal,1);
NRG = zeros(Tfinal,1); CR_NRG = zeros(Tfinal/52,1); CR_EXP = zeros(Tfinal,1);


% Flow Balance
w = zeros(nc,nc); % parameter matrix
w(1,5) = 715; % GWF m3/ton
w(5,6) = 0.2; % Recycle split
w(5,7) = 1-w(5,6);
w(6,3) = 0.2; % Recycle fraction
w(6,2) = 1-w(6,3);

w_flow(1,2,i) = 46339/52; % weekly nonexploitable renewable water
w_flow(1,3,i) = 7890/52; % weekly exploitable renewable water
w_flow(1,5,i) = w(1,4)*P1(i); % GWF
w_flow(1,7,i) = 1E5; 

w_flow(2,1,i) = 0;
w_flow(2,3,i) = 0; 
w_flow(2,7,i) = 0; 

w_flow(3,1,i) = 0;
w_flow(4,5,i) = 300/52; % weekly non renewable groundwater withdrawal
w_flow(3,5,i) = dem_W_total(i) - w_flow(8,5,i);

w_flow(5,1,i) = 0;
w_flow(5,6,i) = w(5,6)*(w_flow(3,5,i) + w_flow(4,5,i) + w_flow(1,5));
w_flow(5,7,i) = w(5,7)*(w_flow(3,5,i) + w_flow(4,5,i) + w_flow(1,5));

if i > 2
    w_flow(6,3,i) = w(6,3)*w_flow(5,6,i-1);
    w_flow(6,2,i) = w(6,2)*w_flow(5,6,i-1);
end

w_flow(7,1,i) = w_flow(1,7,i) + w_flow(2,7,i) + w_flow(5,7,i);


% Stock Balance
net_stock = 0; 
for j = 1:nc
    for k = 1:nc
        net_flow(i,j) = net_flow(i,j) + w_flow(k,j,i) - w_flow(j,k,i);
    end
    if j == 2||j==3
        w_stock(i,j) = w_stock(i,j) + net_flow(i,j);
    end
    net_stock = net_stock + w_stock(i+1,j);
end

if w_total ~= net_stock
    for j = 1:nc
        w_stock(i+1,j) = w_stock(i+1,j) + (w_stock(i,j)/w_total)*(w_total - net_stock);
    end
end

BWF(i) = w_flow(3,5,i) + w_flow(4,5,i);
GWF(i) = w_flow(1,5,i);
REC(i) = w_flow(6,3,i);
NRG(i) = w_stock(i,4);
EXP(i) = w_stock(i,3);
CR_EXP(i) = 100*w_flow(3,5,i)/w_flow(1,3,i);
if i/52 == ceil(i/52)
    for j = i-51:i
        CR_NRG(i/52) = CR_NRG(i/52) + 100*w_flow(4,5,j)/w_stock(i-51,4);
    end
end

csvwrite('Greenwater footprint.csv', GWF);
csvwrite('Bluewater footprint.csv', BWF);
csvwrite('Recycled water.csv', REC);
csvwrite('Non-renewable groundwater reserves.csv', NRG);
csvwrite('Renewable exploitable water.csv', EXP);
csvwrite('Criticality ratio (renewable).csv',CR_EXP);
csvwrite('Criticality ratio (non-renewable).csv',CR_NRG);