clc
clear all

Tfinal = 10400; t_rain = 18;

for  ggsm = 1:3
    twf = readmatrix(['GGSM_' num2str(ggsm-1) '/Total Water Footprint.csv']);

% ========================== DEFINE MODEL VARIABLES PT 1 =========================================================================

    for cc = 1:2
        amp = zeros(t_rain*12,7);
        pmp = zeros(Tfinal*12/52,7);
        pwp = zeros(Tfinal,7);
        pap = zeros(Tfinal/52,7);
        pwe = zeros(Tfinal,7);
        dev = zeros(12,7);
        exp = 135*ones(Tfinal,7);
        meo = zeros(Tfinal/52,7);
        men = zeros(Tfinal/52,7);
        ngo = zeros(Tfinal,7);
        ngn = zeros(Tfinal,7);
        
        
% ========================== IMPORT RAINFALL DATA =========================================================================
        
        for c = 1:7    
            
            rain = readmatrix(['rain_0' num2str(c) '.csv']); % 01 - global, 02 - africa, 03 - asia, 04 - europe, 05 - north america, 06 - oceania, 07 - south america
            sz = size(rain);
            months = sz(1);
            years = sz(2);
            states = 3; % dry, wet and flood
            area = [136 30 45 10 25 8 18]; % million sq.km
        
            exp(:,c) = exp(:,c).*(area(c)/area(1));
            ngn(:,c) = ones(Tfinal,1).*(1E5*area(c)/area(1));
            ngo(:,c) = ones(Tfinal,1).*(1E5*area(c)/area(1));
            
            
% ========================== DEFINE STATE THRESHOLDS =========================================================================
            
            % rmin and rmax are the dry-wet and wet-flood threshold vectors respectively
            % there are 12 entries, one for each month. mean and stddev were analyzed for each month over a period of 18 years. 
            % each value in rmin is (mean - stddev) while each value in rmax is (mean + stddev)
            
            rdry = mean(rain,2) - 2*std(rain,0,2);
            rflood = mean(rain,2) + 2*std(rain,0,2);
            rmin = min(mean(rain,2) - 3*std(rain,0,2));
            rmax = max(mean(rain,2) + 3*std(rain,0,2)); 
            
% ========================== DEFINE MODEL PARAMETERS =========================================================================
            
            expwet = 0.05; % Exploitable water recharge fraction in wet state.
            expflood = 0.01; % Exploitable water recharge fraction in flood state.
            climatechange = [0 1]; % toggle between 0 and 1 to activate climate change
            deltap = 0.0005; % Step change in transition probabilities. +ve deltap implies increasing frequency of floods and droughts. This captures the effects of climate change.
            deltar = 0.0001; % y-o-y change in upper and lower limits of precipitation. This captures extreme precipitation occurrences.
            
% ========================== DEFINE MODEL VARIABLES PT 2 =========================================================================
            
            rf2 = rain(:,1); % Convert matrix to column vector
            for k = 2:years
                    rf2 = vertcat(rf2, rain(:,k));
            end
            
            for k = 1:states
                for j = 1:states
                    nrain(k,j,:) = zeros(months,1); % transitions counter. 3x3 = 9 transitions. nrain(i,j,m) is the number of years in which month m is in state j and month m-1 in state i
                    Mrain(k,j,:) = zeros(months,1); % transition probability matrix. we convert the counter values into probabilities and store them in this.
                    Mrain_default(k,j,:) = zeros(months,1); % tpm copy. used to store original values of the tpm in case the main tpm is modified during simulations.
                end
            end
            
            Tfinal = 200*52; % simulation time frame of the GGSM, in weeks
            rf = zeros(Tfinal*months/52,1); % predicted monthly rainfall
            expr = zeros(Tfinal*months/52,1); % predicted monthly exploitable water recharge vector
            weekly_precip = zeros(Tfinal,1); % predicted weekly rainfall, input to the GGSM
            weekly_exp = zeros(Tfinal,1); % predicted weekly exploitable water recharge, input to the GGSM
            annmean = zeros(Tfinal/52,1);
            nsim = 50; % number of simulations of this model
            rain_mean = zeros(size(rf,1),nsim); % sample mean of all predicted rainfall
            rain_bestfit = zeros(size(rf,1),1);
            in = zeros(11,1); in_cum = zeros(12,1);
            d_exp_old = zeros(Tfinal,1);
            d_exp_new = zeros(Tfinal,1);
            
            
% ========================== UPDATE TRANSITIONS COUNTER =========================================================================
            
            
            for j = 1:years
                for k = 2:months
                    if rain(k,j) < rdry(k) && rain(k-1,j) < rdry(k-1)
                        nrain(1,1,k) = nrain(1,1,k) + 1;
                    else if rain(k,j) < rdry(k) && rain(k-1,j) > rdry(k-1) && rain(k-1,j) < rflood(k-1)
                            nrain(2,1,k) = nrain(2,1,k) + 1;
                    else if rain(k,j) > rdry(k) && rain(k,j) < rflood(k) && rain(k-1,j) < rdry(k-1)
                                nrain(1,2,k) = nrain(1,2,k) + 1;
                    else if rain(k,j) > rdry(k) && rain(k,j) < rflood(k) && rain(k-1,j) > rdry(k-1) && rain(k-1,j) < rflood(k-1)
                                    nrain(2,2,k) = nrain(2,2,k) + 1;
                    else if rain(k,j) > rflood(k) && rain(k-1,j) < rdry(k-1)
                                        nrain(1,3,k) = nrain(1,3,k) + 1;
                    else if rain(k,j) > rflood(k) && rain(k-1,j) > rdry(k-1) && rain(k-1,j) < rflood(k-1)
                                            nrain(2,3,k) = nrain(2,3,k) + 1;
                    else if rain(k,j) > rflood(k) && rain(k-1,j) > rflood(k-1)
                                                nrain(3,3,k) = nrain(3,3,k) + 1;
                    else if rain(k,j) < rdry(k) && rain(k-1,j) > rflood(k-1)
                                                    nrain(3,1,k) = nrain(3,1,k) + 1;
                    else if rain(k,j) > rdry(k) && rain(k,j) < rflood(k) && rain(k-1,j) > rflood(k-1)
                                                        nrain(3,2,k) = nrain(3,2,k) + 1;
                    end
                    end
                    end
                    end
                    end
                    end
                    end
                    end
                    end
                end
            
            %     Looping back to month 1 from month 12
            
                if rain(1,j) < rdry(1) && rain(12,j) < rdry(12)
                        nrain(1,1,1) = nrain(1,1,1) + 1;
                else if rain(1,j) < rdry(1) && rain(12,j) > rdry(12) && rain(12,j) < rflood(12)
                            nrain(2,1,1) = nrain(2,1,1) + 1;
                else if rain(1,j) > rdry(1) && rain(1,j) < rflood(1) && rain(12,j) < rdry(12)
                                nrain(1,2,1) = nrain(1,2,1) + 1;
                else if rain(1,j) > rdry(1) && rain(1,j) < rflood(1) && rain(12,j) > rdry(12) && rain(12,j) < rflood(12)
                                    nrain(2,2,1) = nrain(2,2,1) + 1;
                else if rain(1,j) > rflood(1) && rain(12,j) < rdry(12)
                                        nrain(1,3,1) = nrain(1,3,1) + 1;
                else if rain(1,j) > rflood(1) && rain(12,j) > rdry(12) && rain(12,j) < rflood(12)
                                            nrain(2,3,1) = nrain(2,3,1) + 1;
                else if rain(1,j) > rflood(1) && rain(12,j) > rflood(12)
                                                nrain(3,3,1) = nrain(3,3,1) + 1;
                else if rain(1,j) < rdry(1) && rain(12,j) > rflood(12)
                                                    nrain(3,1,1) = nrain(3,1,1) + 1;
                else if rain(1,j) > rdry(1) && rain(1,j) < rflood(1) && rain(12,j) > rflood(12)
                                                        nrain(3,2,1) = nrain(3,2,1) + 1;
                end
                end
                end
                end
                end
                end
                end
                end
                end 
            end
            
% ========================== UPDATE TRANSITION PROBABILITY MATRIX =========================================================================
            
            sum = zeros(states,months);
            % if month m-1 was in state k then month m has to be in one of the three states. so sum(k,m) is the total number of possibile transitions.
            
            for m = 1:months
                for j = 1:states
                    for k = 1:states
                        for l = 1:states
                            sum(k,m) = sum(k,m) + nrain(k,l,m); 
                        end
                        Mrain(k,j,m) = nrain(k,j,m)/sum(k,m); % the all important transition probability matrix
                    end
                end
            end
            
            Mrain_default = Mrain; % storing original tpm values
            
% ========================== PREDICT RAINFALL =========================================================================
            
            rf(1) = rain(1,1); % seed value 
            
            % The first for loop simulates the model as many times as desired
            
            for i = 1:nsim
                rmin = min(mean(rain,2) - 3*std(rain,0,2));
                rmax = max(mean(rain,2) + 3*std(rain,0,2));
                Mrain = Mrain_default; % Values are reset after every simulation
            
            % The second for loop generates the predicted rainfall vector for each simulation
            
                for m = 1:Tfinal*months/52 - 1
                    j = mod(m,months);
                    if j == 0
                        j = months;
                    end
                    rain_rnd1 = rand; % Random number in (0,1)
                    rain_rnd2 = rand; % Random number in (0,1)
            
            % This is important because it makes the model truly random. No 2 simulations will yield the same prediction.
            
            % 3 cases are considered based on the rainfall amount in the current timestep. This rainfall can either be in the dry,wet or flood state.  
            % To start the loop, initial rainfall amount is the same as the rainfall in month 1 and year 1 of the dataset. 
            % Given the current state, rain_rnd1 is then compared with the 3 possible transition probabilities. 
            % First step is to evaluate dry vs non-dry transition. If rain_rnd1 < dry transition probability then the next state will be dry, otherwise it will be non-dry. 
            % Second step is to repeat this exercise for the 2 non-dry states (wet and flood) by using rain_rnd2. 
            % Note that in step 2 we use relative instead of absolute probabilities because after eliminating the dry transition, total probability is no longer 1.
            
            % Once we determine the state of rainfall in the next timestep, we also need to determine the amount. Each state has a defined interval of rainfall amounts. 
            % Dry - (lowest recorded rainfall, dry-wet threshold),
            % Wet - (dry-wet threshold, wet-flood threshold), 
            % Flood - (wet-flood threshold, highest recorded rainfall)
            % Again, a random number is generated in the interval corresponding to the chosen state. 
            % This number becomes the amount of rainfall in the next timestep.
            
                        if rf(m) < rdry(j)
                            if rain_rnd1 < Mrain(1,1,j)
                                if j == months
                                    j = 0;
                                end
                                rf(m+1) = rmin + (rdry(j+1) - rmin)*rand;
                            else
                            if rain_rnd2 < Mrain(1,2,j)/(Mrain(1,2,j)+Mrain(1,3,j))
                                if j == months
                                    j = 0;
                                end
                                rf(m+1) = rdry(j+1) + (rflood(j+1) - rdry(j+1))*rand;
                            else
                                if j == months
                                    j = 0;
                                end
                                rf(m+1) = rflood(j+1) + (rmax - rflood(j+1))*rand;
                            end
                            end
                        end
                        if j == 0
                            j = months;
                        end
                        if rf(m) > rdry(j) && rf(m) < rflood(j)
                            if rain_rnd1 < Mrain(2,1,j)
                                if j == months
                                    j = 0;
                                end
                                rf(m+1) = rmin + (rdry(j+1) - rmin)*rand;
                            else
                            if rain_rnd2 < Mrain(2,2,j)/(Mrain(2,2,j)+Mrain(2,3,j))
                                if j == months
                                    j = 0;
                                end
                                rf(m+1) = rdry(j+1) + (rflood(j+1) - rdry(j+1))*rand;
                            else
                                if j == months
                                    j = 0;
                                end
                                rf(m+1) = rflood(j+1) + (rmax - rflood(j+1))*rand;
                            end
                            end
                        end
                        if j == 0
                            j = months;
                        end
                        if rf(m) > rflood(j)
                            if rain_rnd1 < Mrain(3,1,j)
                                if j == months
                                    j = 0;
                                end
                                rf(m+1) = rmin + (rdry(j+1) - rmin)*rand;
                            else
                            if rain_rnd2 < Mrain(3,2,j)/(Mrain(3,2,j)+Mrain(3,3,j))
                                if j == months
                                    j = 0;
                                end
                                rf(m+1) = rdry(j+1) + (rflood(j+1) - rdry(j+1))*rand;
                            else
                                if j == months
                                    j = 0;
                                end
                                rf(m+1) = rflood(j+1) + (rmax - rflood(j+1))*rand;
                            end
                            end
                        end
                        if j == 0
                            j = months;
                        end
                        for k = 1:states
                            if Mrain(k,2,j) > 0
                                    Mrain(k,1,j) = Mrain(k,1,j) + deltap*climatechange(cc);
                                    Mrain(k,3,j) = Mrain(k,3,j) + deltap*climatechange(cc);
                                    Mrain(k,2,j) = Mrain(k,2,j) - 2*deltap*climatechange(cc);
                            end
                        end
                        rmax = rmax*(1+deltar*climatechange(cc));
                        rmin = rmin*(1-deltar*climatechange(cc));
                            
                end
                
                rain_mean(:,i) = rf;
                rain_bestfit = rain_bestfit + rain_mean(:,i);
                
            end
                
            rain_bestfit = rain_bestfit./nsim;
        
            for j = 1:Tfinal
                weekly_precip(j) = rain_bestfit(ceil(j*months/52))*7*area(c); % convert monthly data to weekly
            end
        
            
            for i = 1:2400
                j = ceil(i/12);
                annmean(j) = annmean(j) + rain_bestfit(i);
            end
            annmean = annmean./12;
        
% ========================== COMPUTE EXPLOITABLE WATER RECHARGE =========================================================================
            
                
            % Only a fraction of rainfall is usable to humans as exploitable water.
            % Exploitable water recharge starts beyond the dry threshold of rainfall, after which it is directly proportional to rainfall. 
            % The proportionality constant is higher in the wet state because beyond the flood threshold, most rainfall is wasted and not stored.
                
            for j = 1:(Tfinal*months/52)-1
                k = mod(j,months);
                if k == 0
                    k = 12;
                end
                if rain_bestfit(j) < rmin
                    expr(j+1) = 0; 
                else
                if rain_bestfit(j) > rmin && rain_bestfit(j) < rmax
                    expr(j+1) = expwet*rain_bestfit(j); 
                else
                    expr(j+1) = expwet*rmax;
                end
                end
            end
            expr(1) = expr(2);
            
            for j = 1:Tfinal
                weekly_exp(j) = expr(ceil(j*months/52))*7*area(c); % convert monthly data to weekly
            end
            
        
% ========================== ERROR ANALYSIS =========================================================================
        
            bestfit_error = abs((rain_bestfit(1:216,1)-rf2)./rf2); % deviation from actual data
            
            for i = 1:size(rf2,1)
                for j = 1:11
                    if j < 11
                        if bestfit_error(i) < 0.1*j && bestfit_error(i) > 0.1*(j-1)
                            in(j) = in(j) + 1;
                        end
                    else
                        if bestfit_error(i) > 0.1*(j-1)
                            in(j) = in(j) + 1;
                        end
                    end
                end
            end
            
            for i = 1:11
                for j = 1:i
                    in_cum(i+1) = in_cum(i+1) + in(j); % cumulative probability distribution of error margins
                end
            end
            in_cum = in_cum./size(rf2,1); 

% ========================== INTEGRATE WITH GGSM =========================================================================

            amp(:,c) = rf2;
            pmp(:,c) = rain_bestfit;
            pwp(:,c) = weekly_precip;
            pap(:,c) = annmean;
            pwe(:,c) = weekly_exp;
            dev(:,c) = in_cum;
        
    
            for i = 1:Tfinal-1
                
                d_exp_new(i) = pwe(i,c) - twf(i,c);
                
                if d_exp_new(i) < 0
                    ngn(i+1,c) = ngn(i,c) + d_exp_new(i);
                    pwe(i,c) = 0;
                else
                    if d_exp_new(i) > 0
                        pwe(i,c) = d_exp_new(i);
                        ngn(i+1,c) = ngn(i,c);
                    end
        
                end
            end

            for i = 1:Tfinal-1

                d_exp_old(i) = exp(i,c) - twf(i,c);

                if d_exp_old(i) < 0
                    ngo(i+1,c) = ngo(i,c) + d_exp_old(i);
                    exp(i,c) = 0;
                else
                    if d_exp_old(i) > 0
                        exp(i,c) = d_exp_old(i);
                        ngo(i+1,c) = ngo(i,c);
                    end
                end
            end
        
            pwe(Tfinal,c) = pwe(Tfinal-1,c);
            exp(Tfinal,c) = exp(Tfinal-1,c);
            
            for i = 1:Tfinal
                j = ceil(i/52);
                meo(j,c) = meo(j,c) + exp(i,c);
                men(j,c) = men(j,c) + pwe(i,c);
            end
            
            meo(:,c) = meo(:,c)./52;
            men(:,c) = men(:,c)./52;   

            for i = 1:Tfinal
                if ngn(i,c) < 0
                    ngn(i:Tfinal,c) = 0.*ngn(i:Tfinal,c);
                    break
                end
            end
            for i = 1:Tfinal
                if ngo(i,c) < 0
                    ngo(i:Tfinal,c) = 0.*ngo(i:Tfinal,c);
                    break
                end
            end
        end

            

            
% ========================== OUTPUT FILES =========================================================================
    
        writematrix(amp,['Climate Change_' num2str(cc-1) '/actual_monthly_precipitation.csv']);
        writematrix(pmp,['Climate Change_' num2str(cc-1) '/predicted_monthly_precipitation.csv']);
        writematrix(pwp,['Climate Change_' num2str(cc-1) '/predicted_weekly_precipitation.csv']);
        writematrix(pap,['Climate Change_' num2str(cc-1) '/predicted_annual_precipitation.csv']);
        writematrix(pwe,['Climate Change_' num2str(cc-1) '/predicted_weekly_exploitablewater.csv']);
        writematrix(men,['Climate Change_' num2str(cc-1) '/GGSM_' num2str(ggsm-1) '/new_annual_exploitablewater.csv']);
        writematrix(meo,['Climate Change_' num2str(cc-1) '/GGSM_' num2str(ggsm-1) '/old_annual_exploitablewater.csv']);
        writematrix(ngo,['Climate Change_' num2str(cc-1) '/GGSM_' num2str(ggsm-1) '/old_nonrenewable_groundwater.csv']);
        writematrix(ngn,['Climate Change_' num2str(cc-1) '/GGSM_' num2str(ggsm-1) '/new_nonrenewable_groundwater.csv'])
        writematrix(dev,['Climate Change_' num2str(cc-1) '/deviation.csv']);
    end
end

   