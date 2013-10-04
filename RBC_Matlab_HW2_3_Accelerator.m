%% Basic RBC model with partial depreciation and endogenous labor

% Original by Jesus Fernandez-Villaverde at Haverford, July 31, 2013
% Modifiied by Kory Kantenga at University of Pennsylvania, Oct 3, 2013

%% 0. Housekeeping

clear all
close all
clc

tic

%%  1. Calibration

aalpha = 1/3;     % Elasticity of output w.r.t. capital
bbeta  = 0.95;    % Discount factor
ddelta = 0.09;    % Rate of depreciation - added 9/10/13

% Productivity values
vProductivity = [0.9792; 0.9896; 1.0000; 1.0106; 1.0212]';

% Transition matrix
mTransition   = [0.9727, 0.0273, 0.0000, 0.0000, 0.0000;
    0.0041, 0.9806, 0.0153, 0.0000, 0.0000;
    0.0000, 0.0082, 0.9837, 0.0082, 0.0000;
    0.0000, 0.0000, 0.0153, 0.9806, 0.0041;
    0.0000, 0.0000, 0.0000, 0.0273, 0.9727];

%% 2a. Steady State

laborSteadyState = 1/3; %added 9/10/13
capitalSteadyState = laborSteadyState*((((1/bbeta)+ddelta-1)/aalpha)^(1/(aalpha-1))); %modified 9/10/13
outputSteadyState = (capitalSteadyState^aalpha)*(laborSteadyState^(1-aalpha)); %modified 9/10/13
consumptionSteadyState = outputSteadyState-ddelta*capitalSteadyState; %modified 9/10/13

%% 2b. Calibrate Disutility of Labor
ppsi = (1/consumptionSteadyState)*3*(1-aalpha)*((((1/bbeta)+ ddelta - 1)/aalpha)^(aalpha/(aalpha-1))); %disutility of labor - added 9/12/13
fprintf('Calibration\n');
fprintf(' Psi = %2.6f\n', ppsi);
fprintf('Steady State Values\n');
fprintf(' Output = %2.6f, Labor = %2.6f\n',outputSteadyState, laborSteadyState);
fprintf(' Capital = %2.6f, Consumption = %2.6f\n',  capitalSteadyState, consumptionSteadyState);
fprintf('\n');

% We generate the grid of capital
vGridCapital = 0.5*capitalSteadyState:0.00006:1.5*capitalSteadyState; %FIXME: change step to 0.00006
vGridLabor = 0.5*laborSteadyState:0.005:1.5*laborSteadyState; %TODO: change step?

nGridLabor = length(vGridLabor); %added 9/12/13
nGridCapital = length(vGridCapital);
nGridProductivity = length(vProductivity);

%% 3. Required matrices and vectors

mOutput           = zeros(nGridCapital,nGridProductivity,nGridLabor); %modified 9/12/13
mLabor            = zeros(nGridCapital,nGridProductivity,2);
mLaborIndex       = zeros(nGridCapital,nGridProductivity,2);
mValueFunction    = zeros(nGridCapital,nGridProductivity,2);
mValueFunctionNew = zeros(nGridCapital,nGridProductivity,2);
mPolicyFunction   = zeros(nGridCapital,nGridProductivity,2);
mPolicyIndex       = zeros(nGridCapital,nGridProductivity,2);
mConsumptionFunction = zeros(nGridCapital,nGridProductivity,2);
mEulerError        = zeros(nGridCapital,nGridProductivity,2);
expectedValueFunction = zeros(nGridCapital,nGridProductivity);
vTime = zeros(2,1);

%% 4. We pre-build output for each point in the grid

for j = 1:nGridLabor
    
    mOutput(:,:,j) = (vGridCapital'.^aalpha)*vProductivity*(vGridLabor(j)^(1-aalpha));
    
end

%% 5. Main iteration

maxDifference = 10.0;
maxIteration = 500; %FIXME
tolerance = 0.0000001; %FIXME
iteration = 0;

%Standard VFI
fprintf('Standard VFI\n')
nTime = tic;
while (maxDifference>tolerance)&&(iteration<maxIteration)  
    
    expectedValueFunction = mValueFunction(:,:,1)*mTransition';
    
    for nProductivity = 1:nGridProductivity
        
        % We start from previous choice (monotonicity of policy function)
        gridCapitalNextPeriod = 1;
        
        for nCapital = 1:nGridCapital
                        
            valueHighSoFar1 = -1000.0;
            capitalChoice  = vGridCapital(1);
            
            for nCapitalNextPeriod = gridCapitalNextPeriod:nGridCapital
                
                valueHighSoFar2 = -1000.0;
                
                for nLabor = 1:nGridLabor
                
                consumption = mOutput(nCapital,nProductivity,nLabor)+(1-ddelta)*vGridCapital(nCapital)-vGridCapital(nCapitalNextPeriod);
                valueProvisional = (1-bbeta)*(log(consumption)-ppsi*(vGridLabor(nLabor)^2)/2)+bbeta*expectedValueFunction(nCapitalNextPeriod,nProductivity);
                
                   if (valueProvisional>valueHighSoFar2)
                    valueHighSoFar2 = valueProvisional;
                    laborChoice = vGridLabor(nLabor);
                    gridLabor = nLabor;
                   else
                        break; %We break when we have achieved the max
                    end
                end
               
                
                if (valueProvisional>valueHighSoFar1)
                    valueHighSoFar1 = valueProvisional;
                    capitalChoice = vGridCapital(nCapitalNextPeriod);
                    gridCapitalNextPeriod = nCapitalNextPeriod;
                else
                   break; % We break when we have achieved the max
                end    
            end
            
            mConsumptionFunction(nCapital,nProductivity,1) = consumption;
            mValueFunctionNew(nCapital,nProductivity,1) = valueHighSoFar1;
            mPolicyFunction(nCapital,nProductivity,1) = capitalChoice;
            mLabor(nCapital,nProductivity,1) = laborChoice; %Labor function choosing next period capital optimally
            
        end %end for capital grid loop
        
    end %end of productivity grid loop
    
    maxDifference = max(max(abs(mValueFunctionNew(:,:,1)-mValueFunction(:,:,1))));
    mValueFunction(:,:,1) = mValueFunctionNew(:,:,1);
    
    iteration = iteration+1;
    if (mod(iteration,10)==0 || iteration ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference); 
    end
           
end
vTime(1,1) = toc(nTime);
fprintf('Time for VFI = %2.8f', vTime(1,1));
    fprintf(' seconds\n');
    fprintf('\n')

fprintf('\n')

%Reset Iteration Counter & Max Difference
iteration = 0;
maxDifference = 10.0;

%Accelerated VFI
nTime = tic;
fprintf('Accelerated VFI (Skip Max Operator 9 of 10 times)\n')
while (maxDifference>tolerance)&&(iteration<maxIteration) 
    
    expectedValueFunction = mValueFunction(:,:,2)*mTransition';
        
    if (mod(iteration,10)==0)||(iteration==0)||(iteration==1) %Maximize and Bellman Operator
        for nProductivity = 1:nGridProductivity
            
            % We start from previous choice (monotonicity of policy function)
            gridCapitalNextPeriod = 1;
            
            for nCapital = 1:nGridCapital
                
                valueHighSoFar1 = -1000.0;
                capitalChoice  = vGridCapital(1);
                
                for nCapitalNextPeriod = gridCapitalNextPeriod:nGridCapital
                    
                    valueHighSoFar2 = -1000.0;
                    
                    for nLabor = 1:nGridLabor
                        
                        consumption = mOutput(nCapital,nProductivity,nLabor)+...
                            (1-ddelta)*vGridCapital(nCapital)-vGridCapital(nCapitalNextPeriod);
                        valueProvisional = (1-bbeta)*(log(consumption)-ppsi*...
                            (vGridLabor(nLabor)^2)/2)+bbeta*...
                            expectedValueFunction(nCapitalNextPeriod,nProductivity);
                        
                        if (valueProvisional>valueHighSoFar2)
                            valueHighSoFar2 = valueProvisional;
                            laborChoice = vGridLabor(nLabor);
                            gridLabor = nLabor;
                        else
                            break; %We break when we have achieved the max
                        end
                    end %end of labor choice loop
                    
                    
                    if (valueProvisional>valueHighSoFar1)
                        valueHighSoFar1 = valueProvisional;
                        capitalChoice = vGridCapital(nCapitalNextPeriod);
                        gridCapitalNextPeriod = nCapitalNextPeriod;
                    else
                        break; % We break when we have achieved the max
                    end
                end %end of k' loop
                
                mConsumptionFunction(nCapital,nProductivity,2) = consumption;
                mValueFunctionNew(nCapital,nProductivity,2) = valueHighSoFar1;
                mPolicyFunction(nCapital,nProductivity,2) = capitalChoice;
                mPolicyIndex(nCapital,nProductivity,2) = nCapitalNextPeriod;
                mLabor(nCapital,nProductivity,2) = laborChoice; %Labor function choosing next period capital optimally
                mLaborIndex(nCapital,nProductivity,2) = nLabor; %Indexes from labor function
                
            end %end for capital grid loop
            
        end %end of productivity grid loop
        
    else %Bellman Operator only
        for nProductivity = 1:nGridProductivity
            for nCapital = 1:nGridCapital
                consumption = mOutput(nCapital,nProductivity,...
                    mLaborIndex(nCapital,nProductivity,2))...
                    +(1-ddelta)*vGridCapital(nCapital)-...
                    mPolicyFunction(nCapital,nProductivity,2);
                
                mConsumptionFunction(nCapital,nProductivity,2) = consumption;       
                mValueFunctionNew(nCapital,nProductivity,2) = (1-bbeta)*...
                    (log(consumption)-ppsi*(mLabor(nCapital,nProductivity,2)^2)/2)+...
                    bbeta*expectedValueFunction(mPolicyIndex(nCapital,nProductivity,2)...
                    ,nProductivity);
            end
        end
    end
    
    maxDifference = max(max(abs(mValueFunctionNew(:,:,2)-mValueFunction(:,:,2))));
    mValueFunction(:,:,2) = mValueFunctionNew(:,:,2);
    
    iteration = iteration+1;
    if (mod(iteration,10)==0 || iteration ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference);
    end
    
end
vTime(2,1) = toc(nTime);
fprintf('Time for VFI = %2.8f', vTime(2,1));
    fprintf(' seconds\n');
    fprintf('\n')

%% 6. Calculate Euler Equation Error

mInterestRate = zeros(nGridCapital,nGridProductivity,2);
mInterestRatio = zeros(nGridCapital,nGridProductivity,2);
expectedInterestRatio = zeros(nGridCapital,nGridProductivity,2);

for j = 1:2
    for nProductivity = 1:nGridProductivity
        for nCapital = 1:nGridCapital
            mInterestRate(nCapital,nProductivity,j) = 1 - ddelta + aalpha*vProductivity(nProductivity)*...
                (mPolicyFunction(nCapital,nProductivity,j)^(aalpha-1))*...
                (mLabor(nCapital,nProductivity,j)^(1-aalpha));
        end
    end
    
    for nProductivity = 1:nGridProductivity
        for nCapital = 1:nGridCapital
            mInterestRatio(nCapital,nProductivity,j) = ...
                mInterestRate(nCapital,nProductivity,j)/...
                mConsumptionFunction(nCapital,nProductivity,j);
        end
    end
    
    expectedInterestRatio(:,:,j) = mInterestRatio(:,:,j)*mTransition';
    
    for nProductivity = 1:nGridProductivity
        for nCapital = 1:nGridCapital
            mEulerError(nCapital,nProductivity,j) = log10(abs(1 - ...
                mConsumptionFunction(nCapital,nProductivity,j)*bbeta*...
                expectedInterestRatio(nCapital,nProductivity,j)));
        end
    end
    
end

toc

%% 7. Plotting results

figure(1)

subplot(2,2,1)
hold on
plot(vGridCapital,mValueFunction(:,3,1),'r')
plot(vGridCapital,mValueFunction(:,3,2))
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
legend('Standard VFI','Accelerated VFI')
title('Value Function for z=0')
hold off

subplot(2,2,2)
hold on
plot(vGridCapital,mPolicyFunction(:,3,1),'r')
plot(vGridCapital,mPolicyFunction(:,3,2))
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
legend('Standard VFI','Accelerated VFI')
title('Policy Function for z=0')
hold off

subplot(2,2,3)
hold on
plot(vGridCapital,mLabor(:,3,1),'r')
plot(vGridCapital,mLabor(:,3,2))
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
legend('Standard VFI','Accelerated VFI')
title('Labor Function for z=0')
hold off

subplot(2,2,4)
hold on
plot(vGridCapital,mEulerError(:,3,1),'r')
plot(vGridCapital,mEulerError(:,3,2))
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
legend('Standard VFI','Accelerated VFI')
title('Euler Error for z=0')
hold off

print -depsc2 HW2_Q3_accelerator.eps
