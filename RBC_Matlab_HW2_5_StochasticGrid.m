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

%% 2. Steady State

laborSteadyState = 1/3; %added 9/10/13
capitalSteadyState = laborSteadyState*((((1/bbeta)+ddelta-1)/aalpha)^(1/(aalpha-1))); %modified 9/10/13
outputSteadyState = (capitalSteadyState^aalpha)*(laborSteadyState^(1-aalpha)); %modified 9/10/13
consumptionSteadyState = outputSteadyState-ddelta*capitalSteadyState; %modified 9/10/13

%% 3. Calibrate Disutility of Labor
ppsi = (1/consumptionSteadyState)*3*(1-aalpha)*((((1/bbeta)+ ddelta - 1)/aalpha)^(aalpha/(aalpha-1))); %disutility of labor - added 9/12/13
fprintf('Calibration\n');
fprintf(' Psi = %2.6f\n', ppsi);
fprintf('Steady State Values\n');
fprintf(' Output = %2.6f, Labor = %2.6f\n',outputSteadyState, laborSteadyState);
fprintf(' Capital = %2.6f, Consumption = %2.6f\n',  capitalSteadyState, consumptionSteadyState);
fprintf('\n');

%% 4. Grid Loop to Change Grid

%Set up cell array of capital grids and labor grid
nCapitalGrids = 2;

%Deterministic Uniform Grid
aCapitalGrids{1} = 0.5*capitalSteadyState:0.00006:1.5*capitalSteadyState; %FIXME to 0.00006
%Stochastic Truncated Normal Grid
norm = makedist('Normal','mu',capitalSteadyState);
pd = truncate(norm,0.5*capitalSteadyState,1.5*capitalSteadyState);
vStochasticGrid = sort(random(pd,1,length(aCapitalGrids{1})),2);
vStochasticGrid(1,1) = 0.5*capitalSteadyState; %include lower bound
vStochasticGrid(1,length(aCapitalGrids{1})) = 1.5*capitalSteadyState; %include upper bound
aCapitalGrids{2} = vStochasticGrid;

nGridProductivity = length(vProductivity);
vGridLabor = 0.5*laborSteadyState:0.005:1.5*laborSteadyState; %TODO: change step?
nGridLabor = length(vGridLabor);
vTime = zeros(nCapitalGrids,1);

for iCapitalStep = 1:nCapitalGrids
    
    timeVFI = tic;
    vGridCapital = aCapitalGrids{iCapitalStep};
    nGridCapital = length(vGridCapital);
    
    fprintf('Grid Size = %d', nGridCapital);
    fprintf(' points\n');
    
    %% 5. Required matrices and vectors
    mOutput           = zeros(nGridCapital,nGridProductivity,nGridLabor);
    mLabor            = zeros(nGridCapital,nGridProductivity);
    mValueFunction    = zeros(nGridCapital,nGridProductivity);
    mValueFunctionNew = zeros(nGridCapital,nGridProductivity);
    mPolicyFunction   = zeros(nGridCapital,nGridProductivity);
    mConsumptionFunction = zeros(nGridCapital,nGridProductivity);
    expectedValueFunction = zeros(nGridCapital,nGridProductivity);
    
    %% 6. We pre-build output for each point in the grid
    
    for j = 1:nGridLabor
        
        mOutput(:,:,j) = (vGridCapital'.^aalpha)*vProductivity*(vGridLabor(j)^...
            (1-aalpha));
        
    end
    
    %% 7. Main iteration
    
    maxDifference = 10.0;
    maxIteration = 500; %TODO: change?
    tolerance = 10^(-6); %TODO: change?
    iteration = 0;
    
    while (maxDifference>tolerance)&&(iteration<maxIteration)
        
        expectedValueFunction = mValueFunction*mTransition';
        
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
                            (vGridLabor(nLabor)^2)/2)+bbeta*expectedValueFunction(nCapitalNextPeriod,nProductivity);
                        
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
                
                mValueFunctionNew(nCapital,nProductivity) = valueHighSoFar1;
                mPolicyFunction(nCapital,nProductivity) = capitalChoice;
                mLabor(nCapital,nProductivity) = laborChoice; %Labor function choosing next period capital optimally
                mConsumptionFunction(nCapital,nProductivity) = consumption;
                
            end %end for capital grid loop
            
        end %end of productivity grid loop
        
        maxDifference = max(max(abs(mValueFunctionNew-mValueFunction)));
        mValueFunction = mValueFunctionNew;
        
        iteration = iteration+1;
        if (mod(iteration,10)==0 || iteration ==1)
            fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference);
        end
        
    end
    
    fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference);
    fprintf('\n')
    
    vTime(iCapitalStep,1) = toc(timeVFI);
    fprintf('Time for VFI = %2.8f', vTime(iCapitalStep,1));
    fprintf(' seconds\n');
    fprintf('\n')
    
    %% 8. Euler Error Analysis
    
    mInterestRate = zeros(nGridCapital,nGridProductivity);
    mInterestRatio = zeros(nGridCapital,nGridProductivity);
    expectedInterestRatio = zeros(nGridCapital,nGridProductivity);
    mEulerError = zeros(nGridCapital,nGridProductivity);
    
    for nProductivity = 1:nGridProductivity
        for nCapital = 1:nGridCapital
            mInterestRate(nCapital,nProductivity) = 1 - ddelta + aalpha*vProductivity(nProductivity)*...
                (mPolicyFunction(nCapital,nProductivity)^(aalpha-1))*...
                (mLabor(nCapital,nProductivity)^(1-aalpha));
        end
    end
    
    for nProductivity = 1:nGridProductivity
        for nCapital = 1:nGridCapital
            mInterestRatio(nCapital,nProductivity) = ...
                mInterestRate(nCapital,nProductivity)/...
                mConsumptionFunction(nCapital,nProductivity);
        end
    end
    
    expectedInterestRatio(:,:) = mInterestRatio(:,:)*mTransition';
    
    for nProductivity = 1:nGridProductivity
        for nCapital = 1:nGridCapital
            mEulerError(nCapital,nProductivity) = log10(abs(1 - ...
                mConsumptionFunction(nCapital,nProductivity)*bbeta*...
                expectedInterestRatio(nCapital,nProductivity)));
        end
    end
    
    
    
    %% 9. Plotting results
    
    %{
    figure;
    
    subplot(2,2,1)
    plot(vGridCapital,mValueFunction)
    xlim([vGridCapital(1) vGridCapital(nGridCapital)])
    title('Value Function')
    
    subplot(2,2,2)
    plot(vGridCapital,mPolicyFunction)
    xlim([vGridCapital(1) vGridCapital(nGridCapital)])
    title('Policy Function')
    
    subplot(2,2,3)
    plot(vGridCapital,mLabor)
    xlim([vGridCapital(1) vGridCapital(nGridCapital)])
    title('Labor Function')
    
    subplot(2,2,4)
    plot(vGridCapital,mEulerError)
    xlim([vGridCapital(1) vGridCapital(nGridCapital)])
    ylabel('
    title('Euler Equation Error')
    
    %}
    
    %% 10. Store for Comparison when z = 0
    
    aComparison{iCapitalStep} = [vGridCapital' mValueFunction(:,3)...
        mPolicyFunction(:,3) mLabor(:,3) mEulerError(:,3)]; %#ok<SAGROW>
    
end

%% 11. Plot Comparison between Grids when z = 0

%Compare Grids
gridStatistic = zeros(2,5);
for j = 1:2
    gridStatistic(j,1) = mean(aCapitalGrids{j});
    gridStatistic(j,2) = std(aCapitalGrids{j});
    gridStatistic(j,3) = skewness(aCapitalGrids{j});
    gridStatistic(j,4) = kurtosis(aCapitalGrids{j});
    gridStatistic(j,5) = length(aCapitalGrids{j});
end
fprintf('\nGrid Statistics\n')
fprintf(' %1.5f\n',gridStatistic);
fprintf('Grid1 is determinstic, Grid2 is stochastic\n')

figure;

subplot(1,2,1)
hist(aCapitalGrids{1},20);
xlim([0.5*capitalSteadyState 1.5*capitalSteadyState])
title('Histogram for Deterministic Grid')

subplot(1,2,2)
hist(aCapitalGrids{2},20);
xlim([0.5*capitalSteadyState 1.5*capitalSteadyState])
title('Histogram for Stochastic Grid')

figure;

subplot(2,2,1)
%Value Function
hold on
plot(aComparison{1}(:,1),aComparison{1}(:,2))
plot(aComparison{2}(:,1),aComparison{2}(:,2),'r')
xlim([0.5*capitalSteadyState 1.5*capitalSteadyState])
legend('Deterministic', 'Stochastic')
title('Value Function for z=0')
hold off

subplot(2,2,2)
%Policy Function
hold on
plot(aComparison{1}(:,1),aComparison{1}(:,3))
plot(aComparison{2}(:,1),aComparison{2}(:,3),'r')
xlim([0.5*capitalSteadyState 1.5*capitalSteadyState])
legend('Deterministic', 'Stochastic')
title('Policy Function for z=0')
hold off

subplot(2,2,3)
%Policy Function
hold on
plot(aComparison{1}(:,1),aComparison{1}(:,4))
plot(aComparison{2}(:,1),aComparison{2}(:,4),'r')
legend('Deterministic', 'Stochastic')
xlim([0.5*capitalSteadyState 1.5*capitalSteadyState])
title('Labor Function for z=0')
hold off

subplot(2,2,4)
%Policy Function
hold on
plot(aComparison{1}(:,1),aComparison{1}(:,5))
plot(aComparison{2}(:,1),aComparison{2}(:,5),'r')
legend('Deterministic', 'Stochastic')
xlim([0.5*capitalSteadyState 1.5*capitalSteadyState])
title('Euler Error for z=0')
ylabel('Log10|Euler Equation Error|')
hold off

toc