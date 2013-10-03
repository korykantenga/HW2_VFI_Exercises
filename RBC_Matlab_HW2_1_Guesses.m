%% Basic RBC model with partial depreciation and endogenous labor

% Original by Jesus Fernandez-Villaverde at Haverford, July 31, 2013
% Modifiied by Kory Kantenga at University of Pennsylvania, Oct 2, 2013

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
capitalSteadyState = laborSteadyState*((((1/bbeta)+ddelta-1)/aalpha)^(1/(aalpha-1)));
outputSteadyState = (capitalSteadyState^aalpha)*(laborSteadyState^(1-aalpha));
consumptionSteadyState = outputSteadyState-ddelta*capitalSteadyState;

%% 3. Calibrate Disutility of Labor
ppsi = (1/consumptionSteadyState)*3*(1-aalpha)*((((1/bbeta)+ ddelta - 1)/aalpha)^(aalpha/(aalpha-1))); %disutility of labor - added 9/12/13
valueSteadyState = (1-bbeta)*consumptionSteadyState - ppsi*((laborSteadyState^2)/2);
fprintf('Calibration\n');
fprintf(' Psi = %2.6f\n', ppsi);
fprintf('Steady State Values\n');
fprintf(' Output = %2.6f, Labor = %2.6f\n',outputSteadyState, laborSteadyState);
fprintf(' Capital = %2.6f, Consumption = %2.6f\n',  capitalSteadyState, consumptionSteadyState);
fprintf('\n');

% We generate the grid of capital
vGridCapital = 0.5*capitalSteadyState:0.00001:1.5*capitalSteadyState; %TODO: change step to 0.00001 - modified 9/10/13
vGridLabor = 0.5*laborSteadyState:0.005:1.5*laborSteadyState; %TODO: change step?

nGridLabor = length(vGridLabor);
nGridCapital = length(vGridCapital);
nGridProductivity = length(vProductivity);
nValueGuess = 3; %number of initial guesses

%% 4. Required matrices and vectors

vTime             = zeros(nValueGuess,1); %store time of VFI loop
mOutput           = zeros(nGridCapital,nGridProductivity,nGridLabor);
mLabor            = zeros(nGridCapital,nGridProductivity);
mValueFunctionNew = zeros(nGridCapital,nGridProductivity);
mPolicyFunction   = zeros(nGridCapital,nGridProductivity);
expectedValueFunction = zeros(nGridCapital,nGridProductivity);
aValueFunction = zeros(nGridCapital,nGridProductivity,nValueGuess);
aPolicyFunction = zeros(nGridCapital,nGridProductivity,nValueGuess);
aLaborFunction = zeros(nGridCapital,nGridProductivity,nValueGuess);
aConsumptionFunction = zeros(nGridCapital,nGridProductivity,nValueGuess);
aEulerError       = zeros(nGridCapital,nGridProductivity,nValueGuess);

%% 5. We pre-build output for each point in the grid

for j = 1:nGridLabor
    
    mOutput(:,:,j) = (vGridCapital'.^aalpha)*vProductivity*(vGridLabor(j)^(1-aalpha));
    
end


%% 6. Select an Initial Guess

for nGuess = 1:nValueGuess
    
    if (nGuess==1)
        mValueFunction = zeros(nGridCapital,nGridProductivity);
        fprintf('Initial Guess: Constant Function equal to zero\n')
    end
    
    if (nGuess==2)
        mValueFunction = valueSteadyState*ones(nGridCapital,nGridProductivity);
        fprintf('Initial Guess: Constant Function equal to Steady State Value\n')
    end
    
    if (nGuess==3)
        mValueFunction = zeros(nGridCapital,nGridProductivity);
        for z = 1:nGridProductivity
            mValueFunction(:,z) = aValueFunction(:,5,1); %VF for z=1.0212
        end
        fprintf('Initial Guess: Deterministic Solution (z=1)\n')
    end
    
    %% 7. Main iteration
    
    maxDifference = 10.0;
    maxIteration = 1000;
    tolerance = 0.0000001;
    iteration = 0;
    
    timeVFI = tic;
    
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
                
                mValueFunctionNew(nCapital,nProductivity) = valueHighSoFar1;
                mPolicyFunction(nCapital,nProductivity) = capitalChoice;
                mLabor(nCapital,nProductivity) = laborChoice; %Labor function choosing next period capital optimally
                aConsumptionFunction(nCapital,nProductivity,nGuess) = consumption;
                
            end %end for capital grid loop
            
        end %end of productivity grid loop
        
        maxDifference = max(max(abs(mValueFunctionNew-mValueFunction)));
        mValueFunction = mValueFunctionNew;
        aValueFunction(:,:,nGuess) = mValueFunctionNew;
        aPolicyFunction(:,:,nGuess) = mPolicyFunction;
        aLaborFunction(:,:,nGuess) = mLabor;
        
        iteration = iteration+1;
        if (mod(iteration,10)==0 || iteration ==1)
            fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference);
        end
        
    end
    
    fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference);
    fprintf('\n')
    vTime(nGuess,1) = toc(timeVFI);
    fprintf('Time for VFI = %2.8f', vTime(nGuess,1));
    fprintf(' seconds\n');
    fprintf('\n')
    
    %% 8. Plotting results
    
    %{
    figure(nGuess)
    
    subplot(3,1,1)
    plot(vGridCapital,aValueFunction(:,:,nGuess))
    xlim([vGridCapital(1) vGridCapital(nGridCapital)])
    title('Value Function')
    
    subplot(3,1,2)
    plot(vGridCapital,mPolicyFunction)
    xlim([vGridCapital(1) vGridCapital(nGridCapital)])
    title('Policy Function')
    
    subplot(3,1,3)
    plot(vGridCapital,mLabor)
    xlim([vGridCapital(1) vGridCapital(nGridCapital)])
    title('Labor Function choosing Future Capital Optimally')
    %}
    
end

%% 9. Euler Equation Error Analysis

aInterestRate = zeros(nGridCapital,nGridProductivity,nValueGuess);
aInterestRatio = zeros(nGridCapital,nGridProductivity,nValueGuess);
expectedInterestRatio = zeros(nGridCapital,nGridProductivity,nValueGuess);

for nGuess = 1:nValueGuess
    for nProductivity = 1:nGridProductivity
        for nCapital = 1:nGridCapital
            aInterestRate(nCapital,nProductivity,nGuess) = 1 - ddelta + aalpha*vProductivity(nProductivity)*...
                (aPolicyFunction(nCapital,nProductivity,nGuess)^(aalpha-1))*...
                (aLaborFunction(nCapital,nProductivity,nGuess)^(1-aalpha));
        end
    end
    
    for nProductivity = 1:nGridProductivity
        for nCapital = 1:nGridCapital
            aInterestRatio(nCapital,nProductivity,nGuess) = ...
                aInterestRate(nCapital,nProductivity,nGuess)/...
                aConsumptionFunction(nCapital,nProductivity,nGuess);
        end
    end
    
    expectedInterestRatio(:,:,nGuess) = aInterestRatio(:,:,nGuess)*mTransition';
    
    for nProductivity = 1:nGridProductivity
        for nCapital = 1:nGridCapital
            aEulerError(nCapital,nProductivity,nGuess) = log10(abs(1 - ...
                aConsumptionFunction(nCapital,nProductivity,nGuess)*bbeta*...
                expectedInterestRatio(nCapital,nProductivity,nGuess)));
        end
    end
    
end

%% 10. Plot Comparisons

figure;

subplot(2,2,1)
hold on
plot(vGridCapital,aValueFunction(:,:,2)-aValueFunction(:,:,1),'-x')
plot(vGridCapital,aValueFunction(:,:,3)-aValueFunction(:,:,2),'--')
legend('z1','z2','z3','z4','z5')
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
ylim([-.000001 .000001])
title('Deviations from Value Function with Guess Zero')
xlabel('Capital')
hold off

subplot(2,2,2)
hold on
plot(vGridCapital,aPolicyFunction(:,:,2)-aPolicyFunction(:,:,1),'-x')
plot(vGridCapital,aPolicyFunction(:,:,3)-aPolicyFunction(:,:,2),'--')
legend('z1','z2','z3','z4','z5')
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
ylim([-.000001 .000001])
title('Deviations from Policy Function with Guess Zero')
xlabel('Capital')
hold off

subplot(2,2,3)
hold on
plot(vGridCapital,aLaborFunction(:,:,2)-aLaborFunction(:,:,1),'-x')
plot(vGridCapital,aLaborFunction(:,:,3)-aLaborFunction(:,:,2),'--')
legend('z1','z2','z3','z4','z5')
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
ylim([-.000001 .000001])
title('Deviations from Labor Function with Guess Zero')
xlabel('Capital')
hold off

subplot(2,2,4)
hold on
plot(vGridCapital,aEulerError(:,3,1),'-x')
plot(vGridCapital,aEulerError(:,3,2),'r--d')
plot(vGridCapital,aEulerError(:,3,3),'g')
legend('Guess Zero','Steady State Guess','Deterministic Guess')
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
title('Euler Error for z = 0')
ylabel('Log10|Euler Equation Error|')
xlabel('Capital')
hold off

toc





