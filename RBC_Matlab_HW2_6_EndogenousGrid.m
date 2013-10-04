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
nGridProductivity = length(vProductivity);

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
resourcesSteadyState = consumptionSteadyState + capitalSteadyState;

% Labor Grid
vGridLabor = 0.5*laborSteadyState:0.005:1.5*laborSteadyState; %TODO: change step?
nGridLabor = length(vGridLabor);

%% 3. Calibrate Disutility of Labor
ppsi = (1/consumptionSteadyState)*3*(1-aalpha)*((((1/bbeta)+ ddelta - 1)/aalpha)^(aalpha/(aalpha-1))); %disutility of labor - added 9/12/13
valueSteadyState = (1-bbeta)*consumptionSteadyState - ppsi*((laborSteadyState^2)/2);
fprintf('Calibration\n');
fprintf(' Psi = %2.6f\n', ppsi);
fprintf('Steady State Values\n');
fprintf(' Output = %2.6f, Labor = %2.6f\n',outputSteadyState, laborSteadyState);
fprintf(' Capital = %2.6f, Consumption = %2.6f\n',  capitalSteadyState, consumptionSteadyState);
fprintf('\n');

%% 4. Required matrices and vectors

% We generate the grid of capital
vGridCapital = 0.5*capitalSteadyState:0.00006:1.5*capitalSteadyState; %TODO: change step to 0.00001
nGridCapital = length(vGridCapital);
vGridResources = zeros(nGridCapital,nGridProductivity);
for nProductivity = 1:nGridProductivity
    vGridResources(:,nProductivity) = vProductivity(nProductivity)*...
        (laborSteadyState^(1-aalpha))*(vGridCapital.^aalpha)+...
        vGridCapital.*(1-ddelta);
end

mConsumption = zeros(nGridCapital,nGridProductivity);
mResources        = zeros(nGridCapital,nGridProductivity);
mLabor            = zeros(nGridCapital,nGridProductivity);
mValueFunction = valueSteadyState*ones(nGridCapital,nGridProductivity); %#ok<NASGU>
mValueFunctionRes = zeros(nGridCapital,nGridProductivity);
mValueFunctionTilda = zeros(nGridCapital,nGridProductivity);
mValueFunctionNew = zeros(nGridCapital,nGridProductivity);
mValueFunctionTildaNew = zeros(nGridCapital,nGridProductivity);
mValueFunctionTildaDerivative = zeros(nGridCapital,nGridProductivity);
mPolicyFunction = zeros(nGridCapital,nGridProductivity);
expectedValueFunction = zeros(nGridCapital,nGridProductivity);

%% 5. We pre-build output and guess for V-tilda at labor = laborSteadyState

mOutput = (vGridCapital'.^aalpha)*vProductivity*(laborSteadyState^(1-aalpha));

% Use straight line as guess
for nCapital = 1:nGridCapital
    mValueFunctionTilda(nCapital,:) = valueSteadyState + (nCapital/nGridCapital);
end

%% 6. Endogenous Grid iteration (labor fixed at Steady State)

maxDifference = 10.0;
tolerance = 0.0000001;
iteration = 0;

fprintf('Endogenous Grid VFI\n')

while (maxDifference>tolerance)
    
    
    %Derivative V-tilda
    mValueFunctionTildaDerivative(1,:) =(mValueFunctionTilda(2,:)-...
        mValueFunctionTilda(1,:))/nGridCapital;
    mValueFunctionTildaDerivative(nGridCapital,:) =...
        (mValueFunctionTilda(nGridCapital,:)-...
        mValueFunctionTilda(nGridCapital-1,:))/nGridCapital;
    for iCapital = 2:nGridCapital-1
        mValueFunctionTildaDerivative(iCapital,:) =(mValueFunctionTilda(iCapital+1,:)-...
            mValueFunctionTilda(iCapital-1,:))/(2*nGridCapital);
    end
    
    %Solve for consumption
    mConsumption = 1./mValueFunctionTildaDerivative;
    
    %Determine endogenous resources in period t (hence k)
    mResources = bsxfun(@plus,vGridCapital',mConsumption);
    
    for nProductivity = 1:nGridProductivity
        for nCapital = 1:nGridCapital
            mValueFunctionRes(nCapital,nProductivity)=...
                log(mConsumption(nCapital,nProductivity))-...
                (ppsi/2)*(laborSteadyState^(2))+...
                mValueFunctionTilda(nCapital,nProductivity);
        end
    end
    
    mValueFunctionRes(mValueFunctionRes<-1000)=-1000;
    
    for nProductivity = 1:nGridProductivity
        mResources(:,nProductivity) = sortrows(mResources(:,nProductivity));
        F = griddedInterpolant({mResources(:,nProductivity)},...
            mValueFunctionRes(:,nProductivity));
        mValueFunctionRes(:,nProductivity)=F(vGridResources(:,nProductivity));
    end
    
    mValueFunctionTildaNew = bbeta*mValueFunctionRes*mTransition';
    
    maxDifference = max(max(abs(...
        mValueFunctionTildaNew-mValueFunctionTilda)));
    mValueFunctionTilda = mValueFunctionTildaNew;
    
    iteration = iteration+1;
    fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, ...
        maxDifference);
    
end

fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference);
fprintf('\n')

mValueFunction = mValueFunctionTilda;

%% 7. Recover Policy Functions

fprintf('Standard VFI to Recover Policy Functions')

maxDifference = 10.0;
tolerance = 0.0000001;
iteration = 0;

%Reset Output
mOutput = zeros(nGridCapital,nGridProductivity,nGridLabor);

for j = 1:nGridLabor
    
    mOutput(:,:,j) = (vGridCapital'.^aalpha)*vProductivity*(vGridLabor(j)^...
        (1-aalpha));
    
end

while (maxDifference>tolerance)
    
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

figure(1)

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
plot(vGridCapital,mLabor)
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
ylabel('log10|EulerError|')
title('Euler Error')

print -depsc2 HW2_Q6_endogenousgrid.eps

