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
fprintf('Calibration\n');
fprintf(' Psi = %2.6f\n', ppsi);
fprintf('Steady State Values\n');
fprintf(' Output = %2.6f, Labor = %2.6f\n',outputSteadyState, laborSteadyState);
fprintf(' Capital = %2.6f, Consumption = %2.6f\n',  capitalSteadyState, consumptionSteadyState);
fprintf('\n');

%% 3. Initial Guess (interpolated multigrid)

fprintf('Guess V-Tilda using interpolated multigrid\n')

%Multigrid first Step
nCapitalGrids = 1;
for iCapitalStep = 1:nCapitalGrids
    capitalStepSize = 6*(10^-(iCapitalStep+2)); %GridSizes:199,1985,19849
    aCapitalGrids{iCapitalStep} = 0.5*capitalSteadyState:capitalStepSize:1.5*capitalSteadyState; %#ok<SAGROW>
end
vTime = zeros(nCapitalGrids+1,1);

iTime = tic;
for iCapitalStep = 1:1
    
    timeVFI = tic;
    
    %Set old grid and value function if not first loop
    if (iCapitalStep>1)
        vOldGridCapital = aCapitalGrids{iCapitalStep-1};
        mOldValueFunction = mValueFunction;
    end
    
    vGridCapital = aCapitalGrids{iCapitalStep};
    nGridCapital = length(vGridCapital);
    
    fprintf('Grid Step %d', iCapitalStep);
    fprintf('\n');
    fprintf('Grid Size = %d', nGridCapital);
    fprintf(' points\n');
    
    %% 6. Required matrices and vectors
    mOutput           = zeros(nGridCapital,nGridProductivity,nGridLabor);
    mLabor            = zeros(nGridCapital,nGridProductivity);
    mValueFunction    = zeros(nGridCapital,nGridProductivity);
    mValueFunctionNew = zeros(nGridCapital,nGridProductivity);
    mPolicyFunction   = zeros(nGridCapital,nGridProductivity);
    mConsumptionFunction = zeros(nGridCapital,nGridProductivity);
    expectedValueFunction = zeros(nGridCapital,nGridProductivity); %#ok<NASGU>
    
    %% 7. Interpolate new guess for value function
    if (iCapitalStep>1)
        F = griddedInterpolant({vOldGridCapital',vProductivity},mOldValueFunction,'spline');
        mValueFunction = F({vGridCapital',vProductivity});
    end
    
    %% 8. We pre-build output for each point in the grid
    
    for j = 1:nGridLabor
        
        mOutput(:,:,j) = (vGridCapital'.^aalpha)*vProductivity*(vGridLabor(j)^...
            (1-aalpha));
        
    end
    
    %% 9. Main iteration
    
    maxDifference = 10.0;
    maxIteration = 500; %TODO: change?
    tolerance = 10^-(iCapitalStep+3); %FIXME: tolerance level increases with grid increase
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
    fprintf('Time for step VFI = %2.8f', vTime(iCapitalStep,1));
    fprintf(' seconds\n');
    fprintf('\n')
    
end
fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference);
fprintf('\n')


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
mValueFunctionRes = zeros(nGridCapital,nGridProductivity);
mValueFunctionTildaNew = zeros(nGridCapital,nGridProductivity);
mValueFunctionTildaDerivative = zeros(nGridCapital,nGridProductivity);
mPolicyFunction = zeros(nGridCapital,nGridProductivity);
expectedValueFunction = zeros(nGridCapital,nGridProductivity);

%% 5. We pre-build output and guess for V-tilda at labor = laborSteadyState

mOutput = (vGridCapital'.^aalpha)*vProductivity*(laborSteadyState^(1-aalpha));

%Interpolation
vOldGridCapital = aCapitalGrids{nCapitalGrids};
nOldGripCapital = length(vOldGridCapital);
mOldValueFunction = mValueFunction;

F = griddedInterpolant({vOldGridCapital',vProductivity},mOldValueFunction,...
    'spline');
mValueFunctionTilda = F({vGridCapital',vProductivity});



%% 6. Endogenous Grid iteration (labor fixed at Steady State)

maxDifference = 10.0;
tolerance = 0.0000001;
iteration = 0;

fprintf('Endogenous Grid VFI\n')

while (maxDifference>tolerance)
    
    %Derivative V-tilda    
    for nProductivity = 1:nGridProductivity
        tempStore = diff(mValueFunctionTilda(:,nProductivity))/0.00006;
        a = tempStore(1,1);
        mValueFunctionTildaDerivative(:,nProductivity) =...
            vertcat(a,tempStore);
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
    
    for nProductivity = 1:nGridProductivity
        mValueFunctionRes(:,nProductivity) =...
            interp1(vGridResources(:,nProductivity),...
            mValueFunctionRes(:,nProductivity),...
            'linear','extrap');
    end
    
    mValueFunctionTildaNew = bbeta*mValueFunctionRes*mTransition';
    
    maxDifference = max(max(abs(...
        mValueFunctionTildaNew-mValueFunctionTilda)));
    mValueFunctionTilda = mValueFunctionTildaNew;
    
    iteration = iteration+1;
    if (mod(iteration,10)==0 || iteration ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, ...
            maxDifference);
    end
    
end

fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference);
fprintf('\n')

toc

%{
%% 6. Plotting results

figure(1)

subplot(3,1,1)
plot(vGridCapital,mValueFunction)
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