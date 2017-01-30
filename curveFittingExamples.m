% This script demonstrates how to use the curve fitting functions to fit
% morphogen gradient data. It starts by generating noisy concentration
% profiles, and then it fits to these profiles.

% Clear the workspace:
clear all;

% Define the parameter values:
% The amplitude of the concentration profiles:
amplitude = 1;
% The decay length of the exponential function:
lambdaExponential = 0.5;
% The proximal intrinsic decay length of the two-domain and
% two-domain-gradual-sink functions
lambdaProximal = 0.5;
% The distal intrinsic decay length of the two-domain model:
lambdaDistal = 0.1;
% The slope of increase of the consumption rate in the distal domain for the
% two-domain-gradual-sink model:
consumptionRateSlopeDistal = 100;
% The location of the interface boundary for the two-domain and
% two-domain-gradual-sink models:
interfaceBoundary = 1;
% The value of the offset:
offsetVal = 0;

% The magnitude (standard deviation) of the multiplicative and additive
% noises
multNoiseStd = 0.05;
addNoiseStd = 0.05;

% The grid:
x = 0:0.02:3;

% For each of the functions, store the parameter values in an array,
% ensuring to follow the correct order:
parameters.Exponential = [amplitude lambdaExponential offsetVal];
parameters.TwoDomain = [amplitude lambdaProximal lambdaDistal offsetVal];
parameters.TwoDomainGradualSink = [amplitude lambdaProximal consumptionRateSlopeDistal offsetVal];

% Define the landmarks
landmarks.zeroLocation = 0;
landmarks.interfaceBoundaryLocation = interfaceBoundary;

% Specify how to deal with the offset:
offset.value = offsetVal;
offset.mode = 'free';

% Generate the noise-free data for each of the functions:
yNoiseFree.linear = amplitude*(1-x./x(end));
yNoiseFree.exponential = decayingExponential(parameters.exponential,x,landmarks,offset);
yNoiseFree.twoDomain = twoDomain(parameters.twoDomain,x,landmarks,offset);
yNoiseFree.twoDomainGradualSink = twoDomainGradualSink(parameters.twoDomainGradualSink,x,landmarks,offset);

% Generate the additive noise:
addNoise = addNoiseStd*randn(1,length(x));

% Generate the multiplicative noise for each of the functions, and then add
% both the additive and multiplicative noises to the noise-free data:
yNoiseFreeFields = fields(yNoiseFree);
for i = 1:length(yNoiseFreeFields)
    multNoise = multNoiseStd*randn(1,length(x)).*yNoiseFree.(yNoiseFreeFields{i});
    yNoisy.(yNoiseFreeFields{i}) = yNoiseFree.(yNoiseFreeFields{i}) + multNoise + addNoise;
end

% Fit the exponential function to the noisy exponential data, considering
% the offset as a free fitting parameter:
% Offset is free parameter:
offset.mode = 'free';
% Fit only the exponential function:
fitFlags.fitTwoDomainModel = false;
fitFlags.fitTwoDomainGradualSinkModel = false;
% Perform the fit:
fitResults = morphogenGradientCurveFitting(x,yNoisy.exponential,offset,landmarks,fitFlags);
% The curve corresponding to the best-fit parameters:
yBestFit = decayingExponential(fitResults.exponential.P,x,landmarks,offset);
% Plot the fit, and indicate the estimate and confidence interval of the
% decay length:
figure;
plot(x,yNoisy.exponential,x,yBestFit);
text(1,.8,...
    sprintf('best-fit decay length: %.2f\n95%% CI: [%.2f,%.2f]',...
    fitResults.exponential.P(2),...
    fitResults.exponential.ci(2,1),...
    fitResults.exponential.ci(2,2)...
    ));
title('Fitting decaying exponential with offset term as free parameter');

% Fit the two-domain solution to the noisy two-domain data, considering
% the offset as a free fitting parameter:
% Offset is free parameter:
offset.mode = 'free';
% Fit two-domain model:
fitFlags.fitTwoDomainModel = true;
fitFlags.fitTwoDomainGradualSinkModel = false;
% Perform the fit:
fitResults = morphogenGradientCurveFitting(x,yNoisy.twoDomain,offset,landmarks,fitFlags);
% The curve corresponding to the best-fit parameters:
yBestFit = twoDomain(fitResults.twoDomain.P,x,landmarks,offset);
% Plot the fit, and indicate the estimate and confidence interval of the
% decay lengths:
figure;
plot(x,yNoisy.twoDomain,x,yBestFit);
text(1,.8,...
    sprintf('best-fit proximal decay length: %.2f\n95%% CI: [%.2f,%.2f]',...
    fitResults.twoDomain.P(2),...
    fitResults.twoDomain.ci(2,1),...
    fitResults.twoDomain.ci(2,2)...
    ));
text(1,.6,...
    sprintf('best-fit distal decay length: %.2f\n95%% CI: [%.2f,%.2f]',...
    fitResults.twoDomain.P(3),...
    fitResults.twoDomain.ci(3,1),...
    fitResults.twoDomain.ci(3,2)...
    ));
title('Fitting two-domain model with offset term as free parameter');

% Fit the two-domain-gradual-sink solution to the noisy 
% two-domain-gradual-sink data, considering the offset as a free fitting 
% parameter:
% Offset is free parameter:
offset.mode = 'free';
% Fit two-domain-gradual-sink model:
fitFlags.fitTwoDomainModel = false;
fitFlags.fitTwoDomainGradualSinkModel = true;
% Perform the fit:
fitResults = morphogenGradientCurveFitting(x,yNoisy.twoDomainGradualSink,offset,landmarks,fitFlags);
% The curve corresponding to the best-fit parameters:
yBestFit = twoDomainGradualSink(fitResults.twoDomainGradualSink.P,x,landmarks,offset);
% Plot the fit, and indicate the estimate and confidence interval of the
% proximal decay length and the slope of increase for the distal consumption
% rate constant:
figure;
plot(x,yNoisy.twoDomainGradualSink,x,yBestFit);
text(1,.8,...
    sprintf('best-fit proximal decay length: %.2f\n95%% CI: [%.2f,%.2f]',...
    fitResults.twoDomainGradualSink.P(2),...
    fitResults.twoDomainGradualSink.ci(2,1),...
    fitResults.twoDomainGradualSink.ci(2,2)...
    ));
text(1,.6,...
    sprintf('best-fit distal consumption rate slope: %.2f\n95%% CI: [%.2f,%.2f]',...
    fitResults.twoDomainGradualSink.P(3),...
    fitResults.twoDomainGradualSink.ci(3,1),...
    fitResults.twoDomainGradualSink.ci(3,2)...
    ));
title('Fitting two-domain-gradual-sink model with offset term as free parameter');

% Fit the exponential function to the noisy linear data, considering the 
% offset as a free fitting parameter. We will see that this results in very
% unreasonable estimates for the parameter values.
% Offset is free parameter:
offset.mode = 'free';
% Fit only the exponential function:
fitFlags.fitTwoDomainModel = false;
fitFlags.fitTwoDomainGradualSinkModel = false;
% Perform the fit:
fitResults = morphogenGradientCurveFitting(x,yNoisy.linear,offset,landmarks,fitFlags);
% The curve corresponding to the best-fit parameters:
yBestFit = decayingExponential(fitResults.exponential.P,x,landmarks,offset);
% Plot the fit, and indicate the estimate and confidence interval of the
% decay length:
figure;
plot(x,yNoisy.linear,x,yBestFit);
text(1,.8,...
    sprintf('best-fit decay length: %.2f\n95%% CI: [%.2f,%.2f]',...
    fitResults.exponential.P(2),...
    fitResults.exponential.ci(2,1),...
    fitResults.exponential.ci(2,2)...
    ));
title('Fitting exponential function to linear data, with offset term as a free parameter, produces unreliable estimate');

% Fit the exponential function to the noisy linear data, considering the 
% offset as a fixed parameter. We will see that this results in more
% reasonable estimates for the parameter values, though the fit still isn't
% great.
% Offset is fixed parameter:
offset.mode = 'fixed';
% Fit only the exponential function:
fitFlags.fitTwoDomainModel = false;
fitFlags.fitTwoDomainGradualSinkModel = false;
% Perform the fit:
fitResults = morphogenGradientCurveFitting(x,yNoisy.linear,offset,landmarks,fitFlags);
% The curve corresponding to the best-fit parameters:
yBestFit = decayingExponential(fitResults.exponential.P,x,landmarks,offset);
% Plot the fit, and indicate the estimate and confidence interval of the
% decay length:
figure;
plot(x,yNoisy.linear,x,yBestFit);
text(1,.8,...
    sprintf('best-fit decay length: %.2f\n95%% CI: [%.2f,%.2f]',...
    fitResults.exponential.P(2),...
    fitResults.exponential.ci(2,1),...
    fitResults.exponential.ci(2,2)...
    ));
title('Fitting exponential function to linear data, with offset term as a fixed parameter, produces more reasonable, albeit not perfect, parameter estimates');