function fitResults = morphogenGradientCurveFitting(x,y,offset,landmarks,flags)
%MORPHOGENGRADIENTCURVEFITTING  Fit solutions of uniform consumption and
%two-domain models to morphogen gradient data.
% fitResults = MORPHOGENGRADIENTCURVEFITTING(x,y,offset,landmarks,flags)
%
% This function was developed to fit the solutions of various models to
% morphogen gradient data. In particular, it was developed for our paper in
% which we argue that the Dpp morphogen gradient in the fruit fly wing
% forms according to a two-domain model.
%
% This function is capable of fitting the steady-state solutions to three 
% models, and the user can specify which solution(s) to fit. The three 
% options are:
%   1. Uniform consumption model (i.e., simple decaying exponential)
%   2. Two-domain model
%   3. Two-domain-gradual-sink model
% These three models are described in the Supplemental Modeling Notes of 
% the paper.
%
% The user has the option to specify: 
%	1. Which models to fit (specified in 'flags').
%	2. The location of the origin (x = 0) from which fitting should start
%	(the data in x will be zeroed according to the value in
%	'landmarks.zeroLocation').
%   3. Whether the offset term should be kept as a free fitting parameter
%   or constrained.
%
%
% *** INPUT ARGUMENTS ***
%
% 'x': a vector specifying the values on the abscissa
% 'y': a vector specifying the values on the ordinate
% 'offset': a structure specifying how to account for the background
% fluorescence. It contains two fields:
%	1. '.mode': set as 'fixed' offset term is to be constrained; 
%   set as 'free' of the offset term is to be kept as a free fitting
%   parameter
%   2. '.value': specifies the offset if .mode = 'fixed';
%   specifies the initial guess for the offset if '.mode' = 'free'. If
%   '.mode' = 'free' and '.value' is not a field, the program will take the 
%   initial guess to be the minimum of 'y'.
% 'landmarks': a structure specifying the locations of landmarks in the
% tissue. It contains two fields:
%   1. '.zeroLocation': location of the origin (along the vector 'x'). If 
%   '.zeroLocation' is not a field, the origin is taken to occur at 
%   x = 0.
%   2. '.interfaceBoundaryLocation': location of the interface boundary
%   (along the vector 'x'). '.interfaceBoundaryLocation' does not need to
%   be specified if the user is only fitting the decaying exponential
%   function (since that function doesn't depend on an interface boundary).
% 'flags': a structure specifying which models to fit. The decaying
% exponential will always be fitted (since the other fits depend on it to
% generate the initial guess for the parameters). Thus, this structure has
% two fields:
%   1. '.fitTwoDomainModel': boolean variable specifying whether to fit the
%   two-domain model.
%   2. '.fitTwoDomainGradualSinkModel': boolean variable specifying whether
%   to fit the two-domain-gradual-sink model.
%
%
% *** OUTPUT ARGUMENTS ***
% 'fitResults': a structure containing a field for each of the models the
% user opted to fit:
%   1. '.exponential': A structure specifying the results of the decaying
%   exponential fit.
%   2. '.twoDomain': A structure specifying the results of the two-domain
%   fit.
%   3. '.twoDomainGradualSink': A structure specifying the results of the
%   two-domain-gradual-sink fit.
% Each of these fields are also, in turn, a structure which specifies the
% results of the fit, containing the following fields:
%   1. '.P': A vector specifying the best-fit parameter values, ordered as
%   follows:
%       '.exponential.P(1)': amplitude
%       '.exponential.P(2)': decay length
%       '.exponential.P(3)' (iff 'offset.mode' = 'free'): offset
%       '.twoDomain.P(1)': amplitude
%       '.twoDomain.P(2)': proximal intrinsic decay length
%       '.twoDomain.P(3)': distal intrinsic decay length
%       '.twoDomain.P(4)' (iff 'offset.mode' = 'free'): offset
%       '.twoDomainGradualSink.P(1)': amplitude
%       '.twoDomainGradualSink.P(2)': proximal intrinsic decay length
%       '.twoDomainGradualSink.P(3)': distal gradual sink slope
%       '.twoDomainGradualSink.P(4)' (iff 'offset.mode' = 'free'): offset
%   2. '.ci': for n fit parameters, an nx2 matrix specifying the 95%
%   confidence intervals
%   2. '.mse': Mean-squared error
%   3. '.R2': Coefficient of determination (R^2)
%   4. '.warningReturned': Whether Matlab returned a warning for this fit
%   result
%
% 
% *** EXAMPLES ***
%
% Refer to 'morphogenGradientFittingExamples.m'.
%
%
% *** IMPLEMENTATION ***
%
% Assuming that the data contains background, but is otherwise well-fit by
% an exponential, then it can be modeled as 'y = a*exp(-x/lambda) + b',
% where 'a' is the amplitude, 'lambda' is the decay length of the
% exponential, and 'b' is the offset term accounting for the background
% (which can be either fixed at a certain value, or kept free as another
% fitting parameter). To obtain an initial guess in order to start the
% non-linear fitting, let's assume we can get a guess for the value of the
% offset term 'b', and then subtract that value from 'y' to get 'y2'; 
% 'y2 = y - b = a*exp(-x/lambda)'. Then, log-transform this equation:
% 'log(y2) = log(a) - x/lambda'. This equation can be formed as a linear
% fitting problem in order to get estimates for a and lambda. 

% We then use these estimates as the initial guess for the non-linear
% problem of fitting the exponential function to the original data (which
% hasn't been log-transformed). 
%
% If the user specifies to also fit the two-domain model and/or the
% two-domain-gradual-sink model, then we also fit these models using 
% non-linear fitting, using, as the initial guesses for the parameter
% values, the values of the best-fit parameters for the decaying 
% exponential function (obtained as described above above).
% 
%
% ******
% Created by Abed Alnaif, abed.alnaif@gmail.com
% Tested in Matlab R2012b
% ******


% zero 'x' and 'landmarks' according to 'zeroLocation'
if isfield(landmarks,'zeroLocation')
    x = x - landmarks.zeroLocation;
    landmarks.interfaceBoundaryLocation = landmarks.interfaceBoundaryLocation - landmarks.zeroLocation;
    landmarks.zeroLocation = 0;
end

% set maximum number of iterations for nonlinear fit. default is 100.
options.MaxIter = 1000;

%% ***** Fit exponential *****
% First, determine the initial guess by linear regression on
% log-transformed data
% subtract background
if isfield(offset,'value')
    yOverOffset = y-offset.value;
else
    yOverOffset = y-min(y);
end
% If there are any values <= 0, change those values to match the smallest
% positive values. This is to avoid taking the logarithm of something
% that's <= 0, which would create problems in the fitting
yOverOffset(yOverOffset<=0) = min(yOverOffset(yOverOffset>0));
yLog = log(yOverOffset);
% Obtain the initial guess by linear regression
p0Temp = polyfit(x,yLog,1);
% 'p0' specifies the initial guess. p0(1) = amplitude; p0(2) = decay 
% length; p0(3) = offset (if left as a free fitting parameter)
p0(1) = exp(p0Temp(2));
p0(2) = -1/p0Temp(1);
% If the user specifies to keep the offset as a free fitting parameter,
% then we must also set an initial guess for this parameter
if strcmp(offset.mode,'free')
    if isfield(offset,'value')
        p0(3) = offset.value;
    else
        p0(3) = min(y);
    end
end

% reset warning messages
lastwarn('');

% Fit the decaying exponential function. The field 'P' will be a vector
% specifying the values of each of the best-fit parameters, in the same
% order as specified above for 'p0'.
[fitResults.exponential.P,r,~,cov,fitResults.exponential.mse] = nlinfit(x,y,@(p,x) decayingExponential(p,x,nan,offset),p0,options);

% Determine the confidence intervals for the best-fit parameters
fitResults.exponential.ci = nlparci(fitResults.exponential.P,r,'covar',cov);

% calculate R^2
yFit = decayingExponential(fitResults.exponential.P,x,nan,offset);
SStot = sum((y-mean(y)).^2);
SSres = sum((y-yFit).^2);
fitResults.twoDomain.R2 = 1-(SSres/SStot);

% Determine whether a warning message has been returned
if strcmp(lastwarn(),'')
    fitResults.exponential.warningReturned = false;
else
    fitResults.exponential.warningReturned = true;
end

%% ***** Fit two-domain model *****
if flags.fitTwoDomainModel
    % 'p0' specifies the initial guess. p0(1) = amplitude; p0(2) = proximal
    % intrinsic decay length; p0(3) = distal intrinsic decay length; p0(4)
    % = offset (if left as a free fitting parameter)
    % As the initial guess, let's set the amplitude as the amplitude
    % estimated from the exponential fit, and set the initial value for
    % both of the intrinsic decay lengths as the decay length estimated
    % from the exponential fit
    p0(1) = fitResults.exponential.P(1);
    p0(2) = fitResults.exponential.P(2);
    p0(3) = fitResults.exponential.P(2);
    % If the user specifies to keep the offset as a free fitting parameter,
    % then we must also set an initial guess for this parameter. Take it to
    % be the offset from the best-fit exponential
    if strcmp(offset.mode,'free')
        p0(4) = fitResults.exponential.P(3);
    end
    
    % reset warning messages
    lastwarn('');
    
    % Fit the two-domain model
    %  The field 'P' will be a vector specifying the values of each of the 
    % best-fit parameters, in the same order as specified above for 'p0'.
    [fitResults.twoDomain.P,r,~,cov,fitResults.twoDomain.mse] = nlinfit(x,y,@(p,x) twoDomain(p,x,landmarks,offset),p0,options);
    
    % Determine the confidence intervals for the best-fit parameters
    fitResults.twoDomain.ci = nlparci(fitResults.twoDomain.P,r,'covar',cov);
    
    % calculate R^2
    yFit = twoDomain(fitResults.twoDomain.P,x,landmarks,offset);
    SStot = sum((y-mean(y)).^2);
    SSres = sum((y-yFit).^2);
    fitResults.twoDomain.R2 = 1-(SSres/SStot);
    
    % Determine whether a warning message has been returned
    if strcmp(lastwarn(),'')
        fitResults.twoDomain.warningReturned = false;
    else
        fitResults.twoDomain.warningReturned = true;
    end
end

%% ***** Fit two-domain-gradual-sink model *****
if flags.fitTwoDomainGradualSinkModel
    % As the initial guess, let's set the amplitude as the amplitude
    % estimated from the exponential fit, and set the initial value for
    % the proximal intrinsic decay length as the decay length 
    % estimated from the exponential fit. Let's simply set the initial
    % guess for the slope of sink onset as 100 (I wasn't able to figure
    % out a more intelligent way to determine the guess). 
    p0(1) = fitResults.exponential.P(1);
    p0(2) = fitResults.exponential.P(2);
    p0(3) = 100;
    % If the user specifies to keep the offset as a free fitting parameter,
    % then we must also set an initial guess for this parameter. Take it to
    % be the offset from the best-fit exponential
    if strcmp(offset.mode,'free')
        p0(4) = fitResults.exponential.P(3);
    end
    
    % reset warning messages
    lastwarn('');
    
    % Fit the two-domain-gradual-sink model
    %  The field 'P' will be a vector specifying the values of each of the 
    % best-fit parameters, in the same order as specified above for 'p0'.
    [fitResults.twoDomainGradualSink.P,~,residual,~,~,~,jacobian] =...
    	lsqcurvefit(@(p,x) twoDomainGradualSink(p,x,landmarks,offset),p0,x,y,[0 0 0 0],[Inf Inf Inf Inf]);

    % Determine the confidence intervals for the best-fit parameters
    fitResults.twoDomainGradualSink.ci = nlparci(fitResults.twoDomainGradualSink.P,residual,'jacobian',jacobian);
    
    % calculate R^2
    yFit = twoDomainGradualSink(fitResults.twoDomainGradualSink.P,x,landmarks,offset);
    SStot = sum((y-mean(y)).^2);
    SSres = sum((y-yFit).^2);
    fitResults.twoDomainGradualSink.R2 = 1-(SSres/SStot);
    
    % Determine whether a warning message has been returned
    if strcmp(lastwarn(),'')
        fitResults.twoDomainGradualSink.warning_returned = false;
    else
        fitResults.twoDomainGradualSink.warning_returned = true;
    end

end

end

