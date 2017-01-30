function y = decayingExponential(p,x,~,offset)
%DECAYINGEXPONENTIAL  The decaying exponential function.
% y = DECAYINGEXPONENTIAL(p,x,parameters,offset_in)
%
% This function is used to model the decaying exponential function:
%      y = p(1)*exp(-x/p(2))+p(3)
% where x represents the independent variable, p(1) represents the
% amplitude, p(2) represents the decay length, and p(3) represents the
% offset
%
%
% *** INPUT ARGUMENTS ***
%
% 'p': A vector specifying the parameters of the exponential, as specified
% above. p(3) is not necessary if 'offset.mode' = 'fixed'.
% 'x': A vector indicating the points at which to evaluate the function
% '~': A placeholder variable so that format of function stays
% consistent with other functions. It can be simply specified as 'NaN',
% since this function does not depend on any landmarks.
% 'offset': A structure specifying how to account for the function's offset.
% It contains two fields:
%   1. '.mode': Set as 'free' for the offset term to be represented by p(3).
%   Set as 'fixed' for it to be specified by '.value'
%   2. '.value': the value of the offset term if '.mode' = '.fixed'. If
%   '.mode' = 'free', '.value' need not be a field.
%
%
% *** OUTPUT ARGUMENTS ***
%
% 'y': A vector, the same size as 'x', with the value of the function
% evaluated at each element in 'x'.
%
%
% *** EXAMPLES ***
%
% % Plot the decaying exponential with 'offset.mode' = 'free', and with
% amplitude = 1, decay length = 0.2, and offset = 0.1
% x = 0:0.01:1;
% p = [1,0.2,0.1];
% offset.mode = 'free';
% y = decayingExponential(p,x,NaN,offset);
% figure;
% plot(x,y);
%
% % Plot the decaying exponential with 'offset.mode' = 'fixed', and with
% amplitude = 1, decay length = 0.2, and offset = 0.1
% x = 0:0.01:1;
% p = [1,0.2];
% offset.mode = 'fixed';
% offset.value = 0.1;
% y = decayingExponential(p,x,NaN,offset);
% figure;
% plot(x,y);
% 
%
% ******
% Created by Abed Alnaif, abed.alnaif@gmail.com
% Tested in Matlab R2012b
% ******

% Determine value of function offset
if strcmp(offset.mode,'free')
   offsetVal = p(3);
else
   offsetVal = offset.value;
end

% Decaying exponential function
y = p(1)*exp(-x/p(2)) + offsetVal;

end