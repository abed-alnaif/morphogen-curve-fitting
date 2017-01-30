function y = twoDomainGradualSink(p,x,landmarks,offset)
%TWODOMAINGRADUALSINK  The decaying exponential function.
% y = TWODOMAINGRADUALSINK(p,x,parameters,offset_in)
%
% This function represents the steady-state solution to the
% two-domain-gradual-sink model. Refer to the Supplemental Modeling Notes 
% of the paper for more details on this model.
% 
% The parameters of the model are: 
%   p(1) = amplitude
%   p(2) = proximal intrinsic decay length
%   p(3) = distal gradual sink slope
%   p(4) = offset
%   landmarks.interfaceBoundaryLocation = location of interface
%	boundary
%
%
% *** INPUT ARGUMENTS ***
%
% 'p': A vector specifying the parameters of the exponential, as specified
% above. p(4) is not necessary if 'offset.mode' = 'fixed'.
% 'x': A vector indicating the points at which to evaluate the function
% 'landmarks': A structure with one field: '.interfaceBoundaryLocation'
% specifies the location of the interface boundary.
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
% % Plot the two-domain solution with 'offset.mode' = 'free', and with
% % amplitude = 1, proximal decay length = 0.5, distal gradual sink slope 
% % = 1, and offset = 0.1
% x = 0:0.01:2;
% p = [1,0.5,100,0.1];
% offset.mode = 'free';
% landmarks.interfaceBoundaryLocation = 1;
% y = twoDomainGradualSink(p,x,landmarks,offset);
% figure;
% plot(x,y);
%
% % Plot the two-domain solution with 'offset.mode' = 'fixed', and with
% % amplitude = 1, proximal decay length = 0.5, distal gradual sink slope 
% % = 1, and offset = 0.1
% x = 0:0.01:2;
% p = [1,0.5,100,0.1];
% offset.mode = 'fixed';
% offset.value = 0.1;
% landmarks.interfaceBoundaryLocation = 1;
% y = twoDomainGradualSink(p,x,landmarks,offset);
% figure;
% plot(x,y);
% 
%
% ******
% Created by Abed Alnaif, abed.alnaif@gmail.com
% Tested in Matlab R2012b
% ******

% Extract parameter values
m0 = p(1);
lambdaP = p(2);
qm = p(3);
xB = landmarks.interfaceBoundaryLocation;

% Determine value of function offset
if strcmp(offset.mode,'free')
    offsetVal = p(4);
else
    offsetVal = offset.value;
end

% Steady-state solution to two-domain-gradual-sink model
y = ...
    ... % left domain
    (  heaviside(-x+xB).*...
    ( m0.*csch(lambdaP.^(-1).*xB).*((-1).*sinh(lambdaP.^(-1).*(x+(-1).* ...
xB))+sinh(lambdaP.^(-1).*x).*(cosh(lambdaP.^(-1).*xB)+(-1).* ...
lambdaP.*qm.^(1/3).*airy(0,lambdaP.^(-2).*qm.^(-2/3)).^(-1).* ...
airy(1,lambdaP.^(-2).*qm.^(-2/3)).*sinh(lambdaP.^(-1).*xB)) ...
.^(-1)) )  )...
+ ...
... % right domain
(  heaviside(x-xB).*...
( m0.*airy(0,qm.^(-2/3).*(lambdaP.^(-2)+qm.*(x+(-1).*xB))).*(airy(0, ...
lambdaP.^(-2).*qm.^(-2/3)).*cosh(lambdaP.^(-1).*xB)+(-1).* ...
lambdaP.*qm.^(1/3).*airy(1,lambdaP.^(-2).*qm.^(-2/3)).*sinh( ...
lambdaP.^(-1).*xB)).^(-1) )  )...
...
+ offsetVal;

end