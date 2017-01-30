function y = twoDomain(p,x,landmarks,offset)
%TWODOMAIN  The decaying exponential function.
% y = TWODOMAIN(p,x,parameters,offset_in)
%
% This function represents the steady-state solution to the two-domain
% model. Refer to the Supplemental Modeling Notes of the paper for more
% details on this model.
% 
% The parameters of the model are: 
%   p(1) = amplitude
%   p(2) = proximal intrinsic decay length
%   p(3) = distal intrinsic decay length
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
% % amplitude = 1, proximal decay length = 0.5, distal decay length = 0.1, 
% % and offset = 0.1
% x = 0:0.01:2;
% p = [1,0.5,0.1,0.1];
% offset.mode = 'free';
% landmarks.interfaceBoundaryLocation = 1;
% y = twoDomain(p,x,landmarks,offset);
% figure;
% plot(x,y);
%
% % Plot the two-domain solution with 'offset.mode' = 'fixed', and with
% % amplitude = 1, proximal decay length = 0.5, distal decay length = 0.1, 
% % and offset = 0.1
% x = 0:0.01:2;
% p = [1,0.5,0.1];
% offset.mode = 'fixed';
% offset.value = 0.1;
% landmarks.interfaceBoundaryLocation = 1;
% y = twoDomain(p,x,landmarks,offset);
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
    offsetVal = p(4);
else
    offsetVal = offset.value;
end

% Steady-state solution to two-domain model
y = ( (p(1).*(exp((-x+landmarks.interfaceBoundaryLocation)./p(3)).*p(3).*heaviside(x-landmarks.interfaceBoundaryLocation) + ...
(p(3).*cosh((x - landmarks.interfaceBoundaryLocation)./p(2)) - p(2).*sinh((x - landmarks.interfaceBoundaryLocation)./p(2))).*heaviside(-x+landmarks.interfaceBoundaryLocation)))./...
(p(3).*cosh(landmarks.interfaceBoundaryLocation./p(2)) + p(2).*sinh(landmarks.interfaceBoundaryLocation./p(2))) ) + offsetVal;

end