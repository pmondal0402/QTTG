function [spinor_state] = spinor1(theta,phi)
theta = theta*pi/180;
phi = phi*pi/180;

spinor_state = (cos(theta/2)*[1 0]' + exp(1i*phi)*sin(theta/2)*[0 1]');
%spinor_state = spinor_state/spinor_state'*spinor_state;


