function h = ODE_h(in1,u1,in3,alpha)
%ODE_H
%    H = ODE_H(IN1,U1,IN3,ALPHA)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    13-Jan-2019 22:22:05

lambda1 = in1(3,:);
x2 = in1(2,:);
h = [x2;u1;0.0;-lambda1];
