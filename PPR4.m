%Initial Guess
x0= [2,2,2,2,2,2];

%parameters
%{
    1.  1.65       w_bp, l_bp, l_f  [in]
    2.  0.16       h_bp             [in]
    3.  0.44       d_bp             [in]
    4.  25         k                [lb/in]
    5.  1.5        nd               [unitless]
    6.  40000      Sy 6061 Al       [psi]
    7.  7250       Sy PLA           [psi]
    8.  0.1484     delta x          [in]
    9.  0.8516     spring length    [in]
    10. 0.13       diameter screw   [in]
    11. 0.51       diameter pulley  [in]
    12. 30000    Ss 6061 Al       [psi]
%}

p = [1.65, 0.16, 0.44, 25, 1.5, 40000, 7250, 0.1484, 0.8516, 0.13, 0.51, 30000];

%Bound Definitions
lb = [.0001 .0001 .63 p(10)/2 .0001 .0001]; 
ub = [.5 .5 Inf Inf Inf Inf];

options = optimset('Algorithm', 'active-set');
[xopt, fopt] = fmincon(@ExtruderNonLinF, x0, [], [], [], [], lb, ub, @ExtruderNonLinCon, options, p);
disp(fopt)
disp(xopt)

function f = ExtruderNonLinF(x, p)
    % Objective Function
    d1 = x(1); d2 = x(2); l1 = x(3); l2 = x(4); wl = x(5); wf = x(6);
    f = (wl.*wf.*p(1)) + (l1.*d1 + l2.*d2 + d1.*d2).*wl + (p(1).^2) * p(2) - (pi * p(3).^2 * p(2));
end 

function [C, Ceq] = ExtruderNonLinCon(x, p)

    %Constraint Function
    d1 = x(1); d2 = x(2); l1 = x(3); l2 = x(4); wl = x(5); wf = x(6);
    
    %%%Intermediate Calculations
    j10=p(4)*p(8); %spring force [lb]
    j11=j10*(2.*l1+d2)./(2.*l2+d1); %force on filament [lb]
    j12=pi/4*.06889764.^2; %cross sectional area of filament [in^2]
    j13=(3/2).*(j10./(wl.*d1)); %Shear stress on Lever
    j14=(3/2).*(j10./(wf.*wl)); %Shear stress on Flange
    j6=sqrt((3.*j10.*wl./(p(1).*wf.^3).^2)+ (18.*j10.*(p(2)+wl)./(wf.^2)).^2); %sigma f
    j8=6.*j10.*(l1)./(wl.*d1.^2); %sigma l, xx
    j9=6*j10.*l2.*(2*l1+d2)./(l2.*d2.^2.*(2.*l2+d1)); %sigma l, yy
    j7=sqrt(.5.*((j8-j9).^2)+(j8.^2)+j9.^2); %sigma l
    
    %Constraints
    C(1)=j6-p(6)/p(5);
    C(2)=j7-p(6)/p(5);
    C(3)=j11/j12-p(7);
    C(4)=j13-(p(12)/p(5));
    C(5)=j14-(p(12)/p(5));
    Ceq=[];
end