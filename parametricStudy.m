clear;clc;close all;
%parameters
%{
    1.  1.65       w_bp, l_bp, l_f  [in]
    2.  0.16       h_bp             [in]
    3.  0.44       d_bp             [in]
    4.  25         k                [lb/in]
    5.  1.5        nd               [unitless]
    6.  40000      Sy               [psi]
    7.  7250       Sy PLA           [psi]
    8.  0.1484     delta x          [in]
    9.  0.8516     spring length    [in]
    10. 0.13       diameter screw   [in]
    11. 0.51       diameter pulley  [in]
    12. 30000      Ss               [psi]
%}
sample_size = 25;
p = [1.65, 0.16, 0.44, 25, 1.5, 40000, 7250, 0.1484, 0.8516, 0.13, 0.51, 30000, 10000000];
result = zeros([8 sample_size]);
tmp = linspace(0,0,12);
max = 200;
min = -95;
change = linspace(min, max, sample_size);

for j = 1:sample_size
    tmp = p(:);
    tmp(4) = p(4)*(1+change(j)/100);
    result(1, j) = ParametricStudy(tmp);
end

for i = 2:9
    for j = 1:sample_size
        tmp = p(:);
        tmp(i+4) = p(i+4)*(1+change(j)/100);
        result(i, j) = ParametricStudy(tmp);
    end
end

hold on
for i = 1:8
    plot(change, result(i, :), '-');
end

lgd = legend('k', 'S_y', 'S_y PLA ', '\Deltax ', 'l_s_p_r_i_n_g', 'd_s_c_r_e_w', 'd_p_u_l_l_e_y', 'S_s');
lgd.FontSize = 20;
xlabel('% change', 'FontSize', 20)
ylabel('volume', 'FontSize', 20)
title('parametric study', 'FontSize', 20)

figure(2);
hold on;
plot(change, result(1, :), '->');
plot(change, result(2, :), '-+');
plot(change, result(8, :), '-*');
plot(change, result(9, :), '-o');
lgd = legend('Spring Constant', 'Yield Strength', 'Shear Strength', 'Elastic Modulus');
lgd.FontSize = 20;
xlabel('% change', 'FontSize', 20)
ylabel('volume', 'FontSize', 20)
title('parametric study refined', 'FontSize', 20)

for i = 1:8
    fprintf("\n\nparameter %d:\n\n", i);
    for j = 1:sample_size
        fprintf("%.4f \t", result(i,j));
    end
end

function obj = ParametricStudy(p)
    %Initial Guess
    x0= [1,1,1,1,1,1]; 
    
    %Bound Definitions
    lb = [.0001 .0001 .63 0.51/2 .13 .0001]; 
    ub = [Inf Inf Inf Inf Inf Inf];
    
    options = optimset('Algorithm', 'active-set');
    [xopt, fopt] = fmincon(@ExtruderNonLinF, x0, [], [], [], [], lb, ub, @ExtruderNonLinCon, options, p);
    clc;
    if fopt < 0
        fopt = 0;
    end
    obj = fopt;
end
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
    j15=10; %fuser lb
    j16=wl.*(d1.^3)/12; %MOI Lever
    j17=wl.*(wf.^3)/12; %MOI Flange
    j19=l1+d2; %Lever length
    j18=j19./10; %b
    j20=j19-j18; %a

    %Constraints
    C(1)=j6-p(6)/p(5);
    C(2)=j7-p(6)/p(5);
    C(3)=j11/j12-p(7);
    C(4)=j13-(p(12)/p(5));
    C(5)=j14-(p(12)/p(5));
    C(6)=(j15.*j18.*j20).*(j19.^2-j18.^2-j20.^2)./(6.*p(13).*j17.*j19)-.0001;
    C(7)=(j15.*j19.^3)./(3.*p(13).*j16)-.0001;
    Ceq=[];
end