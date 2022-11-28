clc;clear;
% x1.,y1, y2, z1, z2
x_init=[1,1,1,1,1];
upp_bound=[10,Inf,24,10,10];
low_bound=[0,3.16,-Inf,-10,0];
[optim_vect,res_obj,exit_flag,op]=fmincon(@obj,x_init,...
    [],[],[],[],low_bound,upp_bound...
    ,@nonlincon,optimset('TolFun',1e-12,'disp','iter','PlotFcn','optimplotfval','TolX',1e-40));

feval=obj(x_init);
p1=phi_1(x_init);
p2=phi_2(x_init);
function [c,ceq]=nonlincon(x)
c=[];
ceq(1)=phi_1(x);
ceq(2)=phi_2(x);
end

%objective function
function cost=obj(v)
cost=v(1)^2+v(5)+v(2)+exp(-v(3));
end
% constraint 1
function phi= phi_1(v)
phi=v(4)^2+v(1)+v(5)-0.2*v(3)-v(2);
end
% constraint 2
function phi= phi_2(v)
phi=v(2)^(0.5)+v(4)+v(5)-v(3);
end