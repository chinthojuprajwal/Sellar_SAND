clc;clear;
%% problem formulation

% the original problem was to optimise 
%         obj(x1,z1,z2)=x1^2+z2+y1+e^(-y2) w.r.t z1,z2,x1
%         subject to discipline 1: phi_1(z1,x1,z2,y2)=y1=z1^2+x1+z2-0.2y2
%          discipline 2: phi_2(y1,z1,z2)=y2=sqrt(y1)+z1+z2  and limits

% For ATC problem formulation, this is broken down ino

% sa is shared variable set [x1,z1,z2] local to optimisation 1
% sb is shared variable set [x1,z1,z2] local to optimisation 2
% sc is shared variable set [z1,z2] local to optimisation 3
% ya is coupled variable set [y1,y2] local to optimisation 1
% yb is coupled variable set [y1,y2] local to optimisation 2
% yc is coupled variable set [y1,y2] local to optimisation 3

%   optimisation 1

%   optimise obj(x1,y1,z2,y2)+L(sa-sb)+L(ya-yb)+L(sa-sc)+L(ya-yc) 
%   w.r.t sa,ya and respective limits for each variable

%   where L is the augumented lagrangian penality

%   optimisation 2

%   optimise L(phi_1(z1,x1,z2,y2)-y1)+L(sa-sb)+L(ya-yb)
%   w.r.t sb,yb and respective limits for each variable

%   where L is the augumented lagrangian penality


%   optimisation 3

%   optimise L(phi_2(y1,z1,z2)-y2)+L(sa-sc)+L(ya-yc)
%   w.r.t sc,yc and respective limits for each variable

%   where L is the augumented lagrangian penality

%   L(x-y) is defined here as v*|x-y|+w*w*|x-y|*|x-y|
%   v is updated as v=v+w*w*|x-y|
%   w is updated as w=w*beta

% initial targets
ya=[10,10];
yb=ya;yc=ya;
sa=[10,10,10];
sb=sa;sc=[10,10];
ysx=[sa,ya];
% growth factor
beta=1.05;

%initial values, limits for optimiser
xi_obj=[5,5,5,5,5];
xi_obj3=[5,5,5,5];
low_bound=[0,-10,0,3.16,-Inf];
upp_bound=[10,10,10,Inf,24];
low_boundc=[-10,0,3.16,-Inf];
upp_boundc=[10,10,Inf,24];
v1=zeros([1,2]);v2=zeros([1,2]);v3=zeros([1,3]);v4=zeros([1,2]);
w1=ones([1,2]);w2=ones([1,2]);w3=ones([1,3]);w4=ones([1,2]);
v5=zeros([1,3]);v6=zeros([1,2]);v9=zeros([1,1]);
w5=ones([1,3]);w6=ones([1,2]);w9=ones([1,1]);
v10=zeros([1,2]);v11=zeros([1,2]);v12=zeros([1,1]);
w10=ones([1,2]);w11=ones([1,2]);w12=ones([1,1]);
f_new=10000;f_old=100001;i=0;
f_obj1=1;f_obj2=1;f_obj3=1;

    scale1=obj([xi_obj(1),xi_obj(2),xi_obj(3)],[xi_obj(4),xi_obj(5)]);
    scale2=con([xi_obj(1),xi_obj(2),xi_obj(3)],sb,[xi_obj(4),xi_obj(5)],yb,v3,w3,v1,w1)+...
        con([xi_obj(2),xi_obj(3)],sc,[xi_obj(4),xi_obj(5)],yc,v4,w4,v2,w2);
    scale=scale2/scale1;
    
% ATC loop
while f_obj2+f_obj3>0.1 %f_new<f_old
    
    
    %% optimiser 1
    
    %obja is the sum of objective function, augumented langrangian terms
    %for shared and couple variable consistancies between local copies of
    %obj - discipline 1 and obj - discipline 2

    obja=@(x) 100*scale*obj([x(1),x(2),x(3)],[x(4),x(5)])+...
        con([x(1),x(2),x(3)],sb,[x(4),x(5)],yb,v3,w3,v1,w1)+...
        con([x(2),x(3)],sc,[x(4),x(5)],yc,v4,w4,v2,w2);
    %optimiser
    [x_obj1,f_obj1,ef_obj1]=fmincon(obja,xi_obj,[],[],[],[],...
        low_bound,upp_bound,[],optimset('display','off'));
    % assigning results to local copies
    sa=[x_obj1(1),x_obj1(2),x_obj1(3)];ya=[x_obj1(4),x_obj1(5)];
    %updating weights
    v1=v1+w1.*w1.*   (ya-yb);w1=w1*beta; v2=v2+w2.*w2.*   (ya-yc);w2=w2*beta;
    v3=v3+w3.*w3.*   (sa-sb);w3=w3*beta;v4=v4+w4.*w4.*   ([sa(2),sa(3)]-sc);w4=w4*beta;
    
    %% optimiser 2
    %objb is the sum of discipline 1 inequality with penality, augumented langrangian terms
    %for shared and couple variable consistancies between local copies of
    %obj - discipline 1
    objb=@(x)  con(sa,[x(1),x(2),x(3)],ya,[x(4),x(5)],v5,w5,v6,w6)...
        +v9*phi_1([x(1),x(2),x(3)],[x(4),x(5)])+...
        w9*w9*phi_1([x(1),x(2),x(3)],[x(4),x(5)])*phi_1([x(1),x(2),x(3)],[x(4),x(5)]);
    %optimiser
    [x_obj2,f_obj2,ef_obj2]=fmincon(objb,xi_obj,[],[],[],[],low_bound,...
        upp_bound,[],optimset('display','off'));
    %assigning results to local copies
    sb=[x_obj2(1),x_obj2(2),x_obj2(3)];yb=[x_obj2(4),x_obj2(5)];
    % updating penality weights
    v5=v5+w5.*w5.*   (sa-sb);w5=w5*beta; v6=v6+w6.*w6.*   (ya-yb);w6=w6*beta;
    v9=v9+w9*w9*phi_1(sb,yb);w9=w9*beta;
    
    
    %% optimiser 3
    %objb is the sum of discipline 2 inequality with penality, augumented langrangian terms
    %for shared and couple variable consistancies between local copies of
    %obj - discipline 2    
    
    objc=@(x)  con([sa(2),sa(3)],[x(1),x(2)],ya,[x(3),x(4)],v10,w10,v11,w11)...
        +v12*phi_2([x(1),x(2)],[x(3),x(4)])+...
        w12*w12*phi_2([x(1),x(2)],[x(3),x(4)])*phi_2([x(1),x(2)],[x(3),x(4)]);
    %optimiser
    [x_obj3,f_obj3,ef_obj3]=fmincon(objc,xi_obj3,[],[],[],[],low_boundc,upp_boundc,...
        [],optimset('display','off'));
    %assigning results to local copies
    sc=[x_obj3(1),x_obj3(2) ];yc=[x_obj3(3),x_obj3(4)];
    %updating weights
    v10=v10+w10.*w10.*   ([sa(2),sa(3)]-sc);w10=w10*beta; v11=v11+w11.*w11.*   (ya-yc);
    w11=w11*beta;
    v12=v12+w12*w12*phi_2(sc,yc);w12=w12*beta;
    f_old=f_new;
    f_new=f_obj1;
    i=i+1;
end
disp([sa,ya]);
%function to compute object function

function cost=obj(sa,ya)
cost=sa(1)^2+sa(3)+ya(1)+exp(-ya(2));
end

%Function to compute nonlinear inequality constraint based on discipline 1


function phi= phi_1(sb,yb)
phi=    (sb(2)^2+sb(1)+sb(3)-0.2*yb(2)-yb(1));
end

%Function to compute nonlinear inequality constraint based on discipline 2

function phi= phi_2(sc,yc)
phi=   (yc(1)^(0.5)+sc(1)+sc(2)-yc(2));
end

% function that applies augumented lagrangian to local vs other coupling
% and shared variables. Basically, consistancy constraint

function con=con(sa,sb,ya,yb,vs,ws,vy,wy)
con=sum(vy.*   (ya-yb)+(wy.*wy.*(   (ya-yb).^2)))+sum(vs.*   (sa-sb)...
    +(ws.*ws.*((sa-sb).^2)));
end
