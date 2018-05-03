
%% Construct Variant Matrix X0
b0=b0';
r0=r0';
r1=ones(9,Nc-1);
for i=1:Nc-1
   r1(:,i)=r1(:,i).*r0;
end
r1=[r0,r1];
r1=[r1;zeros(11,Nc)];
b1=ones(20,Nc-1);
for i=1:Nc-1
    b1(:,i)=b1(:,i).*b0;
end
b1=[b0,b10];
X0=[r1,b1];

%% Construct UB and LB
ub=ones(9,Nc);      %r
ub=[ub;zeros(11,Nc)];
ub=[ub,ones(20,Nc)];

lb_r=ones(9,Nc);
lb_r=lb_r*rmin;
lb_r=[lb_r;zeros(11,Nc)];
lb_b=ones(20,Nc);
lb_b=lb_b*bmin;
lb=[lb_r,lb_b];

%% Optimization
[X,fval,exitflag,output]=fmincon(...
    @(X0) obj_function(...
        X0,rhol,vl,ql,Ll,Loff,lambdal,lambdaoff,d,beta,w,rhooff,...
        Qc,von,Np,rhomax,rhocrit,tau,kappa,theta,...
        phir,phib,phiw,vf,alpha,A,E,T),...
    [],[],[],[],lb,ub);

%% Split Matrix
b0=X(1:9,1:Nc);
r0=X(1:20,(Nc+1):(2*Nc));