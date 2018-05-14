
%This function is a wrapper of the original obj_function for fmincon()

function pi=wrapper(X0,rhol,vl,ql,Ll,Loff,lambdal,lambdaoff,d,beta,w,rhooff,...
    Qc,von,Np,Nc,rhomax,rhocrit,tau,kappa,theta,wmax,... %qin% %qout%
    phir,phib,phiw,vf,alpha,A,E,T)

    [pi,~,~]=obj_function(X0,rhol,vl,ql,Ll,Loff,lambdal,lambdaoff,d,beta,w,rhooff,...
        Qc,von,Np,Nc,rhomax,rhocrit,tau,kappa,theta,wmax,... %qin% %qout%
        phir,phib,phiw,vf,alpha,A,E,T);
    
end