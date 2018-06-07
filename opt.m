function [b1,r1,V,qo,fval,exitflag,output,pi,pi_tt,pi_w]=opt(...
        b0,r0,rhol,vl,ql,Ll,Loff,...
        lambdal,lambdaoff,d,beta,w,rhooff,...
        Qc,von,Np,rhomax,rhocrit,tau,kappa,theta,wmax,qin,...
        phir,phib,phiw,vf,alpha,A,E,T,Nc,rmin,bmin,s_intv,rmin_new,bmin_new)

    r_con=[];
    for i=1:9
        if rmin_new(i)==rmin
            r_con=[r_con,i]; 
        end
    end
    [~,nr]=size(r_con);

    b_con=[];
    for i=1:20
        if bmin_new(i)==bmin
            b_con=[b_con,i];
        end
    end
    [~,nb]=size(b_con); 
    
    %% Construct Variant Matrix X0
    b0=b0';
    r0=r0';
    width=T*Nc/s_intv+1;
    r1=ones(9,width-1);
    for i=1:width-1
        for j=1:9
            r1(j,i)=r1(j,i)*r0(j);
        end
    end
    r1=[r0,r1];
    r2=[];
    for i=1:nr
        id=r_con(i);
        r2=[r2;r1(id,:)];
    end
    r2=r2(:);

    b1=ones(20,width);
    for i=1:width
        for j=1:20
            b1(j,i)=(b1(j,i)*b0(j)*120-60)/5;
        end
    end
    b2=[];
    for i=1:nb
        id=b_con(i);
        b2=[b2;b1(id,:)];
    end
    b2=b2(:);
    X0=[r2;b2];

    %% Construct UB and LB

    ub_r=ones(nr,width);      %r
    ub_r=ub_r(:);
    ub_b=ones(nb,width);
    ub_b=ub_b*12;
    ub_b=ub_b(:);
    ub=[ub_r;ub_b];

    lb_r=ones(nr,width);
    for i=1:width
        for j=1:nr
            lb_r(j,i)=lb_r(j,i)*rmin; 
        end
    end

    lb_b=ones(nb,width);
    for i=1:width
        for j=1:nb
            lb_b(j,i)=lb_b(j,i)*(120*bmin-60)/5;
        end
    end
    lb=[lb_r(:);lb_b(:)];

    %% Optimization

    IntCon=(nr*width+1):((nb+nr)*width);
    fitnessfcn=@(X) wrapper(X,rhol,vl,ql,Ll,Loff,lambdal,...
                    lambdaoff,d,beta,w,rhooff,...
                    Qc,von,Np,Nc,rhomax,rhocrit,tau,kappa,theta,wmax,qin,... 
                    phir,phib,phiw,vf,alpha,A,E,T,s_intv,r_con,b_con);
    nvars=(nr+nb)*width;
    opt=gaoptimset('Display','iter','PopulationSize',80,'Generations',60,...
                    'CrossoverFraction',0.9);
    nonlcon=@(X) mynonlcon(X,X0,width,nr,nb);
        
    [X,fval,exitflag,output]=ga(fitnessfcn,nvars,[],[],[],[],...
                            lb,ub,nonlcon,IntCon,opt);

    %% Split Matrix

    [pi,pi_tt,pi_w,V,qo]=obj_function(...
            X,rhol,vl,ql,Ll,Loff,lambdal,lambdaoff,d,beta,w,rhooff,...
            Qc,von,Np,Nc,rhomax,rhocrit,tau,kappa,theta,wmax,qin,... 
            phir,phib,phiw,vf,alpha,A,E,T,s_intv,r_con,b_con);
    
    rr=X(1:nr);
    rr=rr';
    bb=X(nr*width+1:nr*width+nb);
    for i=1:nb
        bb(i)=(60+5*bb(i))/120;
    end
    r1=ones(1,9);
    for i=1:nr
        id=r_con(i);
        r1(id)=rr(i);
    end
    b1=ones(1,20);
    for i=1:nb
        id=b_con(i);
        b1(id)=bb(i);
    end
end

function [c,ceq]=mynonlcon(X,X0,width,nr,nb)
    ceq=[];
    c=zeros((nr+nb),1);
    for i=1:nr
        c(i)=abs(X0(i)-X(i))-0.2;
    end
    for i=(nr*width+1):(nr*width+nb)
        c(i)=abs(X0(i)-X(i))-2;
    end
end