function pi=obj_function(X0,rhol,vl,ql,Ll,Loff,lambdal,lambdaoff,d,beta,w,rhooff,...
    Qc,von,Np,rhomax,rhocrit,tau,kappa,theta,... %qin% %qout%
    phir,phib,phiw,vf,alpha,A,E,T)
    
    %split input matrix
    [r,c]=size(X0);
    Nc=r/2;
    r=X0(1:9,1:Nc);
    b=X0(1:20,Nc+1:c);
    
    %% initial values
    % 1 for k and 2 for k+1
    pi=0;
    rhol1=rhol;
    vl1=vl;
    qin1=zeros(1,9);
    %qout1=qout;
    w1=w;
    rhooff1=rhooff;
    von1=von;
    ql1=ql;
    rhol2=zeros(1,20);
    vl2=zeros(1,20);
    ql2=zeros(1,20);
    %qin2=zeros(9,1);
    %qout2=zeros(14,1);
    w2=zeros(1,13);
    rhooff2=zeros(1,14);
    Qo=zeros(1,13);
    
    onramp_next_link=[0;3;0;5;6;9;10;12;13;16;0;18;0];
    %offramp_next_link=[1;0;1;0
    
    node_in_out=[...    %in out1 out2 out3 in_link_no
        0,1,0,0,1;...   %1
        1,0,0,0,2;
        0,2,0,0,3;
        2,0,0,0,4;
        3,3,0,0,5;      %5
        4,4,0,0,8;
        5,9,10,11,9;
        6,12,13,14,11;
        7,7,0,0,12;
        8,8,0,0,15;     %10
        0,9,0,0,16;
        9,0,0,0,17;
        0,10,0,0,18;];  %13
    
    link_node=[...  %prev, next
        -1,1;       %-1=node not exist
        1,2;        %0=no node between this link and the next
        2,3;
        3,4;
        4,5;    %5
        5,0;
        0,0;
        0,6;
        6,7;
        7,-1;   %10
        -1,8;
        8,9;
        9,0;
        0,0;
        0,10;   %15
        10,11;
        11,12;
        12,13;
        13,0;
        0,-1];  %20
        
    
    %% iterations
    for i=1:(Np-1)
       
        %% Qo
        for j=1:13
            Qo=ql1(node_in_out(j,5));
            if node_in_out(j,1)~=0
                Qo=Qo+qin1(node_in_out(j,1));
            end
            %for k=2:4
            %   if node_in_out(j,k)~=0
            %       Qo=Qo+qout1(node_in_out(j,k));
            %   end
            %end
        end
        
       %% On-ramps
        for j=1:13  %each node
            if d(j)~=0
                qo1=d(j)+w1(j)/T;
                link_id=onramp_next_link(j);
                qo2=Qc(link_id)*min(1,(rhomax-rhol1(link_id))...
                    /(rhomax-rhocrit));
                qin1(link_id)=r(j,min(i,Nc))*min(qo1,qo2);
                w2(j)=w1(j)+T*(d(j)-qin1(link_id));
            end
        end
        
       %% Off-ramps
        for j=1:13  %each node
            for k=2:4
                offramp_no=node_in_out(j,k);
                if offramp_no~=0
                    rhooff2(offramp_no)=rhooff1(offramp_no)+...
                        T*Qo(j)*beta(j)/...
                        (Loff(offramp_no)*lambdaoff(offramp_no));
                end
            end
        end
        
       %% Main Links
        for j=1:20  %each link
            %rho
            prev_node=link_node(j,1);
            next_node=link_node(j,2);
            if prev_node==-1    %prev node not exist
                tmp=0;
            elseif prev_node==0 %no nodes between prev link and self
                tmp=ql1(j-1)-ql1(j);
            else
                tmp=Qo(prev_node)*(1-beta(prev_node));
            end
            rhol2(j)=rhol2(j)+tmp;
            
            %v
            %v1
            tmp_b=b(j,min(i,Np-1));
            tmp_v_f=vf*tmp_b;
            tmp_rho_crit=rho_crit*(1+A*(1-tmp_b));
            tmp_alpha=alpha*(E-(E-1)*tmp_b);
            tmp_v1=T*(tmp_v_f*exp(-(1/tmp_alpha)*...
                ((rhol1(j)/tmp_rho_crit).^tmp_alpha)))/tau;
            %v2
            if prev_node==-1
                tmp_v2=0;
            elseif prev_node==0
                tmp_v2=(T/L(j))*vl1(j)*(vl1(j-1)-vl1(j));
            else
                on_link=node_in_out(prev_node,1);
                virtual_a=vl1(j-1)*ql1(j-1)...
                    +qin1(on_link)*von1(on_link);    %Numerator for virtual_v
                virtual_b=ql1(j-1)+qin1(on_link);    %Denominator for virtual_v
                virtual_v=virtual_a/virtual_b;
                tmp_v2=(T/L(j))*vl1(j)*(virtual_v-vl1(j));
            end
            %v3
            if next_node==-1
                tmp_v3=0;
            elseif next_node==0
                tmp_v3=-((theta*T)/(tau*Ll(j-1)))*...
                    (rho1(j+1)-rho1(j))/(rho1(j)+kappa);
            else
                virtual_a=rho1(j+1).^2;
                virtual_b=rho1(j+1);
                for k=2:4
                    out_link=node_in_out(prev_node,k);
                    if out_link~=0
                       virtual_a=virtual_a+rhooff1(out_link).^2;
                       virtual_b=virtual_b+rhooff1(out_link);
                    end
                end
                tmp_v3=virtual_a/virtual_b;
            end
            vl2(j)=vl1(j)+tmp_v1+tmp_v2+tmp_v3;
           
            %q
            ql2(j)=rhol2(j)*vl2(j)*lambdal(j);

        end
        
       %% Travel Time & Queue Length
        pi=pi+sum(rhol2.*Ll.*lambdal)+sum(w2);
        %w
        dw=zeros(9,1);
        for j=1:9
            dw(j)=(w2(j)-wmax(j))^2;
        end
        pi=pi+phiw*sum(dw);
        
       %% Discard old state
        rhol1=rhol;
        vl1=vl2;
        %qin1=qin2;
        %qout1=qout2;
        w1=w2;
        rhooff1=rhooff2;
        von1=von2;
        ql1=ql2;
        rhol2=zeros(1,20);
        vl2=zeros(1,20);
        ql2=zeros(1,20);
        %qin2=zeros(9,1);
        %qout2=zeros(14,1);
        w2=zeros(1,13);
        rhooff2=zeros(1,14);
        Qo=zeros(1,13);
    end
    
    %b
    db=zeros(20,Np-1);
    for i=1:20
       for j=1:Np-1
          db(i,j)=(b(i,j+1)-b(i.j)).^2; 
       end
    end
    
    %r
    dr=zeros(9,Np-1);
    for i=1:9
       for j=1:Np-1
          dr(i,j)=(r(i,j+1)-r(i.j)).^2; 
       end
    end
    
    pi=pi+phir*sum(dr)+phib*sum(db);
    
end