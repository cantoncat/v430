
%original objective function
%output V is suppressed by wrapper() in fmincin()
%output pi is suppressed while estimating future traffic conditions

function [pi,V]=obj_function(X0,rhol,vl,ql,Ll,Loff,lambdal,lambdaoff,d,beta,w,rhooff,...
    Qc,von,Np,Nc,rhomax,rhocrit,tau,kappa,theta,wmax,... %qin% %qout%
    phir,phib,phiw,vf,alpha,A,E,T)
    
    %split input matrix
    %[r,c]=size(X0);
    %Nc=r/2;
    r=X0(1:9,1:Nc);
    b=X0(1:20,Nc+1:(2*Nc));
    
    %% initial values
    V=zeros(Np-1,20);
    V=[vl;V];
    % 1 for k and 2 for k+1 except for V the matrix
    pi=0;
    rhol1=rhol;
    %vl1=vl;
    qin1=zeros(1,9);
    %qout1=qout;
    w1=w;
    rhooff1=rhooff;
    von1=von;
    ql1=ql;
    rhol2=zeros(1,20);
    %vl2=zeros(1,20);
    ql2=zeros(1,20);
    %qin2=zeros(9,1);
    %qout2=zeros(14,1);
    w2=zeros(1,13);
    rhooff2=zeros(1,14);
    Qo=zeros(1,13);
    %qoff1=qoff;
    %qoff2=zeros(1,14);
    %voff1=voff;
    %voff2=zeros(1,14);
    
    onramp_next_link=[0;3;0;5;6;9;10;12;13;16;0;18;0];
    %link_last_onramp=[];
        %0;0;1;0;2;...
        %3;0;0;4;5;...
        %0;6;7;0;0;...
        %8;0;9;0;0]
    
    node_in_out=[...    %in out1 out2 out3 in_link_no
        0,1,0,0,1;...   %1
        1,0,0,0,2;
        0,2,0,0,3;
        2,0,0,0,4;
        3,3,0,0,5;      %5
        4,4,0,0,8;
        5,9,10,11,9;
        6,12,13,14,11;
        7,5,0,0,12;
        8,6,0,0,15;     %10
        0,7,0,0,16;
        9,0,0,0,17;
        0,8,0,0,18;];  %13
    
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
            Qo(j)=ql1(node_in_out(j,5));
            if node_in_out(j,1)~=0
                Qo(j)=Qo(j)+qin1(node_in_out(j,1));
            end
        end
        
       %% On-ramps
        for j=1:13  %each node
            if d(j)~=0
                qo1=d(j)+w1(j)/(T/3600);
                link_id=onramp_next_link(j);
                qo2=Qc(node_in_out(j,1))*min(1,(rhomax-rhol1(link_id))...
                    /(rhomax-rhocrit));
                qin1(link_id)=r(node_in_out(j,1),min(i,Nc))*min(qo1,qo2);
                w2(j)=w1(j)+T*(d(j)-qin1(link_id))/3600;
            end
        end
        
       %% Off-ramps
        for j=1:13  %each node
            for k=2:4
                offramp_no=node_in_out(j,k);
                if offramp_no~=0
                    %rhooff2(offramp_no)=rhooff1(offramp_no)+...
                        %(T/3600)*(Qo(j)*beta(j)-qoff1(offramp_no))/...
                        %(Loff(offramp_no)*lambdaoff(offramp_no));
                    rhooff2(offramp_no)=rhooff1(offramp_no)+...
                        (T/3600)*(Qo(j)*beta(j)-...
                        rhooff1(offramp_no)*60*lambdaoff(offramp_no))/...
                        (Loff(offramp_no)*lambdaoff(offramp_no));
                    %rhooff2(offramp_no)=min(max(rhooff2(offramp_no),0),...
                        %rhomax);
                    %tmp_v1=vf*exp(-(1/alpha)*power((rhooff1(offramp_no)/rhocrit),alpha));
                    %tmp_v1=(T/tau)*(tmp_v1-voff(offramp_no));
                    %prev_link=node_in_out(j,5);
                    %tmp_v2=(T/Loff(offramp_no))*(V(i,prev_link)-voff(offramp_no));
                    %voff2(offramp_no)=voff1(offramp_no)+tmp_v1+tmp_v2;
                    %qoff2(offramp_no)=rhooff2(offramp_no)*...
                        %lambdaoff(offramp_no)*voff2(offramp_no);
                    %qoff2(offramp_no)=min(max(qoff2(offramp_no),0),...
                        %2300*lambdaoff(offramp_no));
                    %rhooff2(offramp_no)=rhooff1(offramp_no)+...
                        %(T/3600)*Qo(j)*beta(j)/...
                        %(Loff(offramp_no)*lambdaoff(offramp_no));
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
                tmp=Qo(prev_node)*(1-beta(prev_node))-ql1(j);
            end
            rhol2(j)=rhol1(j)+((T/3600)/(Ll(j)*lambdal(j)))*tmp;
            rhol2(j)=min(max(rhol2(j),0),rhomax);
            
            %v
            %v1
            tmp_b=b(j,min(i,Nc));
            tmp_v_f=vf*tmp_b;
            tmp_rho_crit=rhocrit*(1+A*(1-tmp_b));
            tmp_alpha=alpha*(E-(E-1)*tmp_b);
            tmp_v1=T*(tmp_v_f*exp(-(1/tmp_alpha)*...
                ((rhol1(j)/tmp_rho_crit).^tmp_alpha))-V(i,j))/tau;
            %tmp_v1=vf*exp(-(1/alpha)*(power((rhol1(j)/rhocrit),alpha)));
            %tmp_v1=(T/tau)*(min(tmp_v1,120*b(j,min(i,Nc)))-V(i,j));
            %v2
            if prev_node==-1
                tmp_v2=0;
            elseif prev_node==0
                tmp_v2=((T/3600)/Ll(j))*V(i,j)*(V(i,j-1)-V(i,j));
            else
                on_link=node_in_out(prev_node,1);
                if ((on_link~=0) && (on_link~=-1))
                    virtual_a=V(i,j-1)*ql1(j-1)...
                        +max(qin1(on_link)*von1(on_link),0);    %Numerator for virtual_v
                    virtual_b=ql1(j-1)+max(qin1(on_link),0);    %Denominator for virtual_v
                    virtual_v=virtual_a/virtual_b;
                    tmp_v2=((T/3600)/Ll(j))*V(i,j)*(virtual_v-V(i,j));
                else
                    % last node does not have an on-ramp
                    tmp_v2=((T/3600)/Ll(j))*V(i,j)*(V(i,j-1)-V(i,j));
                end
            end
            %v3
            if next_node==-1
                tmp_v3=0;
            elseif next_node==0
                tmp_v3=-((theta*T)/(tau*Ll(j)))*...
                    (rhol1(j+1)-rhol1(j))/(rhol1(j)+kappa);
            else
                virtual_a=power((rhol1(j+1)),2);
                virtual_b=rhol1(j+1);
                for k=2:4
                    out_link=node_in_out(next_node,k);
                    if ((out_link~=0) && (out_link~=-1))
                       virtual_a=virtual_a+power((max(rhooff1(out_link),0)),2);
                       virtual_b=virtual_b+max(rhooff1(out_link),0);
                    end
                end
                if virtual_b~=0
                    virtual_rho=virtual_a/virtual_b;
                    virtual_rho=min(max(virtual_rho,0),rhomax);
                else
                    virtual_rho=rhol1(j+1);
                end
                tmp_v3=-((theta*T)/(tau*Ll(j)))*...
                    (virtual_rho-rhol1(j))/(rhol1(j)+kappa);
            end
            V(i+1,j)=min(max(V(i,j)+tmp_v1+tmp_v2+tmp_v3,0),120);
           
            %q
            ql2(j)=rhol2(j)*V(i+1,j);
            ql2(j)=min(max(ql(2),0),2300)*lambdal(j);

        end
        
       %% Travel Time & Queue Length
        pi=pi+sum(rhol2.*Ll.*lambdal)+sum(w2);
        %w
        dw=zeros(13,1);
        for j=1:13
            dw(j)=power((max((w2(j)-wmax(j)),0)),2);
        end
        pi=pi+phiw*sum(dw);
        
       %% Discard old state
        %except V - keeping them for travel time estimation
        rhol1=rhol;
        %vl1=vl2;
        %qin1=qin2;
        %qout1=qout2;
        w1=w2;
        rhooff1=rhooff2;
        %von1=von2;
        ql1=ql2;
        %qoff1=qoff2;
        rhol2=zeros(1,20);
        %vl2=zeros(1,20);
        ql2=zeros(1,20);
        %qin2=zeros(9,1);
        %qout2=zeros(14,1);
        w2=zeros(1,13);
        rhooff2=zeros(1,14);
        Qo=zeros(1,13);
        %qoff2=zeros(1,14);
    end
    
    %b
    db=zeros(20,Nc-1);
    for i=1:20
       for j=1:Nc-1
          db(i,j)=power((b(i,j+1)-b(i,j)),2); 
       end
    end
    
    %r
    dr=zeros(9,Nc-1);
    for i=1:9
       for j=1:Nc-1
          dr(i,j)=power((r(i,j+1)-r(i,j)),2); 
       end
    end
    
    pi=pi+phir*sum(sum(dr))+phib*sum(sum(db));
    
end