
%original objective function
%output V,qo,pi_tt,pi_w is suppressed by wrapper() in fmincin()
%output qo is suppressed by wrapper() in fmincin()

function [pi,pi_tt,pi_w,V,qo]=obj_function(X0,rhol,vl,ql,Ll,Loff,lambdal,lambdaoff,d,beta,w,rhooff,...
    Qc,von,Np,Nc,rhomax,rhocrit,tau,kappa,theta,wmax,qin,... %qout%
    phir,phib,phiw,vf,alpha,A,E,T,s_intv,r_con,b_con)
    
    %split input matrix
    %[r,c]=size(X0);
    %Nc=r/2;1,1
    
    vmin=10;
    delta=0.0122;
    phi=1.99;
    
    [~,nr]=size(r_con);
    [~,nb]=size(b_con);
    
    width=T*Nc/s_intv+1;
    rr=reshape(X0(1:nr*width),[nr,width]);
    bb=reshape(X0(nr*width+1:(nr*width+nb*width)),[nb,width]);
    
    rrr=ones(9,width);
    bbb=ones(20,width);
    idi=1;
    for i=1:9
       if ismember(i,r_con)==1
          for j=1:width
             %id=r_con(idi);
             rrr(i,j)=rr(idi,j);
          end
          idi=idi+1;
       end
    end
    
    idi=1;
    for i=1:20
       if ismember(i,b_con)==1
          for j=1:width
             %id=r_con(idi);
             bbb(i,j)=(60+5*bb(idi,j))/120;
          end
          idi=idi+1;
       end
    end
   
    origin_width=s_intv/T;
    r=ones(9,Nc);
    for i=1:width
       for j=0:origin_width-1
           r(1:9,1+origin_width*(i-1)+j)=rrr(1:9,i);
       end
    end
    b=ones(20,Nc);
    for i=1:width
       for j=0:origin_width-1
           b(1:20,1+origin_width*(i-1)+j)=bbb(1:20,i);
       end
    end
    
    
    %r=X0(1:9,1:Nc);
    %b=X0(1:20,Nc+1:(2*Nc));
    
    %% initial values
    pi_tt=0;
    pi_w=0;
    qo=zeros(1,13);
    V=zeros(Np-1,20);
    V=[vl;V];
    % 1 for k and 2 for k+1 except for V the matrix
    pi=0;
    rhol1=rhol;
    %vl1=vl;
    qin1=qin;
    %qout1=qout;
    w1=w;
    rhooff1=rhooff;
    von1=von;
    ql1=ql;
    rhol2=zeros(1,20);
    %vl2=zeros(1,20);
    ql2=zeros(1,20);
    qin2=zeros(9,1);
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
        %for j=1:13
        %for j=8:13
        for j=1:7
            Qo(j)=ql1(node_in_out(j,5));
            if node_in_out(j,1)~=0
                Qo(j)=Qo(j)+qin1(node_in_out(j,1));
            end
        end
        
       %% On-ramps
        %for j=1:13  %each node
        %for j=8:13
        for j=1:7
            if d(j)~=0
                qo1=d(j)+w1(j)/(T/3600);
                link_id=onramp_next_link(j);
                ramp_id=node_in_out(j,1);
                qo2=Qc(ramp_id)*min(1,(max(rhomax-rhol1(link_id),...
                    0))/(rhomax-rhocrit));
                if i==1
                    qo(j)=r(ramp_id,min(i,Nc))*min(qo1,qo2);
                end
                qin2(ramp_id)=r(node_in_out(j,1),min(i,Nc))*min(qo1,qo2);
                w2(j)=max(w1(j)+T*(d(j)-qin2(ramp_id))/3600,0);
            end
        end
        
       %% Off-ramps
        %for j=1:13
        %for j=8:13  %each node
        for j=1:7
            for k=2:4
                offramp_no=node_in_out(j,k);
                if offramp_no~=0
                    rhooff2(offramp_no)=rhooff1(offramp_no)+...
                        (T/3600)*(Qo(j)*beta(j)-...
                        rhooff1(offramp_no)*60*lambdaoff(offramp_no))/...
                        (Loff(offramp_no)*lambdaoff(offramp_no));
                    rhooff2(offramp_no)=min(max(rhooff2(offramp_no),0),...
                        rhomax);
                end
            end
        end
        
       %% Main Links
        %for j=1:20  %each link
        %for j=11:20
        for j=1:10
            prev_node=link_node(j,1);
            next_node=link_node(j,2);
            %rho
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
            %v2
            if prev_node==-1
                tmp_v2=0;
            elseif prev_node==0
                tmp_v2=((T/3600)/Ll(j))*V(i,j)*(V(i,j-1)-V(i,j));
            else
                on_link=node_in_out(prev_node,1);
                if ((on_link~=0) && (on_link~=-1))
                    virtual_a=V(i,j-1)*ql1(j-1)...
                        +max(qin2(on_link)*von1(on_link),0);    %Numerator for virtual_v
                    virtual_b=ql1(j-1)+max(qin2(on_link),0);    %Denominator for virtual_v
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
            %v4 on-ramp merging
            if ((prev_node~=0) && (prev_node~=-1))
                ramp_id=node_in_out(prev_node,1);
                if ramp_id~=0
                    tmp_v4=-delta*(T/3600)*qin2(ramp_id)*V(i,j);
                    tmp_v4=tmp_v4/(Ll(j)*lambdal(j)*(rhol1(j)+kappa));
                else
                    tmp_v4=0;
                end
            else
                tmp_v4=0;
            end
            %v5 lane-drop
            if next_node~=-1
                dlambda=max(lambdal(j)-lambdal(j+1),0);
                if dlambda~=0
                    tmp_v5=-phi*(T/3600)*dlambda*rhol1(j)*(V(i,j))^2;
                    tmp_v5=tmp_v5/(Ll(j)*lambdal(j)*rhocrit);
                else
                    tmp_v5=0;
                end
            else
                tmp_v5=0;
            end
            
            V(i+1,j)=min(max(V(i,j)+tmp_v1+tmp_v2+tmp_v3...
                        +tmp_v4+tmp_v5,vmin),120);
           
            %q
            q_cap=tmp_rho_crit*tmp_v_f*exp(-1/tmp_alpha);
            ql2(j)=rhol2(j)*V(i+1,j);
            ql2(j)=min(max(ql(2),0),q_cap)*lambdal(j);

        end
        
       %% Travel Time & Queue Length
        %for k=1:20
        for j=6:9
        %for j=11:14
           pi_tt=pi_tt+rhol2(j)*Ll(j)*lambdal(j); 
        end
        for j=6:7
             pi_w=pi_w+w2(j);
        end
        %pi_w=pi_w+sum(w2);
        %w
        dw=zeros(13,1);
        %for j=1:13
        for j=6:7
            dw(j)=power((max((w2(j)-wmax(j))/wmax(j),0)),2);
        end
        pi=pi+phiw*sum(dw);
        
       %% Discard old state
        %except V - keeping them for travel time estimation
        rhol1=rhol2;
        w1=w2;
        rhooff1=rhooff2;
        ql1=ql2;
        qin1=qin2;
        rhol2=zeros(1,20);
        ql2=zeros(1,20);
        qin2=zeros(1,9);
        w2=zeros(1,13);
        rhooff2=zeros(1,14);
        Qo=zeros(1,13);
    end
    
    %b
    db=zeros(20,Nc-1);
    %for i=1:20
    for i=6:9
       for j=1:Nc-1
          db(i,j)=power((b(i,j+1)-b(i,j)),2); 
       end
    end
    
    %r
    dr=zeros(9,Nc-1);
    %for i=1:9
    for i=3:4
       for j=1:Nc-1
          dr(i,j)=power((r(i,j+1)-r(i,j)),2); 
       end
    end
    
    pi=pi+phir*sum(sum(dr))+phib*sum(sum(db));
    pi=pi+pi_tt+pi_w;
    
end