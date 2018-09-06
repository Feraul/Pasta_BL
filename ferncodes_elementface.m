%It is called by "preMPFA.m"

function [V,N,F] = ferncodes_elementface(nflag)
global bedge inedge elem coord esurn1 esurn2 nsurn1 nsurn2

for ii=1:size(elem,1)
    i=1;
    n=elem(ii,4);
    if n~=0
        list=elem(ii,1:4);
    else
        list=elem(ii,1:3);
    end
    for jj=1:length(list)
        if jj<length(list)
            ibedge=find(((bedge(:,1)==list(jj+1) & bedge(:,2)==list(jj))|(bedge(:,1)==list(jj) & bedge(:,2)==list(jj+1))));
            
            if ibedge~=0
                F(ii,i)=ibedge + size(inedge,1); % flag da face 
                i=i+1;
            else
                iedge=find((inedge(:,1)==list(jj+1) & inedge(:,2)==list(jj))|(inedge(:,1)==list(jj) & inedge(:,2)==list(jj+1)));
                F(ii,i)=iedge;
                i=i+1;
            end
        else
            ibedge=find(((bedge(:,1)==list(jj) & bedge(:,2)==list(1))|(bedge(:,1)==list(1) & bedge(:,2)==list(jj))));
            
            if ibedge~=0
                F(ii,i)=ibedge + size(inedge,1);
                i=i+1;
            else
                iedge=find((inedge(:,1)==list(jj) & inedge(:,2)==list(1))|(inedge(:,1)==list(1) & inedge(:,2)==list(jj)));
                F(ii,i)=iedge;
                
                i=i+1;
            end
        end
    end
    
    
end
for No=1:size(coord,1)
    N_element_No=esurn2(No+1)-esurn2(No);
    
    n_pontos=nsurn2(No+1)-nsurn2(No);
    
    %construção do tensor permeabilidade.%
    
    for k=1:N_element_No
        
        n1=elem(esurn1(esurn2(No)+k),1);
        n2=elem(esurn1(esurn2(No)+k),2);
        n3=elem(esurn1(esurn2(No)+k),3);
        n4=elem(esurn1(esurn2(No)+k),4);
        a=zeros(2,1);
        ii=1;
        for jj=[n1,n2,n3,n4]
            if jj~=No && jj==0
                a(ii,1)=jj;
                ii=ii+1;
            elseif jj~=No && jj~=0
                for g=1:n_pontos
                    h=nsurn1(nsurn2(No)+g);
                    if jj==h
                        a(ii,1)=jj;
                        ii=ii+1;
                    end
                end
            end
        end
        list=nsurn1(nsurn2(No)+1:nsurn2(No+1));
        list2=esurn1(esurn2(No)+1:esurn2(No+1));
        
        for g=1:size(list,1)
            h=list(g);
            if length(list)==length(list2)
                if length(list)==k
                    if a(1,1)==list(g)
                        r=size(list,1)+1;
                    elseif a(2,1)==h
                        s=length(list);
                    end
                else
                    if a(1,1)==h
                        r=g;
                    elseif a(2,1)==h
                        s=g;
                    end
                end
            else
                if a(1,1)==h
                    r=g;
                elseif a(2,1)==h
                    s=g;
                end
            end
        end
        
        %Verify if the vertex belongs to "begde"
%         pointboundvtx = any(No == bedge(:,1:2),2);
%         if any(pointboundvtx) == 0
        if nflag(No,1) == 50000
            
            if r==n_pontos & s==n_pontos
                
                
                m=a(2,1);
                n=a(1,1);
                
            elseif r<s
                m=a(1,1);
                n=a(2,1);
            else
                m=a(2,1);
                n=a(1,1);
            end
            
            
        else
            if r<s
                m=a(1,1);
                n=a(2,1);
            else
                m=a(2,1);
                n=a(1,1);
            end
        end
        
        
        [Tt]=faces_no(bedge,inedge,No,m,n);
        V(:,k,No)= Tt';
        
        
    end
    
    m1=1;
    vetor1=nsurn1(nsurn2(No)+1:nsurn2(No+1));
    for j= [vetor1']
        ibedge=find(((bedge(:,1)==No & bedge(:,2)==j)|(bedge(:,1)==j & bedge(:,2)==No)));
        if ibedge~=0
            
            N(No,m1)=ibedge+size(inedge,1);
            m1=m1+1;
        else
            iedge=find(((inedge(:,1)==j & inedge(:,2)==No)|(inedge(:,1)==No & inedge(:,2)==j)));
            
            N(No,m1)=iedge;
            m1=m1+1;
        end
    end

end

end