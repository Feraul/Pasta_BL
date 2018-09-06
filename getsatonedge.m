%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 12/01/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%   

%--------------------------------------------------------------------------
%Additional comments: 
%"taylorterms" is acctualy the grad of saturation (2nd order applications). 
%When we use a third-order scheme, we do not have only a grad. We have also 
%a hessian matrix. Thus, for a more generalized understanding we named the 
%matrix like that.

%--------------------------------------------------------------------------

function [swonedge] = getsatonedge(elemeval,vertices,verticescoord,...
    taylorterms,Sw,limiterflag,order,constraint,flagknownvert,...
    satonvertices,mlplimiter)
%Define the strategy recovery according the "smethod" flag
%First-Order Upwind 
if order == 1
    %Recovery "Sw" by a simply strategy Upwind.
    swonedge = Sw(elemeval(1));
    
%Higher-Order Schemes (Second Order).
elseif order == 2
    %Calculate the second-order saturation value.
    swonedge = get2ndorder(elemeval,vertices,verticescoord,...
        taylorterms(:,1:2),Sw,limiterflag,flagknownvert,satonvertices,...
        mlplimiter);

end  %End of IF ("order")

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION:
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Function "get2ndorder"
%--------------------------------------------------------------------------

function [swonedge] = get2ndorder(elemeval,vertices,verticescoord,...
    taylorterms,Sw,limiterflag,flagknownvert,satonvertices,mlplimiter)
%Define global parameters:
global centelem recovtype;

%Initialize parameter "k" (It is "-1" for Upwind approximation)  
%and "epsilon" is the machine limit (1e-16)
k = limiterflag{7};
epsilon = 1e-16;
%Initialize points coordinate:
pointl = centelem(elemeval(1),1:2);
%element inside the domain
if length(elemeval) > 1
    pointr = centelem(elemeval(2),1:2);
%element on the boundary
else
    pointr = pointl;
end  %End of IF (elements)

%Initialize "phi" and "grad". Each row is a element gradient.
phi = ones(length(elemeval),1);
grad = zeros(length(elemeval),2);

%--------------------------------------------------------------------------
%Verify if there is an Element-Based Limiter and Which one it is.
        
%Get the type of Element-based limiter
cvbtype = limiterflag{4};
%Get the tune parameter
delta = limiterflag{5};

%There is BOTH Edge-Based Limiter and a Element-Based Limiter
if strcmp(limiterflag{1},'on') && strcmp(limiterflag{3},'on') && ...
        strcmp(recovtype,'ggg') == 0
    %Define Element-Based Limiter
    %Swept the amount of elements shared by the edge evaluated
    for ishare = 1:length(elemeval)
        %Get Element-based limiter
        phi(ishare) = elemlimiter(elemeval(ishare),Sw,delta,flagknownvert,satonvertices);
        %Get the gradient (reconstructed by MUSCL strategy) for each 
        %element evaluated (in each matrix row)
        grad(ishare,:) = ...
            phi(ishare)*taylorterms(elemeval(ishare),:);
    end  %End of FOR

%There is ONLY an Edge-Based Limiter
elseif strcmp(limiterflag{1},'on') && strcmp(limiterflag{3},'off')
    %Define Element-Based Limiter
    %Swept the amount of elements shared by the edge evaluated
    for ishare = 1:length(elemeval)
        %Get the gradient (reconstructed by MUSCL strategy) for each 
        %element evaluated (in each matrix row)
        grad(ishare,:) = taylorterms(elemeval(ishare),:);
    end  %End of FOR

end  %End of IF
        
%--------------------------------------------------------------------------
%Verify if there is an Edge-Based Limiter and Which one it is.
        
%Get the type of Edge-based limiter
eblimtype = limiterflag{2};

%There is a Edge-Based Limiter and the edge is inside the domain.
if strcmp(limiterflag{1},'on') && length(elemeval) > 1
    %Calculate the gradients centered and upwind.
    %Calculate the vector distance between colocation points.
    dlr = pointr - pointl; 

    %Calculate the gradient centered ("gradcent")
    gradcent = Sw(elemeval(2)) - Sw(elemeval(1));  
    %Calculate the gradient upwind ("gradupwd"). First "grad" matrix row. 
    gradupwd = 2*dot(grad(1,:),dlr) - gradcent; 
            
    %Define "r" (smoothness factor). "epsilon" avoids division by 0
    r = gradcent/(gradupwd + epsilon);
            
    %Define Edge-Based Limiter
    qsi = max(0,edgelimiter(r,eblimtype));
    %Calculate the gradient limited according "k" value.
    gradslimited = qsi*0.5*((1 - k)*gradupwd + (1 + k)*gradcent);

    %Recovery "Sw" on the edge just by using the Edge-Based Limiter
    swonedge = Sw(elemeval(1)) + 0.5*gradslimited;

%There is a Edge-Based Limiter and the edge is over the boundary.
elseif strcmp(limiterflag{1},'on') && length(elemeval) == 1
    %Put First order in recovery process. 
    swonedge = Sw(elemeval(1));
end  %End of IF
        

            


