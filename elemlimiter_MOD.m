%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: Element-based Limiter 
%Type of file: FUNCTION
%Criate date: 11/07/2013 (national workers strik)
%Modify data: 13/05/2015
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: Produce limiter value for 2nd to 6th order in order to guarantee 
%monotonic solution for saturation equation. 

%--------------------------------------------------------------------------
%

%--------------------------------------------------------------------------

function[phi] = elemlimiter_MOD(elemeval,Sw,delta,taylorterms,cvbtype,...
    flagknownvert,satonvertices,constraint)

%Define global parameters:
global coord elem centelem timelevel;

%Define tolerance. It is like a computational zero.
tol = 1e-12;

%Choose Element-Based Limiter according "cvbtype" value
%Obtained from Woodfield et al. (2004)
if strcmp(cvbtype,'wf')
    %Get the amount of elements (face and full). All the elements 
    %surrounding "elemeval"
    [null,esurefull] = getsurelem(elemeval);
    %Attribute to "neighb"
    neighb = esurefull;

    %Get "satelem", "maxsat" and "minsat"
    satelem = Sw(elemeval);
    %Get the bigger and lower saturation value
    maxsat = max(Sw(neighb));
    minsat = min(Sw(neighb));
        
    %Calculate the limiter
    %It avoid division by "0"
    if maxsat - minsat <= 1e-20
        phi = 1;
    %It is safe to division by "0"
    else
        %Define the admintional "gama"
        gama = (satelem - minsat)/(maxsat - minsat);
        %Verify conditions:
        if gama >= 1 || gama <= 0
            phi = 0;
        elseif delta <= gama && gama <= (1 - delta)
            phi = 1;
        elseif 0 < gama && gama < delta
            phi = gama/delta;
        elseif (1 - delta) < gama && gama < 1
            phi = (1 - gama)/delta;
        end  %End of IF (conditions)
    end %End of IF 

%--------------------------------------------------------------------------
%MLP limiter (and variants)

%Obtained from Park and Kim (2012). Look for MLP-u1.
else
    %Get the amount of vertices for the element evaluated
    amountedges = sum(elem(elemeval,1:4) ~= 0);
    %Initialize "qsi"
    qsi = zeros(amountedges,1);
    %Define the vertices
    vertices = elem(elemeval,1:amountedges);
    %Initialize "qsi3rd", "qsi4th", "qsi5th" and "qsi6th" 
    qsi3rd = qsi;
    qsi4th = qsi;
    qsi5th = qsi;
    qsi6th = qsi;

    %Saturation on the centroid
    satoncent = Sw(elemeval);

    %Get some additional parameters if "mlp_vk" is turned on.
    if strcmp(cvbtype,'mlp_vk')
        %Define additional parameters (K, n, etc):
        K = 100;
        K1 = 0.01;
        K2 = 0.01;
        n = 1.5;
        %Initialize "deltax"
        deltax = zeros(amountedges,1);
        %Swept all vertices
        for i = 1:amountedges
            %Put the evaluated vertex as the 1st of the vector "vertices". 
            %It uses "shiftchosen" for that.
            auxverti = shiftchoosen(vertices,vertices(i),'pos');
            %Get "deltax"
            deltax(i) = ...
                norm(coord(auxverti(2),1:2) - coord(auxverti(1),1:2));
        end  %End of FOR
        %Get a mean value for "deltax"
        deltax = mean(deltax);
    end  %End of IF
    
    %Swept all edges
    for i = 1:amountedges
        %Get the amount of elements surrounding each vertex that define the 
        %edge evaluated
        %Define the actual vertex
        actualvtx = vertices(i);
        vertices_aux = shiftchoosen(vertices,i,'pos');
        %Define the next vertex
        nextvtx = vertices_aux(2);
        %Get "esurn" for each vertex:
        [esurn_actual,] = getsurnode(actualvtx);
        [esurn_next,] = getsurnode(nextvtx);
        
        %Verify if the node belongs to "bedge" and if the proble has a
        %known saturation as boundary condition.
        if flagknownvert(actualvtx) == 1 || flagknownvert(nextvtx) == 1 
            %Define saturation on the vertex evaluated
            satonverteval = max(satonvertices([actualvtx nextvtx]));
            %Get the bigger saturatio value among those in the surrounding
            maxsatval = max([satonverteval esurn_actual esurn_next]);
            %Get the lower saturatio value among those in the surrounding
            minsatval = min([satonverteval esurn_actual esurn_next]);
        %The vertex does not belong to known boundary 
        else
            %Get the bigger saturatio value among those in the surrounding
            maxsatval = max(Sw([esurn_actual esurn_next]));
            %Get the lower saturatio value among those in the surrounding
            minsatval = min(Sw([esurn_actual esurn_next]));
        end  %End of IF.
        
        %Verify the accuracy of SOME values.
        %1. Volume average value:
        maxsatval = maxsatval*(abs(maxsatval) > tol);
        minsatval = minsatval*(abs(minsatval) > tol);

        %Get the projection value on the vertex
        vertexcoord = mean(coord([actualvtx nextvtx],1:2),1);
        %Define the distance (centroid -> vertex)
        distcv = vertexcoord - centelem(elemeval,1:2);
        %Chose the type of prjection (second-order or third-order)
        gettaylorsize = size(taylorterms,2);
            
        %------------------------------------------------------------------
        %Evaluate according to "order":
        
%         %------------
%         %FOURTH ORDER (grads, hessian and 3rd deriv - "taylorterms" has 
%         %nine columns)
%         if gettaylorsize >= 9
%             %Initialize "projval_4th"
%             projval_4th = zeros(length(esurn),1);
%             %Swept all elements surrounding the vertex evaluated.
%             for isurnode = 1:length(esurn)
%                 %Define "pointeresurn"
%                 pointeresurn = esurn(isurnode);
%                 
%                 %Get the gradient of saturation for the vertex evaluated:
%                 gradsat = taylorterms(pointeresurn,1:2);
%                 %Get hessian:
%                 hessianterms = taylorterms(pointeresurn,3:5);
%                 %Fill the "hessian" matrix.
%                 hessian = [hessianterms(1) hessianterms(2); ...
%                     hessianterms(2) hessianterms(3)];
%                 %Get third derivative:
%                 der3terms = taylorterms(pointeresurn,6:9);
%             
%                 %Define the terms of the constraint equation. 
%                 %See Clain et al. 2011
%                 consteqterm = constraint(pointeresurn,1:9); 
%         
%                 %Define the distance (centroid -> vertex)
%                 distcv4th = vertexcoord - centelem(pointeresurn,1:2);
% 
%                 %Get "pnprojslope" (see Park and Kim, 2012; ICCFD7 - 4204)
%                 pnprojslope = dot(gradsat,distcv4th) - ...
%                     consteqterm(1:2)*gradsat';
%                 %Get "p1filtered" (see Park and Kim, 2012; ICCFD7 - 4204)
%                 p1filtered = ((distcv4th*hessian)*distcv4th'/2) + ...
%                     (((distcv4th(1)^3)*der3terms(1) + ...
%                     3*(distcv4th(1)^2)*distcv4th(2)*der3terms(2) + ...
%                     3*(distcv4th(2)^2)*distcv4th(1)*der3terms(3) + ...
%                     (distcv4th(2)^3)*der3terms(4))/6) - ...
%                     (consteqterm(3:9)*[[0.5 1 0.5].*hessianterms ...
%                     [(1/6) 0.5 0.5 (1/6)].*der3terms]');
%                 
%                 %Projection of Fourth order on the vertex evaluated.
%                 projval_4th(isurnode) = Sw(pointeresurn) + pnprojslope + ...
%                     p1filtered;
% 
%                 %Attribute "pnprojslope" and "p1filtered" to 
%                 %"pnprojslope_4th" and "p1filtered_4th"
%                 if pointeresurn == elemeval
%                     pnprojslope_4th = pnprojslope;
%                     p1filtered_4th = p1filtered;
%                     projval = projval_4th(isurnode);
%                 end  %End of IF
%             end  %End of FOR (swept "esurn")
%                 
%             %Verify if the node belongs to "bedge" and if the proble has a
%             %known saturation as boundary condition.
%             if flagknownvert(vertices(i)) == 1
%                 %Define saturation on the vertex evaluated
%                 satonverteval = satonvertices(vertices(i));
%                 %Get the bigger saturatio value among those in the sur.
%                 maxprojval = max([satonverteval; projval_4th]);
%                 %Get the lower saturatio value among those in the sur.
%                 minprojval = min([satonverteval; projval_4th]);
%                 %Get "projval" from boundary value
%                 projval = satonvertices(vertices(i));
%             %The vertex does not belong to known boundary 
%             else
%                 %Get the bigger saturatio value among those in the sur.
%                 maxprojval = max(projval_4th);
%                 %Get the lower saturatio value among those in the sur.
%                 minprojval = min(projval_4th);
%             end  %End of IF.
%             
%             %Evaluate the parameter "qsi4th"
%             %1. Projection value:
%             maxprojval = maxprojval*(abs(maxprojval) > tol);
%             minprojval = minprojval*(abs(minprojval) > tol);
%             %Verify the accuracy of all values.
%             %2. The projected value (4th order) value:
%             projval = projval*(abs(projval) > tol);
%             %3. 2nd order projection ("pnprojslope_4th"):
%             pnprojslope_4th = pnprojslope_4th*(abs(pnprojslope_4th) > tol);
%             %4. very higher order projection ("pnprojslope_4th"):
%             p1filtered_4th = p1filtered_4th*(abs(p1filtered_4th) > tol);
% 
%             %Define augmented MLP criteria:
%             qsi4th(i) = ...
%                 ((projval >= minprojval && minprojval >= minsatval) && ...
%                 (projval <= maxprojval && maxprojval <= maxsatval));
%             
%             %Verify other condition
% %             if qsi4th(i) == 0 && (projval - minsatval > 0) 
% %                 qsi4th(i) = (pnprojslope_4th > 0 && p1filtered_4th < 0) || ...
% %                     (abs((projval - satoncent)/(satoncent + 1e-16)) <= ...
% %                     1e-3);
% %             elseif qsi4th(i) == 0 && (projval - maxsatval < 0)
% %                 qsi4th(i) = (pnprojslope_4th < 0 && p1filtered_4th > 0) || ...
% %                     (abs((projval - satoncent)/(satoncent + 1e-16)) <= ...
% %                     1e-3);
% %             end  %End of IF (aditional condition)
% 
%             if qsi4th(i) == 0
%                 qsi4th(i) = ...
%                     ((sign(pnprojslope_4th) ~= sign(p1filtered_4th) && ...
%                     (projval > minsatval && projval < maxsatval)) || ...
%                     (abs((projval - satoncent)/(satoncent + 1e-16)) <= ...
%                     1e-3));
%             end  %End of IF
%             
%             %Attribute the limiter value of 4th order to 3rd order limiter.
%             %qsi3rd(i) = qsi4th(i);
%         end  %End of IF (fourth order)
%         
%         %-----------
%         %THIRD ORDER (grads and hessian - "taylorterms" has five columns)
%         if gettaylorsize >= 5 %&& qsi3rd(i) == 0
%             %Initialize "projvalvho"
%             projvalvho = zeros(length(esurn),1);
%             %Swept all elements surrounding the vertex evaluated.
%             for isurnode = 1:length(esurn)
%                 %Define "pointeresurn"
%                 pointeresurn = esurn(isurnode);
% 
%                 %Get the gradient of saturation:
%                 gradsat = taylorterms(pointeresurn,1:2);
%                 %Get hessian for each elament of "esurn"
%                 %Define the Hessian matrix (3rd to 5th "taylorterms" 
%                 %columns):
%                 hessianterms = taylorterms(pointeresurn,3:5);
%                 %Fill the "hessian" matrix.
%                 hessian = [hessianterms(1) hessianterms(2); ...
%                     hessianterms(2) hessianterms(3)];
%                 %Define the terms of the constraint equation. 
%                 %See Clain et al. 2011 (Eq. 11).
%                 consteqterm = constraint(pointeresurn,1:5); 
% 
%                 %Define the distance (centroid -> vertex)
%                 distcvho = vertexcoord - centelem(pointeresurn,1:2);
% 
%                 %Get "pnprojslope" (see Park and Kim, 2012; ICCFD7 - 4204)
%                 pnprojslope = dot(gradsat,distcvho) - ...
%                     consteqterm(1:2)*gradsat';  
%                 %Get "p1filtered" (see Park and Kim, 2012; ICCFD7 - 4204)
%                 p1filtered = (((distcvho*hessian)*distcvho')/2) - ...
%                     (consteqterm(3:5)*([0.5 1 0.5].*hessianterms)');
%                 
%                 %Projection of third order on the vertex evaluated.
%                 projvalvho(isurnode) = Sw(pointeresurn) + ...
%                     pnprojslope + p1filtered; 
% 
% %                     (dot(gradsat,distcvho) + ...
% %                     ((distcvho*hessian)*distcvho')/2) - ...
% %                     (consteqterm*[gradsat [0.5 1 0.5].*hessianterms]');
% 
%                 %Attribute "pnprojslope" and "p1filtered" to 
%                 %"pnprojslope_3rd" and "p1filtered_3rd"
%                 if pointeresurn == elemeval
%                     pnprojslope_3rd = pnprojslope;
%                     p1filtered_3rd = p1filtered;
%                 end  %End of IF
%             end  %End of FOR ("esurn")
% 
%             %Get the projection for the "elemeval"
%             projval = projvalvho(logical(esurn == elemeval));
%             
%             %--------------------------------------------------------------
%             %Get the maximum and minimum projections on the evaluated
%             %vertex:
%             
%             %Verify if the node belongs to "bedge" and if the proble has a
%             %known saturation as boundary condition.
%             if flagknownvert(vertices(i)) == 1
%                 %Define saturation on the vertex evaluated
%                 satonverteval = satonvertices(vertices(i));
%                 %Get the bigger saturatio value among those in the sur.
%                 maxprojval = max([satonverteval; projvalvho]);
%                 %Get the lower saturatio value among those in the sur.
%                 minprojval = min([satonverteval; projvalvho]);
%                 %Get "projval" from boundary value
%                 projval = satonvertices(vertices(i));
%             %The vertex does not belong to known boundary 
%             else
%                 %Get the bigger saturatio value among those in the sur.
%                 maxprojval = max(projvalvho);
%                 %Get the lower saturatio value among those in the sur.
%                 minprojval = min(projvalvho);
%             end  %End of IF.
%             
%             %Evaluate the parameter "qsi3rd"
%             %1. Projection value:
%             maxprojval = maxprojval*(abs(maxprojval) > tol);
%             minprojval = minprojval*(abs(minprojval) > tol);
%             %2. Projected value on the vertex evaluated (from the element 
%             %evaluated)
%             projval = projval*(abs(projval) > tol);
%             %3. 2nd order projection ("pnprojslope_4th"):
%             pnprojslope_3rd = pnprojslope_3rd*(abs(pnprojslope_3rd) > tol);
%             %4. very higher order projection ("pnprojslope_4th"):
%             p1filtered_3rd = p1filtered_3rd*(abs(p1filtered_3rd) > tol);
%             
%             %Define augmented MLP criteria:
%             qsi3rd(i) = ...
%                 ((projval >= minprojval && minprojval >= minsatval) && ...
%                 (projval <= maxprojval && maxprojval <= maxsatval));
% 
%             %Verify other condition
% %             if qsi3rd(i) == 0 && (projval - minsatval > 0) 
% %                 qsi3rd(i) = (pnprojslope_3rd > 0 && p1filtered_3rd < 0) || ...
% %                     (abs((projval - satoncent)/(satoncent + 1e-16)) <= ...
% %                     1e-3);
% %             elseif qsi3rd(i) == 0 && (projval - maxsatval < 0)
% %                 qsi3rd(i) = (pnprojslope_3rd < 0 && p1filtered_3rd > 0) || ...
% %                     (abs((projval - satoncent)/(satoncent + 1e-16)) <= ...
% %                     1e-3);
% %             end  %End of IF (aditional condition)
%             %Verify other condition
% %             if qsi3rd(i) == 0 && gettaylorsize > 5 
% %                 qsi3rd(i) = ...
% %                     ((sign(pnprojslope_3rd) ~= sign(p1filtered_3rd) && ...
% %                     (projval > minsatval && projval < maxsatval)) || ...
% %                     (abs((projval - satoncent)/(satoncent + 1e-16)) <= ...
% %                     1e-3));
% %             end  %End of IF
%         end  %End of IF (third order)
    
        %------------
        %SECOND ORDER (only grads - "grads" has two columns)
        %Verify if the saturation on the vertex evaluated is known
        if flagknownvert(actualvtx) == 1 || flagknownvert(nextvtx) == 1
            projval = max(satonvertices([actualvtx nextvtx])) - satoncent;
        %The vertex does not belong to known boundary 
        else
            %Projection on the vertex
            projval = dot(taylorterms(elemeval,1:2),distcv);
        end  %End of IF
        
        %Verify if there is a null projection and the very higher order is
        %turned off.
        if (abs(projval) > tol && gettaylorsize <= 5) || ...
                (gettaylorsize >= 9 && qsi3rd(i) == 0 && ...
                abs(projval) > tol)
            %Get the ratio "rv":
            %Define a boolean condition
            booleanrv = projval > 0;
            %Defien "deltasatval"

            deltasatval = (maxsatval - satoncent)*booleanrv + ...
                (1 - booleanrv)*(minsatval - satoncent);
            %Get the ratio "rv":
            rv = deltasatval/projval;
            
            if rv < 1
                deltasatval = (max(min(maxsatval,0.665),satoncent) - satoncent)*booleanrv + ...
                    (1 - booleanrv)*(minsatval - satoncent);
                rv = deltasatval/projval;
            end
            
            %Calculate "qsi" (It Choces according flag "cvbtype")
            switch cvbtype
                %Original MLP (Park and Kin, 2012)
                case 'mlp'
                    qsi(i) = min(1,rv);
                %MLP u2, Venkatakrishnan (see Park and Kin, 2012)
                case 'mlp_vk'
                    %Define "deltaplus" according "projval" sign.
                    booleanrv = projval > 0;
                    %Defien "deltaplus"
                    deltaplus = (maxsatval - satoncent)*booleanrv + ...
                        (1 - booleanrv)*(minsatval - satoncent);
                    
                    %Define "deltaless"
                    deltaless = projval;
                    %Delta sw (difference between the max and min 
                    %saturation values):
                    deltasonvertex = maxsatval - minsatval;
                    %Calculate "teta"
                    teta = deltasonvertex/(K2*(deltax^n));
                    %Calculate epsilon^2:
                    epsilonsqre = (K1*(deltasonvertex^2))/(1 + teta);  %K*(deltax^3);%    
                    
                    rv = deltaplus/(deltaless + 1e-16);
                    %Use Venkathakrishnen approximation
%                     qsi(i) = (((deltaplus^2) + epsilonsqre) + ...
%                         2*deltaless*deltaplus)/...
%                         ((deltaplus^2) + 2*(deltaless^2) + ...
%                         deltaless*deltaplus + epsilonsqre);
                    qsi(i) = (rv^3 + 3*rv)/(rv^3 + (rv^2) + 4);%(rv^2 + 2*rv)/(rv^2 + rv + 3);
                %MLP-Van Albada
                case 'mlp_va'
                    qsi(i) = (rv^2 + rv)/(rv^2 + 1);
                %MLP-Van Albada 2 (Löhner, 2001)
                case 'mlp_va2'
                    qsi(i) = (2*rv)/(rv^2 + 1);
                %MLP-Van Leer 1 (TVD second order for any value of "rv")
                case 'mlp_vl'
                    qsi = 2*rv/(rv + 1);
                %MLP-Van Leer 2 (first-order for "rv" bigger than 1)
                case 'mlp_vl2'
                    qsi(i) = 4*rv/((rv + 1)^2);
                %MLP-New Limiter 2 (Mandal, 2008)
                case 'mlp_nl2'
                    qsi(i) = (rv^3 + rv)/(rv^3 + 1);
                %MLP-New Limiter 2 (Mandal, 2008)
                case 'mlp_nl3'
                    qsi(i) = (rv^2 + 3*rv)/(rv^2 + rv + 2);
                case 'mlp_tst'
                    qsi(i) = 0.5*(rv + 1)*...
                        min((4*rv*(3*rv + 1)/(11*(rv^2) + 4*rv + 1)),...
                        (4*(rv + 3)/((rv^2) + 4*rv + 11)));
            end  %End of SWITCH
        
        %There is a null projection
        else
            qsi(i) = 1;
        end  %End of IF
    end  %End of FOR (swept the vertices)

    %Define the limiter ("phi")
    %Choose according to size of "taylorterms"
    %Second order:
    if gettaylorsize == 2
        phi = min(qsi);
    %Third order limited
    elseif gettaylorsize == 5
        phi = [min(qsi) min(qsi3rd)];
    %Fourth order limited
    elseif gettaylorsize == 9
        phi = [min(qsi) min(qsi3rd) min(qsi4th)];        
    %Fifth order limited
    elseif gettaylorsize == 14
        phi = [min(qsi) min(qsi3rd) min(qsi4th) min(qsi5th)];        
    %Sixth order UNlimited
    elseif gettaylorsize == 20 
        phi = [min(qsi) min(qsi3rd) min(qsi4th) min(qsi5th) min(qsi6th)]; 
    end  %End of IF (MLP limiter)
end  %End of IF (limiter choice)



% line = [4 12 20 28 36 44 52 60 68 76]; 
%     %84 92 100 108 116 124 132 140 148 156];
%     if timelevel == 132 && ismember(elemeval,line)
%         elemeval
%         phi
%     end