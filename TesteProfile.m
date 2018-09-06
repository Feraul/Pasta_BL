% profile on -history;
%------------------
%Comando "ismember"

% a = 1000*rand(1,1300);
% for i = 1:22500
% %     b = ismember(i,a);
%     b = sum(i == a);
% end

%------------------
%Comando "setdiff"

% a = [2 4 5 6 9 4 3 5 6 7 34 65 9 5 67 4 3 12 1 3];
% for i = 1:22500
% %     b = setdiff(3,a);
%     b = a(logical(a ~= 3));
% end



[coord,centelem,elem,esurn1,esurn2,nsurn1,nsurn2,bedge,inedge,normals,...
    esureface1,esureface2,esurefull1,esurefull2,elemarea,dens,visc,...
    satlimit,pormap,bcflag,courant,totaltime,numcase,phasekey,pmethod,...
    smethod,xyrz,r0,symaxe,keymsfv,coarseratio,auxcvfactor,interptype,...
    nonlinparam,multdopt,goefreeopt,order,timeorder,recovtype,lsneightype,...
    lsexp,keygravity,g,keycapil,ncaplcorey,filepath,resfolder,kmap,wells,...
    limiterflag] = preprocessor;

tic
% profile off
% profsave(profile('info'),'teste02')
for i = 1:size(elem,1)
%         g = all(ismember(bedge(:,1:2),...
%             [elem(i,1) elem(i,2)]),2)  ;                     
          g = all(elem(i,1) == bedge(:,1:2) | elem(i,2) == bedge(:,1:2),2);
end
toc
pause

        i = all(elem(ielem,inode) == bedge(:,1:2) | elem(ielem,next) == ...
            bedge(:,1:2),2);
