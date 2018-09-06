%--------------------------------------------------------------------------
%Define the element:

%Face 1, vertices 659 - 205 (horizontal inferior):
%numb = 2470;    %half-edge far (point 1)
%numb = 653;    %half-edge near (point 1)

%Face 2, vertices 205 - 661 (vertical right):
%numb = 656;    %half-edge far (point 1)
%numb = 2478;    %half-edge near (point 1)

%Face 3, vertices 661 - 663 (horizontal superior):
%numb = 2477;    %half-edge far (point 1)
%numb = 2487;    %half-edge near (point 1)

%Face 3, vertices 663 - 659 (horizontal superior):
%numb = 2486;    %half-edge far (point 1)
numb = 2471;    %half-edge near (point 1)


%--------------------------------------------------------------------------
%Data Raio Menor:

%point 1, vertex 181 (diagonal):
%numb = 557;    %half-edge far (point 1)
%numb = 560;    %half-edge near (point 1)

%point 2, vertex 589 (intermed.):
%numb = 2189;  %half-edge far (point 2)
%numb = 2192;  %half-edge near (point 2)

%point 3, vertex 529 (paralel.):
%numb = 1952;  %half-edge far (point 3)
%numb = 1951;  %half-edge near (point 3)

%--------------------------------------------------------------------------
%Data Raio Medio:

%point 1, vertex 221 (diagonal):
%numb = 717;    %half-edge far (point 1)
%numb = 720;    %half-edge near (point 1)

%point 2, vertex 653 (interm.):
%numb = 2448;  %half-edge far (point 2)
%numb = 2447;  %half-edge near (point 2)

%point 3, vertex 539 (paralel.):
%numb = 1989;  %half-edge far (point 3)
%numb = 1992;  %half-edge near (point 3)

%--------------------------------------------------------------------------
%Data Raio Maior:

%point 1, vertex 834 (diagonal):
%numb = 3169;    %half-edge far (point 1)
%numb = 3172;    %half-edge near (point 1)

%point 2, vertex 663 (intermed.):
%numb = 2485;  %half-edge far (point 2)
%numb = 2488;  %half-edge near (point 2)

%point 3, vertex 545 (paralel.):
%numb = 2013;  %half-edge far (point 3)
%numb = 2016;  %half-edge near (point 3)

%--------------------------------------------------------------------------
vpi = 0.2;
amountime = 1001;
mobratio = 'M1';
filepath = 'C:\\Users\\Marcio\\Doutorado\\Programas';
foldername = 'BenchTwophase4_3';
propertie = 'flowrate';

%fname = sprintf('%s\\%s\\@weights_M1\\',char(filepath),char(foldername));
fname = sprintf('%s\\%s\\@flowrate_M1\\',char(filepath),char(foldername));
%fname = sprintf('%s\\%s\\@fw_M1\\',char(filepath),char(foldername));

for i = 1:amountime
    %Time level whose the results will appear
    step = num2str(i);
%    name = [fname 'weights_' mobratio ' ' step '.dat'];
    name = [fname 'flowrate_' mobratio ' ' step '.dat'];
%    name = [fname 'fw_' mobratio ' ' step '.dat'];
    readfile = fopen(name);

    getdata = textscan(readfile,'%f',1,'HeaderLines',numb - 1);
    vector(i) = getdata{1};
    storevpi(i) = i*vpi/amountime;
    fclose(readfile);
end

result = horzcat(storevpi',vector');

strnumb = num2str(numb);
fname = sprintf('%s\\%s\\',char(filepath),char(foldername));
namefile = [fname 'evaluation_' propertie '_' mobratio '_' strnumb '.dat'];
write = fopen(namefile,'w');
for i = 1:amountime
    fprintf(write,'%26.16E %26.16E\r\n',result(i,:));
end

fclose(write);

% fprintf(geometry,'%s\r\n',char(head(i))');
% fprintf(write,'%12.10f\r\n',storeweights);
% 
% fclose(write);
