function VAZAOAGUA_writefile(numb,mobratio,propertie)
%--------------------------------------------------------------------------
%Define the element:

%--------------------------------------------------------------------------
vpi = 0.2;
amountime = 1001;
filepath = 'C:\\Users\\Marcio\\Doutorado\\Programas';
foldername = 'BenchTwophase4_3';

fname = sprintf('%s\\%s\\@%s_%s\\',char(filepath),char(foldername),char(propertie),char(mobratio));

for i = 1:amountime
    %Time level whose the results will appear
    step = num2str(i);
    name = [fname propertie '_' mobratio ' ' step '.dat'];
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
