function  D2export(Re,X,Y,NGX,NGY,varName1,varargin)
global Datpath
result = [X Y];
varName2 = 'VARIABLES=x,y,';
varName = [varName2,varName1];
formatSpec ='%2.15f %2.15f\n';
a = '%2.15f ';

for i = 1:size(varargin,2)
var1 = varargin{i};
    result = [result var1];
    formatSpec = [a,formatSpec];
end
if issparse(result)
    result = full(result);
end
fname = [char(strcat('../',Datpath,'/')),'D2Date_',num2str(Re,'%i'),'.dat'];   %指定文件名
fid = fopen(fname,'wt');
fprintf(fid,varName);
fprintf(fid,'ZONE I=%g J=%g F=POINT\n',NGX,NGY);
fprintf(fid,formatSpec,result');
fclose(fid);
% save(char(strcat('../',Matpath,'/','D2Data_Multy')))