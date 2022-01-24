function inputsfeap(NombreFichero,ValorVariable)
% Función que crea ficheros input de FEAP.


%% Creación de los ficheros de entrada

fid     = fopen(NombreFichero,'w');
Linea   = ['PARAmeters\n','Lz = ',num2str(ValorVariable),'\n\n'];


fprintf(fid,'FEAP\n');
fprintf(fid,'0 0 0 3 6 27\n\n');
fprintf(fid,'NOPRint\n\n');
fprintf(fid,'INCLude,ipara\n\n');
fprintf(fid,'INCLude,iparpul\n\n');
fprintf(fid,'INCLude,icorr\n\n');
fprintf(fid,Linea);
fprintf(fid,'INCLude,iinte\n\n');
fprintf(fid,'INCLude,imesh\n\n');
fprintf(fid,'INCLude,imate\n\n');
fprintf(fid,'INCLude,iboun\n\n');
fprintf(fid,'END\n\n');
fprintf(fid,'TIE,1.e-16\n\n');
fprintf(fid,'OPTImize\n\n');
fprintf(fid,'batch\n \t tplo\nend\ndisp,501,4\ndisp,501,5\ndisp,nf,5\n\n');
fprintf(fid,'INCLude,iresol3\n\n');
fprintf(fid,'STOP\n');

