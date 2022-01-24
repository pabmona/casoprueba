function crearipulse(PuntosPulso,Dtpd,dtc)

fid=fopen('iparpul','w');
fprintf(fid,'PARAmeter \n');
fprintf(fid,'dc = %e \n',dtc);
fprintf(fid,'fc = %e \n',dtc+(Dtpd*1.8));
fprintf(fid,'ft = %e \n',dtc+Dtpd);
fprintf(fid,'it = %e \n',dtc);
fprintf(fid,'f2 = %e \n',dtc+(Dtpd*10));
fprintf(fid,'\n');
fclose(fid);

m   = length(PuntosPulso);
pro = interp1(linspace(0,Dtpd,m),PuntosPulso,(0:dtc:Dtpd)');

save('ipulse','pro','-ascii')
