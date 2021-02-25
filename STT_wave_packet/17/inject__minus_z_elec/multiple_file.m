clear;

F = ones(5);
for ii=1:5
     file_name=['file' sprintf('%d',ii) '.txt'];
     fileID=fopen(file_name, 'w+');
     fprintf(fileID,'%2.0f',F);
     fclose(fileID);
  end
