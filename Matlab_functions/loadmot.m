function res = loadmot(file)

[path,filename,extension] = fileparts(file);
cd(path);
mot2txt = strrep([filename,'.mot'],'.mot','.txt');
copyfile([filename,'.mot'],mot2txt);
res = readmatrix([filename,'.txt']);
delete([filename,'.txt'])
cd ..\..\..\

end