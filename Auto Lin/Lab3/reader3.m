function [ M ] = reader3( f )
%Imprime les r�sultats du labo 3 et les compare aux mod�les th�oriques

%R�sultats exp�rimentaux
fileid=fopen(f);
[A]=textscan(fileid,'%f %f %f %f','HeaderLines',1)
fclose(fileid);

M=[A{1} A{2} A{3} A{4}]

end

