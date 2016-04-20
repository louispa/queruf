function [ M ] = reader( f )
%Imprime les résultats du labo 1 et les compare aux modèles théoriques

%Résultats expérimentaux
fileid=fopen(f);
[A]=textscan(fileid,'%f %f %f %f %f','HeaderLines',1);
fclose(fileid);

M=[A{1} A{2} A{3} A{4} A{5}];

end

