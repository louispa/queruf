function [ M ] = reader( f )
%Imprime les résultats du labo 2 et les compare aux modèles théoriques

%Résultats expérimentaux
fileid=fopen(f);
[A]=textscan(fileid,'%f %f %f %f %f %f');
fclose(fileid);

M=[A{1} A{2} A{3} A{4} A{5} A{6}];

end

