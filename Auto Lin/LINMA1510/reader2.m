function [ M ] = reader2( f )
%Imprime les résultats du labo 2 et les compare aux modèles théoriques
%L'en tête est de la forme:
%blabla || blabla || blabla || blabla || blabla ||
%blabla
%(2 lignes)
%Résultats expérimentaux
fileid=fopen(f);
[A]=textscan(fileid,'%f %f %f %f %f %f','HeaderLines',2);
fclose(fileid);

M=[A{1} A{2} A{3} A{4} A{5} A{6}];

end

