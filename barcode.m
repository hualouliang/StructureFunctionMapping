function lambda = barcode(C)
% 
% Compute compute the barcode of its filtration of connectivity matrix via Kruskal's algorithm for 
% finding a minimum spanning tree (MST) 
%
% lambda = barcode(C); 

%
% Hualou Liang at Drexel University, 2015
%

G=sparse(C);
[ST, pred] = graphminspantree(G, 'Method', 'Kruskal');
[I,J,V] = find(ST);
lambda=sort(V); 

