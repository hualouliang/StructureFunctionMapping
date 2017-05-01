% Demo codes for structure-function topological mapping and persistent homolgy
% 
% Hualou Liang at Drexel University, 2015
%


% load structural and functional matrices (66 x 66)
load data_S_F

% run the topological mapping, use binary structural matrix
[SSEs, corvals, Betti0s,SSE_bs] = topo_mapping(F,S~=0);

% plot SSE based on Betti curves
plot(SSE_bs,'--o')
xlabel('Maximum path length');
ylabel('SSE_\beta')

% % plot regular SSE based on Frobenius norm 
figure
plot(SSEs,'--o')
xlabel('Maximum path length');
ylabel('SSE')

% plot correlation between target and estimated funct matrices
figure
plot(corvals,'--o')
xlabel('Maximum path length');
ylabel('Correlation between actual and predicted FC')

% plot the barcodes
figure
N = size(F, 1);
for i=1:10,
    % plot([0; Betti0s(:, i)],[N:-1:1],'Color',color(i,:),'LineWidth',2);
    plot([0; Betti0s(:, i)],[N:-1:1],'LineWidth',2);
    hold on
end
dF=barcode(1-F); % target
plot([0 dF'],[N:-1:1], 'Color','r', 'LineWidth',2);
xlabel('Filtration value'); 
ylabel('Number of connected components'); % Betti number - beta0
legend('F1', 'F2','F3','F4','F5','F6','F7','F8','F9','F10','F')
