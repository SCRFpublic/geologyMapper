
%% Compile Mex files

of = cd ('toolbox_graph');
compile_mex;
cd(of);
of = cd('toolbox_graph');
compile_mex;
cd(of);

%% Binary geology example

clear;
load BWGanges.mat;

nNodes = 50;
plotRes = 1;
topoMapper(BW(1:2:end,1:2:end), nNodes, plotRes);


%% 2D heterogeneous porous media example

clear;
load spe10data.mat

layer = 10;
nNodes = 20;
plotRes = 1;

% log K
logK2D  = permFieldUpperNess(:,1:200,layer);
% milliDarcy
K2D     = exp(logK2D);

figure

subplot(2,1,1)
imagesc(logK2D);
title('Log K')
colorbar
set(gca, 'FontSize', 18)

subplot(2,1,2)
imagesc(K2D);
title('K [mD]')
colorbar
set(gca, 'FontSize', 18)

G2D = permTopoMapper(K2D, nNodes, plotRes);


%% Compute distances between multiple graph realizations of two 3D cases

% indices
isel    = 11:50;
jsel    = 51:150;
ksel    = 1:20;
% number of nodes
nNodes  = 30;

% Two 3D cubes
K1 = exp(permFieldTarbert(isel, jsel, ksel));
K2 = exp(permFieldUpperNess(isel, jsel, ksel));


% examples for each case
permTopoMapper(K1, nNodes, 1);
permTopoMapper(K2, nNodes, 1);


% multiple realizations of the same 3D cubes
nReal = 10;

laplacianSpectrumReducedGraphs1 = zeros(nReal, nNodes);
laplacianSpectrumFullGraphs1 = zeros(nReal, nNodes*40);
laplacianSpectrumReducedGraphs2 = zeros(nReal, nNodes);
laplacianSpectrumFullGraphs2 = zeros(nReal, nNodes*40);

for i = 1:nReal

    % map graphs
    G1 =  permTopoMapper(K1, nNodes, 0);
    G2 =  permTopoMapper(K2, nNodes, 0);
    
    % compute eigenvalues
    laplacianSpectrumReducedGraphs1(i,:) = computeLaplacianEV(G1.Wred)';
    laplacianSpectrumReducedGraphs2(i,:) = computeLaplacianEV(G2.Wred)';
    laplacianSpectrumFullGraphs1(i,:)   = computeLaplacianEV(G1.W)';
    laplacianSpectrumFullGraphs2(i,:)   = computeLaplacianEV(G2.W)';

end


% plot log-transformed eigenvalues of full graph
eps1 = 10;

figure
plot (log(laplacianSpectrumFullGraphs1(1,:) + eps1)','r')
hold on
plot (log(laplacianSpectrumFullGraphs2(1,:)' + eps1)','b')
plot (log(laplacianSpectrumFullGraphs1' + eps1),'r')
hold on
plot (log(laplacianSpectrumFullGraphs2' + eps1),'b')
legend('Tarbert','Upper Ness')
title('Full Graph Spectrum')
xlabel('Eigenvalue Number')
ylabel('log [\lambda]')
set(gca, 'FontSize', 18)

eps2 = 100;
% plot log-transformed eigenvalues of reduced graph
figure
plot (log(laplacianSpectrumReducedGraphs1' + eps2),'r')
hold on
plot (log(laplacianSpectrumReducedGraphs2' + eps2),'b')
plot (log(laplacianSpectrumReducedGraphs1' + eps2),'r')
hold on
plot (log(laplacianSpectrumReducedGraphs2' + eps2),'b')
legend('Tarbert','Upper Ness')
title('Reduced Graph Spectrum')
xlabel('Eigenvalue Number')
ylabel('log [\lambda]')
set(gca, 'FontSize', 18)

% plot distributions of full graphs
EVfull = [laplacianSpectrumFullGraphs1; laplacianSpectrumFullGraphs2];
showDistributions(EVfull, 2, eps1, {'r', 'b'})
title('Full Graph Spectrum Distributions')

 
% plot distributions of reduced graphs
EVred = [laplacianSpectrumReducedGraphs1; laplacianSpectrumReducedGraphs2];
showDistributions(EVred, 2, eps2, {'r', 'b'})
title('Reduced Graph Spectrum Distributions')


