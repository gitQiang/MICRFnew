function [netfiles,benesss,adjfiles,nodefiles]=getFiles()

netfiles=cell(6,1);
netfiles{1}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/STRINGnetmap.txt';
netfiles{2}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/hotnet/iRefIndexm.txt'; 
netfiles{3}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/StringNew_HPRD_mnet.txt'; 
netfiles{4}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/network_inference/brainspan_net_cor.txt'; 
netfiles{5}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/network_inference/brainspan_net_top5.txt'; 
netfiles{6}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/network_inference/ComCo_PrePPI.txt';

benesss=cell(6,1);
benesss{1}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/Network_betweenness/Betweenness_edge_STRING.txt';
benesss{2}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/Network_betweenness/Betweenness_edge_iRef.txt'; 
benesss{3}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/Network_betweenness/Betweenness_edge_HPRD.txt'; 
benesss{4}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/Network_betweenness/Betweenness_edge_corr1.txt'; 
benesss{5}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/Network_betweenness/Betweenness_edge_coexp5.txt'; 
benesss{6}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/data/Network_betweenness/Betweenness_edge_Co_PrePPI.txt';

nodefiles=cell(2,1);
nodefiles{1}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/CRF_inputmeta.txt';
nodefiles{2}='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/hotnet_inputmeta.txt';

adjfiles=cell(6,1);
adjfiles{1}='adj_STRING.mat'; 
adjfiles{2}='adj_iRef.mat'; 
adjfiles{3}='adj_HPRD.mat'; 
adjfiles{4}='adj_braincor.mat'; 
adjfiles{5}='adj.mat'; 
adjfiles{6}='adj_Co_PrePPI.mat';  

end
