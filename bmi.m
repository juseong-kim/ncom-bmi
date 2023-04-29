%% Load spike data
load('SpikeRaster.mat')

%% Raster plots for all 30 neurons
close all
figure
NN = 5; % number of neurons to plot
NN1 = 1:NN:30;
for k = 1:NN:30
    figure('visible','off') % Figure (k+4)/5
    for j = k:k+NN-1
        N = size(spk_rast(j).neu_fr, 1); % number of trials to plot
        h(j-k+1) = subplot(NN,1,j-k+1);
        for i = 1:N % plot point at spike time for each trial on ith row (y=i)
            triali = spk_rast(j).neu_fr(i,:);
            plot(triali,i*ones(size(triali)),'k.');
            xlabel('Time (ms)'); ylabel('Trial Number');
            ylim([0 N+0.5]);
            hold on; 
        end
        title(sprintf("Neuron %d",j));
        hold off;
    end
    box off;
    sgtitle(sprintf("Raster Plot of Neurons %d-%d",k,k+NN-1),'FontWeight','bold');
    set(gcf,'Position',[10 10 1500 900]);
    exportgraphics(gcf,sprintf("Figures/Raster/raster%d_%d.png",k,k+NN-1),'BackgroundColor','none')
end

%% Cosine tuning for Neuron i to Target Location i
close all
subplot_locations = [6 4 1 3 5 7 9 8]';
subplot_w = 0.4; subplot_h = 0.17; % 0.3347 0.1243
subplot_positions = [
        0.6003    0.4553    subplot_w    subplot_h;   % E
        0.5203    0.6280    subplot_w    subplot_h;   % NE
        0.3002    0.8007    subplot_w    subplot_h;   % N
        0.0800    0.6280    subplot_w    subplot_h;   % NW
        0.0000    0.4553    subplot_w    subplot_h;   % W
        0.0800    0.2827    subplot_w    subplot_h;   % SW
        0.3002    0.1100    subplot_w    subplot_h;   % S
        0.5203    0.2827    subplot_w    subplot_h;   % SE
    ];

dest_dir = "Figures/PreLinReg";
if ~exist(dest_dir,'dir')
    mkdir(dest_dir);
end

fprintf("\nPlotted neurons: ");
Nn = 30; % number of neurons
Ntg = 8; %
frs = zeros(Nn, Ntg);
for neuron = 1:Nn % neuron number (1-30)
    figure('visible','off')
    set(gcf,'Position',[10 10 1500 900]);
    for target = 1:Ntg % target number (1-E, 2-NE, ..., 8-SE)
        subplot_i = subplot_locations(target);
        subplot(5,2,subplot_i);
        fr = plot_fr(spk_rast,neuron,target,true); % TODO output 8 firing rates (one/target)
        frs(neuron, target) = fr;
        set(gca, 'OuterPosition', subplot_positions(target,:));
    end
    exportgraphics(gcf,dest_dir+sprintf("/Neuron%d.png",neuron));
    fprintf("%d ", neuron);
end

%% Plot firing rate vs target location
targets_label = ["E (0°)", "NE (45°)", "N (90°)", "NW (135°)", ...
                        "W (180°)", "SW (225°)", "S (270°)", "SE (315°)"];

if ~exist('frs','var')
    frs = zeros(30,8);
    for neuron=1:30
        for target = 1:8
            fr = plot_fr(spk_rast,neuron,target,false);
            frs(neuron,target) = fr;
        end
    end
end


dest_dir = "Figures/PreLinReg_FT";
if ~exist(dest_dir,'dir')
    mkdir(dest_dir);
end

fprintf("\nPlotted neurons: ");
for neuron = 1:30
    figure('visible','off')
    plot(1:size(frs,2), frs(neuron,:));
    title("Firing Rate vs Target Location",sprintf("Neuron %d",neuron));
    xlabel('Target Location');
    xticklabels(targets_label);
    ylabel('Firing rate (Hz)');
    grid on
    exportgraphics(gcf,dest_dir+sprintf("/Neuron%d.png",neuron))
    fprintf("%d ", neuron);
end

%% Cosine curve vs true firing rates
X = [   1 1 0;                    % E (3)
            1 1/sqrt(2) 1/sqrt(2);    % NE (2)
            1 0 1;                     % N (1)
            1 -1/sqrt(2) 1/sqrt(2);    % NW (8)
            1 -1 0;                   % W (7)
            1 -1/sqrt(2) -1/sqrt(2);  % SW (6)
            1 0 -1;                   % S (5)
            1 1/sqrt(2) -1/sqrt(2);   % SE (4)
        ];

B = zeros(30,3);
NFr = zeros(30,8);
Neurons = zeros(30,5);

for neuron = 1:30
    y = frs(neuron, :)';
    b = regress(y,X);
    B(neuron, :) = b';
    Tx = X(:,2);
    Ty = X(:,3);
    Fr = b(1) + b(2).*Tx + b(3).*Ty;
    NFr(neuron, :) = Fr';
    Neurons(neuron, :) = [b(1) sqrt(b(2)^2+b(3)^2) atan2(b(3),b(2))*180/pi b(2) b(3)]; % B0, DR, PD_A, PD_Vx, PD_Vy
end

%% Step 4
close all
targets_label = ["E (0°)", "NE (45°)", "N (90°)", "NW (135°)", ...
                        "W (180°)", "SW (225°)", "S (270°)", "SE (315°)"];
targets = [0 45 90 135 180 225 270 315];
loc_order = zeros(30,8);
loc_order(5,:) = [8 1:7]; % neuron 5
loc_order(8,:) = [4:8 1:3]; % neuron 8
loc_order(14,:) = [3:8 1:2]; % neuron 14
loc_order(18,:) = 1:8; % neuron 18
loc_order(27,:) = [3:8 1:2]; % neuron 27
loc_order(1,:) = [2:8 1]; % neuron 1
loc_order(11,:) = 1:8; % neuron 11

dest_dir = "Figures/CosineCurve";
if ~exist(dest_dir,'dir')
    mkdir(dest_dir);
end

for neuron=[5 8 14 18 27 1 11]
    figure('visible','off')
    [f, gof] = fit(targets(:,loc_order(neuron,:))', NFr(neuron, loc_order(neuron,:))', "gauss2");
    plot(f, targets, NFr(neuron, :))
    title("Firing Rate vs Target Location",sprintf("Neuron %d, R^2= %d",neuron,gof.rsquare));
    xlabel('Target Location');
    % xticklabels(targets_label);
    ylabel('Firing rate (Hz)'); 
    grid on
    exportgraphics(gcf,dest_dir+sprintf("/Neuron%d.png",neuron))
end

%% Cursor Trajectory Reconstruction
close all
targets_label = ["E (0°)", "NE (45°)", "N (90°)", "NW (135°)", ...
                        "W (180°)", "SW (225°)", "S (270°)", "SE (315°)"];
file_names = ["E" "NE" "N" "NW" "W" "SW" "S" "SE"];
dest_dir = "Figures/CursorTrajectory";
if ~exist(dest_dir,'dir')
    mkdir(dest_dir);
end
for direction = 1:8
    figure
    FR = NFr(:,direction);
    Rl = (NFr(:,direction) - Neurons(:,1))./Neurons(:,2);
    Ra = Neurons(:,3)*pi/180;
    [uu, vu] = pol2cart(mean(Ra), mean(abs(Rl)));
    [u,v] = pol2cart(Ra,Rl);
    [ue, ve] = pol2cart((direction-1)*45*pi/180, 1);
    compass(u,v); hold on
    compass(mean(u),mean(v),'r'); hold on
    compass(ue,ve, 'g');
    title(sprintf("Target Location = %s",targets_label(direction)));
    exportgraphics(gcf,dest_dir+sprintf("/%s.png",file_names(direction)));
end

