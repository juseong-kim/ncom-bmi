function fr = plot_fr(spk_rast,neuron,target, plots)
    targets = [
            1 0;                    % E (3)
            1/sqrt(2) 1/sqrt(2);    % NE (2)
            0 1;                     % N (1)
            -1/sqrt(2) 1/sqrt(2);    % NW (8)
            -1 0;                   % W (7)
            -1/sqrt(2) -1/sqrt(2);  % SW (6)
            0 -1;                   % S (5)
            1/sqrt(2) -1/sqrt(2);   % SE (4)
            ];  
    targets_label = ["E (0°)", "NE (45°)", "N (90°)", "NW (135°)", ...
                        "W (180°)", "SW (225°)", "S (270°)", "SE (315°)"];
    x = targets(target,1);
    y = targets(target,2);
    ti = (spk_rast(neuron).tDir(:,1) == x) & (spk_rast(neuron).tDir(:,2) == y);
    data = spk_rast(neuron).neu_fr(ti, :);

    Nt = size(data, 1); % number of trials
    fr_times = NaN(size(data));
    frs = NaN(size(data));
    fr_trials = zeros(1,Nt); % average firing rate for each trial

    for i = 1:Nt % plot point at spike time for each trial on ith row (y=i)
        triali = data(i,:);
    
        fired = ~isnan(triali);
        triali_fired = triali(fired);
    
        first = find(fired,1);
        fired(first) = 0;
    
        fr_times(i,:) = triali.*fired;
    
        fr = 1./diff(triali_fired);
        frs(i,fired) = fr;
   
        if plots
            plot(triali_fired,i*ones(size(triali_fired)),'k.');
        end

        xlabel('Time (ms)'); ylabel('Trial Number');
        hold on; 
    end
    
    % Merge and sort time and value of firing rates
    fr_times(fr_times == 0) = NaN;
    t_1d = fr_times(:); frs_1d = frs(:);
    t_frs = t_1d(~isnan(t_1d)); % time
    d_frs = frs_1d(~isnan(frs_1d)); % data
    td_frs = [t_frs d_frs];
    td_frs = sortrows(td_frs);
    
    % Firing rates
    yyaxis right
%     plot(td_frs(:,1),td_frs(:,2),'r.');
    edges = 0:0.05:1.4;
    N = histcounts(td_frs(:,1),edges);
    frs_hz = N/edges(2)/Nt;
    fr = mean(N(N>=prctile(N,50))); % get mean of values above 50th percentile
    if plots
        plot(edges(1,1:size(edges,2)-1),frs_hz, 'r-');
        yline(prctile(N,50),'-','50th Percentile of Firing Rates');
    end

    ylabel('Firing rate (Hz)', 'Color', 'R');
    legend off

    title(sprintf("Neuron %d, %s",neuron,targets_label(target)));
end