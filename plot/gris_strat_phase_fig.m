% GRIS_STRAT_PHASE_FIG
% 
% Joe MacGregor
% Last updated: 01/14/14

%%
    set(0, 'DefaultFigureWindowStyle', 'default') %#ok<UNRCH>
%%
    set(0, 'DefaultFigureWindowStyle', 'docked')

%%
figure('position', [200 200 1600 800])
imagesc(wavenum_rad, (1e6 .* twtt_old), (10 .* log10(abs(data_fft))), [-70 -30])
colormap(jet)
hold on
vl                          = vline(0, 'w--');
set(vl, 'linewidth', 3)
set(gca, 'fontsize', 20)
xlabel('Doppler frequency (Hz)')
ylabel('Traveltime ({\mu}s)')
colorbar('fontsize', 20)
text(2.05, 1.6, '(dB)', 'fontsize', 20)
box on
axis ij

%%
figure('position', [200 200 1600 800])
imagesc(wavenum_rad, (1e6 .* twtt_old), (10 .* log10(abs(data_fft))), [-70 -30])
colormap(jet)
hold on
vl                          = vline(0, 'w--');
set(vl, 'linewidth', 3)
set(gca, 'fontsize', 20)
xlabel('Doppler frequency (Hz)')
ylabel('Traveltime ({\mu}s)')
colorbar('fontsize', 20)
text(2.05, 1.6, '(dB)', 'fontsize', 20)
box on
axis ij