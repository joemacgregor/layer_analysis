function power_loss         = slope_loss(freq, speed, dist_stack, num_stack, dip_reflect)
% SLOPE_LOSS Echo-intensity loss due to unfocused stacking of a sloping reflection.
% 
% POWER_LOSS = SLOPE_LOSS(FREQ,SPEED,DIST_STACK,NUM_STACK,DIP_REFLECT)
% calculates the power loss (POWER_LOSS, dB) due to horizontal coherent
% (unfocused) stacking of a sloping reflection. FREQ (Hz) is the radar
% center frequency, SPEED (m/s) is the speed of radio-wave propagation
% through the medium, DIST_STACK (m) is the along-track distance over which
% individual traces are stacked, NUM_STACK (#) is the number of stacked
% traces and DIP_REFLECT (rad) is the relection's dip relative to
% horizontal.
% 
% SLOPE_LOSS follows Equations 1-5 and Figure 3 of:
% 
% N. Holschuh, K. Christianson and S. Anandakrishnan (2014), Power loss in
% dipping internal reflectors, imaged using ice-penetrating radar, Ann.
% Glaciol., 55(67), 49-56.
% 
% Joe MacGregor (UTIG, joemac@ig.utexas.edu)
% Last updated: 08/13/14

if (nargin ~= 5)
    error('slope_loss:nargin', 'Number of inputs not equal to five.')
end
if (~isnumeric(freq) || ~isscalar(freq) || (freq <= 0))
    error('slope_loss:freq', 'FREQ must be a numeric positive scalar.')
end
if (~isnumeric(speed) || ~isscalar(speed) || (speed <= 0))
    error('slope_loss:speed', 'SPEED must be a numeric positive scalar.')
end
if (~isnumeric(dist_stack) || ~isscalar(dist_stack) || (dist_stack <= 0))
    error('slope_loss:dist_stack', 'DIST_STACK must be a numeric positive scalar.')
end
if (~isnumeric(num_stack) || ~isscalar(num_stack) || (num_stack <= 0))
    error('slope_loss:num_stack', 'NUM_STACK must be a numeric positive scalar.')
end
if (~isnumeric(dip_reflect) || ~isscalar(dip_reflect) || (dip_reflect <= 0))
    error('slope_loss:freqw', 'DIP_REFLECT must be a numeric positive scalar.')
end
if (nargout > 1)
    error('slope_loss:nargout', 'Number of outputs greater than one.')
end

twtt_diff                   = (2 / speed) .* ((dist_stack / num_stack) * sin(dip_reflect)); % change in two-way traveltime, s
int_vec                     = repmat(0:(num_stack - 1), 1e2, 1); % integer matrix
twtt                        = linspace(0, (1 / freq))'; % dummy two-way traveltime vector for one full period, s
twtt                        = twtt(:, ones(1, num_stack)); % matrix expanded
amp                         = sin((2 * pi * freq) .* (twtt - (int_vec .* twtt_diff))); % reflection amplitude expressed as a sine wave
amp(twtt < (int_vec .* twtt_diff)) ...
                            = 0; % zero out sine "pulse" before first arrival
power_loss                  = 10 .* log10(sum(mean(amp, 2) .^ 2)); % stack horizontally, square amplitudes to get to power, sum trace and then convert to dB
power_ref                   = 10 .* log10(sum(sin((2 * pi * freq) .* twtt(:, 1)) .^ 2)); % power total at zero dip
power_loss                  = power_loss - power_ref; % normalize power loss by power at zero dip