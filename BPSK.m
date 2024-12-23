clear;
close all;

fb = 1e3;
bit_duration = 1/fb;
Gain = 1000;
fds = Gain*fb;
Tc = 1/fds;
fs = 3*fds;
sps = fs*bit_duration;
T = 1;
numbit = T*fb;
data_bit = randi([0 1], numbit,1 );
bpsk_vec = 2*data_bit - 1;
bpsk_signal = repelem(bpsk_vec, sps);
t1 = (0:length(bpsk_signal)-1)/fs;

%% Design rise cosine 
span = 6;
rolloff = 0.25;
rcos_filter = rcosdesign(rolloff, span, sps);
filtered_signal = conv(bpsk_signal, rcos_filter);
% Plot original and filtered signals
figure;
plot(t1, bpsk_signal);
hold on;
t2 = (0:length(filtered_signal)-1)/fs;
plot(t2, filtered_signal);
legend('Original BPSK', 'Filtered Signal');
title('BPSK Signal Before and After Raised Cosine Filtering');
xlabel('Time (s)');
ylabel('Amplitude');

%% DS
spsdss = fs*Tc;
d = ceil(length(filtered_signal)/spsdss);
rb = randi([0 1],d,1);
rb_s = 2*rb-1;
ds_samp = repelem(rb_s, spsdss);
ds_samp_eqsize = ds_samp(1:length(filtered_signal));


%% 
Spreaded_sig = ds_samp_eqsize.*filtered_signal;
Spreaded_sig = Spreaded_sig;
figure;
plot(t2, Spreaded_sig);
legend('Spreaded Signal');
title('Spreaded BPSK Signal After Raised Cosine Filtering');
xlabel('Time (s)');
ylabel('Amplitude');

%% Shift
tau = 0.2;
samples_shift = round(tau * fs);
shifted_signal = [zeros(samples_shift, 1); Spreaded_sig(1:end-samples_shift)];
[correlation, lag] = xcorr(shifted_signal, Spreaded_sig);
lag_time = lag / fs;
% Plot Original and Shifted Signals
figure;
subplot(2, 1, 1);
t_spreaded = (0:length(Spreaded_sig) - 1) / fs; 
plot(t_spreaded, Spreaded_sig, 'b', 'DisplayName', 'Original Signal');
hold on;
t_shifted_signal = (0:length(shifted_signal) - 1) / fs; 

plot(t_shifted_signal, shifted_signal, 'r', 'DisplayName', ['Shifted Signal (\tau = ', num2str(tau), ' s)']);
legend show;
title('Original and Shifted Signals');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot Correlation
subplot(2, 1, 2);
plot(lag_time, correlation, 'k');
title('Cross-Correlation');
xlabel('Lag (s)');
ylabel('Correlation');