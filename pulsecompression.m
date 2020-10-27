%% 雷达系统导论第一次实验作业（脉冲压缩仿真）
clear all 
close all
clc
%% 参数设定
R0 = 3000;           %目标初始距离
c = 3e8;      
B = 10e6;            %带宽
f0 = 10e9;           %载频
Tp = 10e-6;          %脉冲持续时间
prt = 100e-6;
fs = 100e6;          %采样率
tao = 2*R0/c;        %时延
k = B/Tp;            %调频斜率
fft_num = fs*prt;    %FFT采样点数
n = 0:fft_num-1;     %采样个数
tn = n/fs;           %以1/fs为采样间隔，tn代表从0时刻到第n次采样所经过的时间
f = (-fft_num/2:fft_num/2-1) * (fs/fft_num) ;
%% 回波信号生成
s = rectpuls(tn-tao-Tp/2,Tp).*exp(1j*pi*k*(tn-tao).^2).*exp(-1j*2*pi*f0*tao)%;%雷达回波信号
figure(1);
subplot(211)
plot(tn*1e6,s)
axis([0,50,-1,1])
xlabel('t/μs');
ylabel('Amplitude');
title('(a)');
sfft = fft(s);         %进行FFT
subplot(212)
plot(f*1e-6,fftshift(abs(sfft) / max(abs(sfft))));
grid on
set(gca, 'XTick', [-10:5:10]);  
axis([-5,15,0,1]);
xlabel('f/MHz');
ylabel('Amplitude');
title('(b)');
%% 参考信号生成
refer = rectpuls(tn-Tp/2,Tp).*exp(1j*pi*k*tn.^2);  
figure(2);
subplot(211)
plot(tn*1e6,real(refer));%参考信号FFT的时域图
axis([0,50,-1,1]);
xlabel('t/μs');
ylabel('Amplitude');
title('(a)');
refer_fft = fft(refer);
subplot(212)
plot(f*1e-6,fftshift(abs(refer_fft)) / max(abs(refer_fft)));%参考信号FFT的频域图
grid on
set(gca, 'XTick', [-10:5:10]);  
axis([-5,15,0,1]);
xlabel('f/MHz');
ylabel('Amplitude');
title('(b)');
%% 脉冲压缩处理（频域相乘法）
e = zeros(1,190000);
sfft1 = [sfft(1:5000),e,sfft(5001:10000)];%频域补零
refer_fft1 = [refer_fft(1:5000),e,refer_fft(5001:10000)];
fo = sfft1.*conj( refer_fft1 ); 
fo_ifft = ifft(fo); %回波信号和参考信号的FFT点乘，再经过反傅里叶变换，结果为脉冲压缩
fo_ifft = abs(fo_ifft);
fo_ifft_db = 20*log10(fo_ifft / max(fo_ifft));%幅度-dB
figure(3);  
tn1 = 0:1/(20*fs):prt-1/(20*fs);
plot(tn1*c/2,fo_ifft_db);
grid on;
set(gca, 'YTick', [-13.2:9.2:-4,0]);  
axis([2900,3100,-40,0]);
xlabel('Range/m');
ylabel('Amplitude/dB');
