%% �״�ϵͳ���۵�һ��ʵ����ҵ������ѹ�����棩
clear all 
close all
clc
%% �����趨
R0 = 3000;           %Ŀ���ʼ����
c = 3e8;      
B = 10e6;            %����
f0 = 10e9;           %��Ƶ
Tp = 10e-6;          %�������ʱ��
prt = 100e-6;
fs = 100e6;          %������
tao = 2*R0/c;        %ʱ��
k = B/Tp;            %��Ƶб��
fft_num = fs*prt;    %FFT��������
n = 0:fft_num-1;     %��������
tn = n/fs;           %��1/fsΪ���������tn�����0ʱ�̵���n�β�����������ʱ��
f = (-fft_num/2:fft_num/2-1) * (fs/fft_num) ;
%% �ز��ź�����
s = rectpuls(tn-tao-Tp/2,Tp).*exp(1j*pi*k*(tn-tao).^2).*exp(-1j*2*pi*f0*tao)%;%�״�ز��ź�
figure(1);
subplot(211)
plot(tn*1e6,s)
axis([0,50,-1,1])
xlabel('t/��s');
ylabel('Amplitude');
title('(a)');
sfft = fft(s);         %����FFT
subplot(212)
plot(f*1e-6,fftshift(abs(sfft) / max(abs(sfft))));
grid on
set(gca, 'XTick', [-10:5:10]);  
axis([-5,15,0,1]);
xlabel('f/MHz');
ylabel('Amplitude');
title('(b)');
%% �ο��ź�����
refer = rectpuls(tn-Tp/2,Tp).*exp(1j*pi*k*tn.^2);  
figure(2);
subplot(211)
plot(tn*1e6,real(refer));%�ο��ź�FFT��ʱ��ͼ
axis([0,50,-1,1]);
xlabel('t/��s');
ylabel('Amplitude');
title('(a)');
refer_fft = fft(refer);
subplot(212)
plot(f*1e-6,fftshift(abs(refer_fft)) / max(abs(refer_fft)));%�ο��ź�FFT��Ƶ��ͼ
grid on
set(gca, 'XTick', [-10:5:10]);  
axis([-5,15,0,1]);
xlabel('f/MHz');
ylabel('Amplitude');
title('(b)');
%% ����ѹ������Ƶ����˷���
e = zeros(1,190000);
sfft1 = [sfft(1:5000),e,sfft(5001:10000)];%Ƶ����
refer_fft1 = [refer_fft(1:5000),e,refer_fft(5001:10000)];
fo = sfft1.*conj( refer_fft1 ); 
fo_ifft = ifft(fo); %�ز��źźͲο��źŵ�FFT��ˣ��پ���������Ҷ�任�����Ϊ����ѹ��
fo_ifft = abs(fo_ifft);
fo_ifft_db = 20*log10(fo_ifft / max(fo_ifft));%����-dB
figure(3);  
tn1 = 0:1/(20*fs):prt-1/(20*fs);
plot(tn1*c/2,fo_ifft_db);
grid on;
set(gca, 'YTick', [-13.2:9.2:-4,0]);  
axis([2900,3100,-40,0]);
xlabel('Range/m');
ylabel('Amplitude/dB');
