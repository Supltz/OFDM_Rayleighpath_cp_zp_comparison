function [rx_c_de,Ber,Rx_data2,De_Bit] = errorate(Tx_cd,P_f_station,pilot_seq,data_station,spread_data,code,P_data)
SNR=0:0.5:25;
N_fft=128;
N_cp=16;
N_sc=100;
M=4;
tblen=84;
%% 并串转换
Tx_data=reshape(Tx_cd,[],1);%由于传输需要

%% 信道（通过多经瑞利信道、或信号经过AWGN信道）
Ber=zeros(1,length(SNR));
fs = 2000;                % Sample rate (Hz)
pathDelays = [0 0.0001];  % Path delays (s)
pathPower = [0 -6];       % Path power (dB)
fD = 5;                   % Maximum Doppler shift (Hz)
rchan = comm.RayleighChannel('SampleRate',fs, ...
'PathDelays',pathDelays,'AveragePathGains',pathPower, ...
'MaximumDopplerShift',fD,'Visualization','Impulse and frequency responses');
for jj=1:length(SNR)
    rx_channel=awgn(Tx_data,SNR(jj),'measured');%切换瑞利信道
     rx_channel=rchan(rx_channel);

%% 串并转换
    Rx_data1=reshape(rx_channel,N_fft+N_cp,[]);

%% 去掉保护间隔、循环前缀
    Rx_data2=Rx_data1(N_cp+1:end,:);

%% FFT
    fft_data=fft(Rx_data2);

%% 信道估计与插值（均衡）
    data3=fft_data(1:N_fft,:);
    Rx_pilot=data3(P_f_station(1:end),:); %接收到的导频
    h=Rx_pilot./pilot_seq;
    H=interp1( P_f_station(1:end)',h,data_station(1:end)','linear','extrap');%分段线性插值：插值点处函数值由连接其最邻近的两侧点的线性函数预测。对超出已知点集的插值点用指定插值方法计算函数值

%% 信道校正
    data_aftereq=data3(data_station(1:end),:)./H;
%% 并串转换
    data_aftereq=reshape(data_aftereq,[],1);
    data_aftereq=data_aftereq(1:length(spread_data));
    data_aftereq=reshape(data_aftereq,N_sc,length(data_aftereq)/N_sc);

%% 解扩
    demspread_data = despread(data_aftereq,code);       % 数据解扩

%% QPSK解调
    demodulation_data=pskdemod(demspread_data,M,pi/M);
    De_data1 = reshape(demodulation_data,[],1);
    De_data2 = de2bi(De_data1);
    De_Bit = reshape(De_data2',1,[]);

%% （解交织）
%% 信道译码（维特比译码）
    trellis = poly2trellis(7,[133 171]);
    rx_c_de = vitdec(De_Bit,trellis,tblen,'trunc','hard');   %硬判决

%% 计算误码率
    [err, Ber(jj)] = biterr(rx_c_de(1:length(P_data)),P_data);%译码后的误码率
end
