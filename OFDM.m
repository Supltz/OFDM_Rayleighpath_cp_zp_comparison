clc;
clear;

N_sc=100;      %系统子载波数（不包括直流载波）、number of subcarrierA
N_fft=128;            % FFT 长度
N_cp=16;             % 循环前缀长度、Cyclic prefix
N_symbo=N_fft+N_cp;        % 1个完整OFDM符号长度
N_c=53;             % 包含直流载波的总的子载波数、number of carriers
M=4;               %4PSK调制
SNR=0:0.5:25;         %仿真信噪比?
N_frm=20;            % 每种信噪比下的仿真帧数、frame
Nd=30;               % 每帧包含的OFDM符号数?
P_f_inter=6;      %导频间隔
data_station=[];    %导频位置
L=7;                %卷积码约束长度
tblen=12*L;           %Viterbi译码器回溯深度
stage = 3;          % m序列的阶数?
ptap1 = [1 3];      % m序列的寄存器连接方式
regi1 = [1 1 1];    % m序列的寄存器初始值


%% 基带数据数据产生
P_data=randi([0 1],1,N_sc*Nd*N_frm);


%% 信道编码（卷积码、或交织器）
%卷积码：前向纠错非线性码
%交织：使突发错误最大限度的分散化
trellis = poly2trellis(7,[133 171]);       %(2,1,7)卷积编码
code_data=convenc(P_data,trellis);


%% qpsk调制
data_temp1= reshape(code_data,log2(M),[])';             %以每组2比特进行分组，M=4
data_temp2= bi2de(data_temp1);                             %二进制转化为十进制
modu_data=pskmod(data_temp2,M,pi/M);              % 4PSK调制

%% 扩频
%————————————————————————————————————————————————————————%
%扩频通信信号所占有的频带宽度远大于所传信息必需的最小带宽
%根据香农定理，扩频通信就是用宽带传输技术来换取信噪比上的好处，这就是扩频通信的基本思想和理论依据。
%扩频就是将一系列正交的码字与基带调制信号内积
%扩频后数字频率变成了原来的m倍。码片数量 = 2（符号数）* m（扩频系数）
%————————————————————————————————————————————————————————%

code = mseq(stage,ptap1,regi1,N_sc);     % 扩频码的生成
code = code * 2 - 1;        %将1、0变换为1、-1
modu_data=reshape(modu_data,N_sc,length(modu_data)/N_sc);
spread_data = spread(modu_data,code);        % 扩频
spread_data=reshape(spread_data,[],1);

%% 插入导频
P_f=3+3*1i;                       %Pilot frequency
P_f_station=[1:P_f_inter:N_fft];%导频位置（导频位置很重要，why?）
pilot_num=length(P_f_station);%导频数量

for img=1:N_fft                        %数据位置
    if mod(img,P_f_inter)~=1          %mod(a,b)就是求的是a除以b的余数
        data_station=[data_station,img];
    end
end
data_row=length(data_station);
data_col=ceil(length(spread_data)/data_row);

pilot_seq=ones(pilot_num,data_col)*P_f;%将导频放入矩阵
data=zeros(N_fft,data_col);%预设整个矩阵
data(P_f_station(1:end),:)=pilot_seq;%对pilot_seq按行取

if data_row*data_col>length(spread_data)
    data2=[spread_data;zeros(data_row*data_col-length(spread_data),1)];%将数据矩阵补齐，补0是虚载频~
end;

%% 串并转换

data_seq=reshape(data2,data_row,data_col);
data(data_station(1:end),:)=data_seq;%将导频与数据合并

%% IFFT
ifft_data=ifft(data); 

%% 插入保护间隔、循环前缀
Tx_cd=[ifft_data(N_fft-N_cp+1:end,:);ifft_data];%把ifft的末尾N_cp个数补充到最前面
Tx_zero=[zeros(size(ifft_data(N_fft-N_cp+1:end,:)));ifft_data];%加保护间隔

[rx_c_de,Ber,x1,y1] = errorate(Tx_cd,P_f_station,pilot_seq,data_station,spread_data,code,P_data);
[rx_z_de,Ber1,x2,y2] = errorate(Tx_zero,P_f_station,pilot_seq,data_station,spread_data,code,P_data);




figure(1);
 semilogy(SNR,Ber,'r-o');
 hold on
 semilogy(SNR,Ber1,'b-o');
 legend('循环前缀','0前缀');
 xlabel('SNR');
 ylabel('BER');
 title('瑞利信道下误比特率曲线');
 hold off;

 figure(2)
 subplot(2,1,1);
 x=0:1:50;
 stem(x,P_data(1:51));
 ylabel('amplitude');
 title('发送数据（以前50个数据为例)');
 legend('调制前');

 subplot(2,1,2);
 x=0:1:50;
 stem(x,rx_c_de(1:51));
 ylabel('amplitude');
 title('译码后的数据（以前50个数据为例)');
 legend('解调后');