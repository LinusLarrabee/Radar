close all;
clear all;
clc;
% ===================================================================================%
%                                    雷达参数                                        %
% ===================================================================================%
C=3.0e8;  %光速(m/s)
f0=30e6;  %雷达射频
Lambda=C/f0;%雷达工作波长  为10m
PulseNumber=16;   %回波脉冲数 
BandWidth=4.0e6;  %发射信号带宽
TimeWidth=10.0e-6; %发射信号时宽
PRT=1e-3;   % 雷达发射脉冲重复周期(s),  最大探测距离为1/2 * 3 *10^5 = 150000M
Pf0=1/PRT;%雷达发射脉冲重复频率（Hz） 
fs=100.0e6;  %采样频率
NoisePower=-20;%(dB);%噪声功率（目标为0dB）
fc = f0-0.5*BandWidth; %雷达的起始频率

% ---------------------------------------------------------------%
SampleNumber=fix(fs*PRT);%计算一个脉冲周期的采样点数200 000 ；
TotalNumber=SampleNumber*PulseNumber;%总的采样点数200 000*16=3200 000；
BlindNumber=fix(fs*TimeWidth);%计算一个脉冲周期的盲区-遮挡样点数84；

%===================================================================================%
%                                    目标参数                                       %
%===================================================================================%
 TargetNumber=4;%目标个数
SigPower(1:TargetNumber)=[1 0 1 1];%目标功率,无量纲，幅度为[1 1 0.5 1]；
 TargetDistance(1:TargetNumber)=[30000  67500 12000 55020];%目标距离,  
 DelayNumber(1:TargetNumber)=fix(fs*2*TargetDistance(1:TargetNumber)/C);% 把目标距离换算成采样点（距离门）
TargetVelocity(1:TargetNumber)=[0 10000 1875 937.5];%目标径向速度 单位m/s  
TargetFd(1:TargetNumber)=2*TargetVelocity(1:TargetNumber)/Lambda; %计算目标多普勒// 雷达工作波长Lambda=0.1911,TargetFd=[523,1047,0,3234]
%====================================================================================%
%                                   产生线性调频信号                                  %
%====================================================================================%
 number=fix(fs*TimeWidth);%回波的采样点数=脉压系数长度=暂态点数目+1  //number=84；
if rem(number,2)~=0
   number=number+1;%保证number为偶数，若为奇数加1变偶；
end   

for i=0:fix(number)-1
   Chirp(i+1)=cos(2*pi*fc*(i/fs)+pi*(BandWidth/TimeWidth)*(i/fs)^2);%exp(j*pi*u*(t^2));
end
figure(1);subplot(2,1,1),plot(Chirp);
Chirp_fft = abs(fft(Chirp(1:number)));
subplot(2,1,2),plot((0:fs/ number:fs/2-fs/ number),abs(Chirp_fft(1: number/2)));
title('码内信号频谱');
%====================================================================================%
%                                   码内正交解调                                     %
%====================================================================================%
n = 0:fix(number)-1;
M = 131126;     %131126点fft
local_oscillator_i=cos(n*f0/fs*2*pi);%i路本振信号
local_oscillator_q=sin(n*f0/fs*2*pi);%q路本振信号
fbb_i=local_oscillator_i.*Chirp;%i路解调   先进行一个码元的求解脉冲压缩系数
fbb_q=local_oscillator_q.*Chirp;%q路解调
window=chebwin(51,40); %切比雪夫窗函数
[b,a]=fir1(50,2*BandWidth/fs,window); %b a 分别是经滤波器定义的分子分母系数向量
fbb_i=[fbb_i,zeros(1,25)];  %因为该FIR滤波器又25个采样周期的延迟，为了保证所有的有效信息全部通过滤波器，因此在信号后面扩展了25个0
fbb_q=[fbb_q,zeros(1,25)];
fbb_i=filter(b,a,fbb_i);   %I 路 Q路信号经过低通滤波器
fbb_q=filter(b,a,fbb_q);
fbb_i=fbb_i(26:end);%截取有效信息
fbb_q=fbb_q(26:end);%截取有效信息  
fbb=fbb_i+j*fbb_q;
fbb_fft_result = fft(fbb);
figure(2);subplot(2,1,1),plot(fbb_i);
xlabel('t(单位：秒)');title('雷达发射信号码内解调后I路信号');
subplot(2,1,2),plot(fbb_q);
xlabel('t(单位：秒)');title('雷达发射信号码内解调后Q路信号');
figure(3)
plot((0:fs/number:fs/2-fs/number),abs(fbb_fft_result(1:number/2)));
xlabel('频率f(单位 Hz)');title('雷达发射信号码内解调信号的频谱');
%====================================================================================%
%                                   码内脉冲压缩                                     %
%====================================================================================%

coeff=conj(fliplr(fbb));%产生脉压系数  //Chirp是一个1*84双精度的向量，fliplr沿x轴实现矩阵左右翻转(匹配滤波器h(t)等于输入信号的反转共轭延迟，这里没考虑延迟)，conj用来求复数的共轭；

%====================================================================================%
%                                   生成目标回波                                     %
%====================================================================================%
%-------------------------产生目标回波串-----------------------------------------------------------------------------------------%
SignalAll=zeros(1,TotalNumber);%所有脉冲的信号,先填0
for k=1:TargetNumber% 依次产生各个目标
   SignalTemp=zeros(1,SampleNumber);% 一个脉冲480
   SignalTemp(DelayNumber(k)+1:DelayNumber(k)+number)=sqrt(SigPower(k))*Chirp;%线性调频矩形脉冲信号的复包络；
   %目标一：41—124（41+83） 目标二和目标三：108-191（108+83）目标四：305-388（305+83）导致目标二和目标三是进行叠加
   %一个脉冲的1个目标（未加多普勒速度）DelayNumber=[40 107 107 304];
   %number=84 表示的是对发射信号的采样点数;SigPower(1:TargetNumber)=[1 1 0.25 1]目标功率,无量纲;
   Signal=zeros(1,TotalNumber);
   for i=1:PulseNumber   %回波脉冲数 
       Signal((i-1)*SampleNumber+1:i*SampleNumber)=SignalTemp;%16个周期重复赋值；
   end
   FreqMove=cos(2*pi*TargetFd(k)*(0:TotalNumber-1)/fs);
   %目标的多普勒速度*时间=目标的多普勒相移
   %这样会造成目标二与目标三的不同，虽然在图上它们显示的区间是一样的，其幅度和相移也产生了变换
   %TargetFd=[523,1047,0,3234]
   Signal=Signal.*FreqMove;%加多普勒频移；
   SignalAll=SignalAll+Signal;
end
figure(4);
subplot(2,1,1);plot(real(SignalAll),'r-');title('没有加噪声和闭锁期条件下，检测到的目标信号的实部');

%====================================================================================%
%                                   产生系统噪声信号                                  %
%====================================================================================%
%SystemNoise=normrnd(0,10^(NoisePower/10),1,TotalNumber)+j*normrnd(0,10^(NoisePower/10),1,TotalNumber);
SystemNoise=normrnd(0,10^(NoisePower/10),1,TotalNumber);
% %====================================================================================%
% %                                   总的回波信号                                     %
% %====================================================================================%
 %Echo=SignalAll+SystemNoise;% +SeaClutter+TerraClutter;加高斯噪声；
  Echo=SignalAll;% +SeaClutter+TerraClutter;加高斯噪声；
% for i=1:PulseNumber   %在接收机闭锁期,接收的回波为0
%       Echo((i-1)*SampleNumber+1:(i-1)*SampleNumber+number)=0;
%       %从1到84回波信号设置为零，即闭锁期（这里的84为信号采样点数）
%       %会导致目标一的部分回波信号40-84为零，即闭锁期不接受只发射信号；
% end
% figure(3);
% subplot(2,1,1);plot(real(Echo));title('总回波信号的实部,闭锁期为0');grid on;
% subplot(2,1,2);plot(imag(Echo));title('总回波信号的虚部,闭锁期为0');grid on;%图2始于40点，图3始于85点

%====================================================================================%
%                                   对回波目标串进行正交解调                          %
%====================================================================================%
Echo_number=PulseNumber*PRT*fs;
n=0:Echo_number-1;
local_oscillator_i=cos(n*f0/fs*pi);%I路本振信号
local_oscillator_q=cos(n*f0/fs*pi);%Q路本振信号
s_echo_i=local_oscillator_i.* Echo;%I路解调
s_echo_q=local_oscillator_q.* Echo;%Q路解调
window=chebwin(51,40);%这是采50阶cheby窗的FIR低通滤波器
[b,a]=fir1(50,2*BandWidth/fs,window);
s_echo_i=[s_echo_i,zeros(1,25)];
s_echo_q=[s_echo_q,zeros(1,25)];
s_echo_i=filter(b,a,s_echo_i);
s_echo_q=filter(b,a,s_echo_q);
s_echo_i=s_echo_i(26:end);%截取有效信息
s_echo_q=s_echo_q(26:end);%截取有效信息
s_echo_mf=s_echo_i+j*s_echo_q;
figure(4)
subplot(2,1,1),plot((0:1/fs:PulseNumber*PRT-1/fs),s_echo_i);
xlabel('t(unit:s)'); title('雷达回波信号解调后的I路信号');

subplot(2,1,2),plot((0:1/fs:PulseNumber*PRT-1/fs),s_echo_q);
xlabel('t(unit:s)'); title('雷达回波信号解调后的q路信号');

%====================================================================================%
%                      对已经正交解调的回波信号脉冲压缩                               %
%====================================================================================%
coeff_fft=fft(coeff,M);
for i=1:PulseNumber
     s_echo_fft_result=fft(s_echo_mf(1,(i-1)*PRT*fs+1:i*PRT*fs),M);
     s_pc_fft=s_echo_fft_result.*coeff_fft;
     s_pc_result(i,:)=ifft(s_pc_fft,M);    
end

figure(6);
plot(abs(s_pc_result(1,:)));%一个周期内三峰值点理论上为40 107 304

s_pc_result_1=s_pc_result;
s_pc_result_1=reshape((s_pc_result_1)',1,PulseNumber*M);   %%%%%%%%%%注意，这里由于reshape函数的算法，需要将矩阵转置才能首尾连在一起
figure(5),subplot(2,1,1),plot((0:1/fs:PulseNumber*M/fs-1/fs),abs(s_pc_result_1)),
%N_echo_frame*T_frame-ts
xlabel('t(单位:s)'),title('脉冲压缩处理后结果（实部）');
subplot(2,1,2),plot((0:1/fs:PulseNumber*M/fs-1/fs),imag(s_pc_result_1)),
xlabel('t(单位:s)'),title('脉冲压缩处理后结果（虚部）');

%====================================================================================%
%                      MTI                                                           %
%====================================================================================%

for i=1:PulseNumber-1  %滑动对消，少了一个脉冲
   mti(i,:)=s_pc_result(i+1,:)-s_pc_result(i,:);
end
figure(7);
mesh(abs(mti));title('MTI  result');



%====================================================================================%
%                      MTD                                                           %
%====================================================================================%
% mtd=zeros(PulseNumber,SampleNumber);
% for i=1:SampleNumber
%    buff(1:(PulseNumber))= s_pc_resultmt(1:PulseNumber,i);
%    buff_fft=fftshift(fft(buff)); %用fftshift将零频搬移到中间 这样可以方便观察速度正负
%    mtd(1:PulseNumber,i)=buff_fft(1:PulseNumber)';
% end
% x=0:1:SampleNumber-1;
% y=-8:1:7;%通道这样设后读出的通道数乘单位值则是速度值。
% figure(8);mesh(x,y,abs(mtd));title('MTD  result');

mtd=zeros(PulseNumber,SampleNumber);
for i=1:SampleNumber
   buff(1:(PulseNumber-1))= mti(1:(PulseNumber-1),i);
   buff_fft=fftshift(fft(buff)); %用fftshift将零频搬移到中间 这样可以方便观察速度正负
   mtd(1:PulseNumber-1,i)=buff_fft(1:PulseNumber-1)';
end
x=0:1:SampleNumber-1;
y=-7:1:8;%通道这样设后读出的通道数乘单位值则是速度值。
figure(8);mesh(x,y,abs(mtd));title('MTD  result');