% MTI/MTD????
%% ???????? ????
%  ?????16???????????????/??????MTI/MTD??    

%  ????????????????????XYZ?????????034
%  ?????[3000 8025 9000+(Y*10+Z)*200 8025]?4???
%  ?????[50 0 (Y*10+X+Z)*6 100]
close all; %??????
clear ; %??????
clc;
%% ????
C=3.0e8;  %??(m/s)
RF=3e7;  %???? 1.57GHz
Lambda=C/RF;    %??????
PulseNumber=16;   %?????
BandWidth=4.0e6;  %?????? ??B=1/???????? 
TimeWidth=10.0e-6; %??????
PRT=1e-3;   % ??????????(s),1000us??1/2*1000*300=150000????????
PRF=1/PRT;
Fs=1e8;  %????
NoisePower=-20;%(dB);%????????0dB?
Fc = 30e6;
% ---------------------------------------------------------------%
SampleNumber=fix(Fs*PRT);%??????????????
TotalNumber=SampleNumber*PulseNumber;%???????
BlindNumber=fix(Fs*TimeWidth);%???????????-??????
%% ???? 
TargetNumber=4;%????
SigPower(1:TargetNumber)=[5 1 100 10.25];%????,???
TargetDistance(1:TargetNumber)=[30000 80250 15800 80250];
%????,??m   ?????[3000 8025 9000+(Y*10+Z)*200 8025]
DelayNumber(1:TargetNumber)=fix(Fs*2*TargetDistance(1:TargetNumber)/C);
% ???????????????? fix???0????
TargetVelocity(1:TargetNumber)=[50 10000 0 100];
%?????? ??m/s   ?????[50 0 (Y*10+X+Z)*6 100]
TargetFd(1:TargetNumber)=2*TargetVelocity(1:TargetNumber)/Lambda; 
%?????????2v/?

%% ???????? 
number=fix(Fs*TimeWidth);%???????=??????=?????+1
if rem(number,2)~=0  %rem??
   number=number+1;
end   %?number????

Chirp = zeros(1,number); 
for i=-fix(number/2):fix(number/2)-1
    Ft = Fc*i/Fs+(1/2)*(BandWidth/TimeWidth)*(i/Fs)^2; %????????
    Chirp(i+fix(number/2)+1)=exp(1i*2*pi*Ft);
    %Chirp(i)=cos(pi*(2*Fc*i/Fs+(BandWidth/TimeWidth)*(i/Fs)^2));
    %exp(j*fi)*???????Chirp
end

coeff=conj(fliplr(Chirp));
%?Chirp?????????????????
% ????????????? h(t) = x(^*)(t_0-t) ? t_0 = 0 ???
figure(1);
subplot(2,1,1);
%???????
plot(real(Chirp));axis([0 90 -1.5 1.5]);title('??????');
%??????
subplot(2,1,2);
plot((0: Fs/number: Fs-Fs/number),abs(fft(Chirp)));
%Chirp_fft = fft(real(Chirp));
%plot((0: Fs/number: Fs/2-Fs/number),abs(Chirp_fft(1: number/2)));
%% ???????
% ???3???????
SignalAll=zeros(1,TotalNumber);%???????,??0
for k=1:TargetNumber-1 % ????????
   SignalTemp=zeros(1,SampleNumber);% ??PRT
   SignalTemp(DelayNumber(k)+1:DelayNumber(k)+number)=...
       sqrt(SigPower(k))*Chirp;
   %?????1????????????(DelayNumber(k)+1):(DelayNumber(k)+number)
   Signal=zeros(1,TotalNumber);
   for i=1:PulseNumber % 16?????
      Signal((i-1)*SampleNumber+1:i*SampleNumber)=SignalTemp; 
      %?????16?SignalTemp????
   end
   FreqMove=exp(1i*2*pi*TargetFd(k)*(0:TotalNumber-1)/Fs);
   %????????*??=????????
   Signal=Signal.*FreqMove;%?????????16???1???
   SignalAll=SignalAll+Signal;%?????????16???4???
end
%% ???4???????-------%
   fi=pi/3;
   SignalTemp=zeros(1,SampleNumber);% ????
   SignalTemp(DelayNumber(4)+1:DelayNumber(4)+number)=...
       sqrt(SigPower(4))*exp(1i*fi)*Chirp;
   %?????1????????????
   Signal=zeros(1,TotalNumber);
   for i=1:PulseNumber
      Signal((i-1)*SampleNumber+1:i*SampleNumber)=SignalTemp;
   end
   FreqMove=exp(1i*2*pi*TargetFd(4)*(0:TotalNumber-1)/Fs);
   %????????*??=????????
   Signal=Signal.*FreqMove;
   SignalAll=SignalAll+Signal;

figure(2);
subplot(2,1,1);plot(real(SignalAll),'r-');title('???????');...
    grid on;zoom on;
subplot(2,1,2);plot(imag(SignalAll));title('???????');grid on;zoom on;

%% ???????? 
SystemNoise=normrnd(0,10^(NoisePower/10),1,TotalNumber)...
    +1i*normrnd(0,10^(NoisePower/10),1,TotalNumber);
%???0?????10^(NoisePower/10)???
%% ??????
Echo=SignalAll+SystemNoise;% +SeaClutter+TerraClutter?????????
for i=1:PulseNumber   %???????,??????0
      Echo((i-1)*SampleNumber+1:(i-1)*SampleNumber+number)=0; %??????0
end
figure(3);%???????????
subplot(2,1,1);plot(real(Echo),'r-');title('????????,????0');
subplot(2,1,2);plot(imag(Echo));title('????????,????0');
%% ?????
SignalReal = real(Chirp);

%% ????=================================%
pc_time0=conv(Echo,coeff);%pc_time0?Echo?coeff???
pc_time1=pc_time0(number:TotalNumber+number-1);%????? number-1?
figure(4);%?????????
subplot(2,1,1);plot(abs(pc_time0),'r-');title('?????????,????');
%pc_time0?????
subplot(2,1,2);plot(abs(pc_time1));title('?????????,????');
%pc_time1?????
% ================================????=================================%
Echo_fft=fft(Echo,524288);
%????TotalNumber+number-1?FFT,?????????,???8192??FFT
coeff_fft=fft(coeff,524288);
pc_fft=Echo_fft.*coeff_fft;
pc_freq0=ifft(pc_fft);
figure(5);
subplot(2,1,1);plot(abs(pc_freq0(1:TotalNumber+number-1)));
title('?????????,?????');
subplot(2,1,2);
plot(abs(pc_time0(1:TotalNumber+number-1)-...
    pc_freq0(1:TotalNumber+number-1)),'r');
title('??????????');
pc_freq1=pc_freq0(number:TotalNumber+number-1);
%????? number-1?,??????(8192-number+1-TotalNumber)
% ================??????????????=================================%
for i=1:PulseNumber
      pc(i,1:SampleNumber)=pc_freq1((i-1)*SampleNumber+1:i*SampleNumber);
      %??PRT??????480???????
end
figure(6);
plot(abs(pc(1,:)));title('?????????,?????');

% ================MTI???????,???????????---?????%
for i=1:PulseNumber-1  %???????????
   mti(i,:)=pc(i+1,:)-pc(i,:);
end
figure(7);
mesh(abs(mti));title('MTI  result');

% ================MTD???????,???????????????==%
mtd=zeros(PulseNumber,SampleNumber);
for i=1:SampleNumber
   buff(1:PulseNumber)=pc(1:PulseNumber,i);
   buff_fft=fft(buff);
   mtd(1:PulseNumber,i)=buff_fft(1:PulseNumber);
end
  figure(8);mesh(abs(mtd));title('MTD  result');
  
 %% ??????
 
coeff_fft_c=zeros(1,2*524288); 
for i=1:8192 
    coeff_fft_c(2*i-1)=real(coeff_fft(i)); 
    coeff_fft_c(2*i)=imag(coeff_fft(i)); 
end 
echo_c=zeros(1,2*TotalNumber); 
for i=1:TotalNumber 
    echo_c(2*i-1)=real(Echo(i)); 
    echo_c(2*i)=imag(Echo(i)); 
end 
%% ????DSP?????????????
[fo,message] = fopen('/Users/sunhao/Desktop/Project0.1/coeff_fft_c.dat'...
    ,'wt');%?????? 
if fo < 0
    error('Failed to open coeff_fft_c.dat becauese: %s',message);
end
for i=1:2*8192 
    fprintf(fo,'%f,\r\n',coeff_fft_c(i)); 
end 
fclose(fo); 
 
[fo,message]=fopen('/Users/sunhao/Desktop/Project0.1/echo_c.dat'...
    ,'wt');%16???? 
if fo < 0
    error('Failed to open echo_c.dat becauese: %s',message);
end
for i=1:2*TotalNumber 
    fprintf(fo,'%f,\r\n',echo_c(i)); 
end 
fclose(fo); 