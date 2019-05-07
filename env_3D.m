
%老哥，你就把它当做黑箱子，输入值分别为1000个点中的每个点的采样间隔时间，采样点的幅度，
%哪一个角度的（i_n=1-40对应于角度值为30度到150度）的采样总点数
%输出值为扫描点的幅度，（每个网格扫了1000次）以及扫描点此时的延迟时间
function [delay_N_s,A_N_s]=env_3D(T_s,A,i_n)
%构建实际三维场景，车高设为1,5m，上下倾角定义为88度（22*pi/45）和60度,车的起始点设置为（25,0,1.5）处,有物体的话，rcs设置为1000，其余处为1
%计算可以得到二维地面的最大范围为2.6—43m
%----------构建场景----------
delay_N_s=1:1000000;
A_N_s=1:1000000;
x_I=50;
y_I=50;
x=1:x_I;
y=1:y_I;
z=zeros(x_I,y_I);
RCS=ones(x_I,y_I);
for m=1:x_I
    for n=1:y_I
        if (m>44)&&(n>27)
        if (m<47)&&(n<32)
            z(m,n)=1.5;
            RCS(m,n)=100;
        end
        end
        if (m>10)&&(n>10)
        if (m<30)&&(n<30)
            z(m,n)=1.5;
            RCS(m,n)=100;
        end
        end
    end
end

mesh(x,y,z);

%----------扫描过程----------
scan_angle=11*pi/60:pi/60:5*pi/6;
scan_t=0:1.0e-6:3.9e-5;
at_r=zeros(x_I,y_I);
at_angle=zeros(x_I,y_I);
delay_t=zeros(x_I,y_I);
distance=zeros(x_I,y_I);
h_0=zeros(40,500);
h_1=zeros(40,500);
store_car=1:40;
c=3.0e8;
k=1;

for i=1:40
for m=1:x_I
    for n=1:y_I
        if z(m,n)==0
        format long;
        distance(m,n)=double(sqrt((m-x_I/2)^2+double(n)^2+2.25));
        delay_t(m,n)=2*distance(m,n)/c;
        end
        if z(m,n)>0
        format long;
        distance(m,n)=(1.5-(1.5*tan(22*pi/45)-double(sqrt((m-x_I/2)^2+n^2)))/tan(22*pi/45))/cos(22*pi/45);
        delay_t(m,n)=2*distance(m,n)/c;
        end
    end
end
end
%扫点
for m=1:x_I
        for n=1:x_I
            at_r(m,n)=n/(m-x_I/2);
            at_angle(m,n)=atan(at_r(m,n));
            if at_angle(m,n)<0
                at_angle(m,n)=at_angle(m,n)+pi;
            end
            
        end
end
for i=1:40
    for m=1:x_I
        for n=1:y_I
            if abs(scan_angle(i)-at_angle(m,n))<pi/180
               h_0(i,k)=RCS(m,n);
               h_1(i,k)=delay_t(m,n)+scan_t(i);
               k=k+1;
            end
        end
    end
    k=1;
end

for i=1:40
    for k=1:500
        if (h_0(i,k)>1)||(k==500)
            store_car(i)=k;
            break
        end
    end
end
n=1;
for i=1:40
    if i==i_n
    for k=1:store_car(i)
        for j=1:1000
            delay_N_s(n)=T_s*j+h_1(i,k);
            A_N_s(n)=A*h_0(i,k);
            n=n+1;
        end
    end
    n=1;
    end
end

