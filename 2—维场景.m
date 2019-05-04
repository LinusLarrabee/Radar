%导入二维场景
I=double(imread('1.jpg'));
[m_I,n_I]=size(I);

%定义参数
v_t=30;
v=zeros(m_I,n_I);
scan_t=0:1.0e-7:4.0e-6;
scan_angle=pi/6:pi/60:5*pi/6;
delay_t=zeros(m_I,n_I);
distance=zeros(m_I,n_I);
A=zeros(m_I,n_I);
delay=zeros(1,41);
B=zeros(1,41);
A0=zeros(1:41);


%计算距离
for m=1:m_I
    for n=1:n_I
        format long;
        distance(m,n)=sqrt(double(m).*double(m)+double(n-(n_I+1)/2).*double(n-(n_I+1)/2));
    end
end

%假设：灰度值大于100的点赋予相应的速度30m/s
for m=1:m_I
    for n=1:n_I
        if I(m,n)>100
            v(m,n)=v_t;
        end
    end
end

%计算每一次扫描结束后得到的回波信号的时间延迟，和相应点的变化的幅度变化
%假设：灰度值大于225后，默认判定为障碍物，障碍物之后的点不再继续扫描
%延迟时间=扫描到此角度时的时间+2倍距离除光速的时间
for i=1:41
    for n=1:n_I
        for m=1:m_I
            if scan_angle(i)-atan((n-(n_I+1)/2)/m)<pi/180
               format long;
               delay_t(m,n)=double(2)*double(distance(m,n))/double(3.00000e8)+scan_t(i);
               A(m,n)=1*double(I(m,n));
               if I(m,n)>225
                  break;
               end
            end
        end
     end
end

%返回计算最远的延迟和单次的扫描幅度之和，分别为B(1，i)和A0(1,i)
for i=1:41
    for n=1:n_I
        for m=1:m_I
            if scan_angle(i)-atan((n-(n_I+1)/2)/m)<pi/180
                format long
                A0(1,i)=A0(1,i)+A(m,n);
                if B(1,i)<delay_t(m,n)
                    B(1,i)=delay_t(m,n);
                end
            end
        end
    end
end
 










