%function target_reflection = reflect(sample_modulated)
answer = load('sample_modulated.mat');
sample_modulated = answer.sample_modulated;
target_reflection_delay = zeros(9,9,50000); 
target_reflection = zeros(9,9,50000); 
distance = calculate_dis();
c = 3e8;
t_delay = 2*distance/c;
fs = 5e11;
Ts = 1/fs;
t_delay_cell = t_delay/Ts;
t_delay_cell = round(t_delay_cell);
for i = 1:9
    for j = 1:9
        target_reflection_delay(i,j,:) = circshift(sample_modulated,t_delay_cell(i,j));
    end
end
lambda = 3e-3;
RCS = ones(9,9);
RCS(4,5)=1000;
realgaussian = randn(50000,1);
imgaussian = randn(50000,1);
gaussian_final = 10e-5*(realgaussian + 1i*imgaussian);
%gaussian_final = (gaussian_final)';
for i = 1:9
    for j = 1:9
       target_reflection (i,j,:)= target_reflection_delay(i,j,:).*sqrt(1./distance(i,j).^4).*sqrt((lambda^2*RCS(i,j)/(4*pi)^3));
    end
end
target_reflection_final = zeros(81,50000);
count = 1;
for i = 1:9
    for j = 1:9
        a = target_reflection(i,j,:);
        b = squeeze(a);
        b = b';
        target_reflection_final(count,:)= b;
        count = count+1;
    end
end
for i = 1:81
    target_reflection_final(i,:) = target_reflection_final(i,:)+ (gaussian_final)';
end





