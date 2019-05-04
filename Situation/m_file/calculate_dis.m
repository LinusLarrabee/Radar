function distance = calculate_dis()
distance = zeros(9,9);
for i = 1:9
    for j = 1:5
        distance(i,j) = sqrt((1.5+(i-1)*0.5)^2 + ((5-j)*0.5)^2);
    end
end
for i = 1:9
    for j = 1:4
        distance(i,10-j) = distance(i,j);
    end
end
end
