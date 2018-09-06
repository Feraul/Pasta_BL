
a = 0:0.01:1;
b = zeros(length(a),1);
for i = 1:length(a)
    if a(i) < 0.2
        b(i) = 0.7 + 4*(a(i)^2);
    elseif a(i) >= 0.2 && a(i) < 0.4 
        b(i) = 0.6;
    elseif a(i) >= 0.4 && a(i) <= 0.6 
%        b(i) = 0.6 - 5.*((a(i) - 0.6).^2) + 8.*((a(i) - 0.6).^3);
        b(i) = 0.6 - 5.*((a(i) - 0.6).^2) + 8.*((a(i) - 0.6).^3);
    elseif a(i) > 0.6 && a(i) <= 0.8
        b(i) = 0.4;
    elseif a(i) > 0.8
        b(i) = 0.1;
    end
end

plot(a,b)
ylim([0 1]);
