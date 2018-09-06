r0 = 0.2;
R = 1.2;
pr0 = 1;
pR = 0;

r = zeros(34,1);
p = zeros(34,1);
r(1) = 0.2;
dr = 1/64;
p(1) = 1;

for i = 1:32
    r(i + 1) = r(i) + dr;
    p(i + 1) = pr0 + ((log(r(i + 1)/r0))/(log(R/r0)))*(pR - pr0);
    
    dr = 1/32;
end
p(34) = 0;
r(34) = 1.2;
r

plot(r,p)