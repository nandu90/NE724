clear all
po = 3411;
to = 365*24*60*60;
ts = 0:1:100;

for i=1:size(ts,2)
    p(i) = 0.1*po*((ts(i) + 10)^-0.2 - (to + ts(i) + 10)^-0.2 + 0.87*(to + ts(i) + 2*10^7)^-0.2 - 0.87*(ts(i) + 2*10^7)^-0.2); 
end

plot(ts(1,:),p(1,:))

value = 0.07*3411

hold on
file = 'RUN/power.txt';
delimiterln = ' ';
A = importdata(file,delimiterln);
plot(A(:,1)-2.0, A(:,2))