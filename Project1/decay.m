clear all
% po = 3411;
% to = 365*24*60*60;
% ts = 0:1:100;
% 
% for i=1:size(ts,2)
%     p(i) = 0.1*po*((ts(i) + 10)^-0.2 - (to + ts(i) + 10)^-0.2 + 0.87*(to + ts(i) + 2*10^7)^-0.2 - 0.87*(ts(i) + 2*10^7)^-0.2); 
% end
% 
% plot(ts(1,:),p(1,:))
% 
% value = 0.07*3411

hold on
file = 'RUN/data.txt';
delimiterln = ' ';
A = importdata(file,delimiterln);


plot(A(:,1)-2.0, A(:,2),'lineWidth',2)

xlabel('Time (s)','FontSize',12,'FontWeight','bold')
ylabel('Power (MW)','FontSize',12,'FontWeight','bold')
grid on
print('~/Dropbox/NE724/report/figure/power', '-dpng', '-r300');

for i=1:size(A,1)
   tsat(i) = 652.744; 
end

figure
plot(A(:,1), A(:,3),'lineWidth',2)
hold on
plot(A(:,1), A(:,4),'lineWidth',2)
plot(A(:,1), A(:,5),'lineWidth',2)
plot(A(:,1), A(:,6),'lineWidth',2)
plot(A(:,1), A(:,7),'lineWidth',2)
plot(A(:,1), A(:,8),'lineWidth',2)
plot(A(:,1), tsat(:),'--','lineWidth',2)
xlabel('Time (s)','FontSize',12,'FontWeight','bold')
ylabel('Temperature (F)','FontSize',12,'FontWeight','bold')
grid on
legend('Clad_1','Clad_2','Clad_3','Clad_4','Core Inlet','Core Exit','T_{sat}')
print('~/Dropbox/NE724/report/figure/temp', '-dpng', '-r300');
hold off

figure
plot(A(:,1), A(:,9),'lineWidth',2)
hold on
plot(A(:,1), A(:,10),'lineWidth',2)
plot(A(:,1), A(:,11),'lineWidth',2)
xlabel('Time (s)','FontSize',12,'FontWeight','bold')
ylabel('Flow Rate (lbm/hr)','FontSize',12,'FontWeight','bold')
grid on
legend('m_{core}','m_{loop1}','m_{loop2}')
print('~/Dropbox/NE724/report/figure/mass', '-dpng', '-r300');
hold off

