delta=0:0.002:0.5;
N=8;
A=10;
Rk1=zeros(size(delta));
for i=1:N
    e1=A*exp(1i*(2*pi/N)*(delta+1)*(i-1));
    Rk1=Rk1+e1;
end

Rk2=zeros(size(delta));
for i=1:N
    e2=A*exp(1i*(2*pi/N)*(delta)*(i-1));
    Rk2=Rk2+e2;
end

Rk3=zeros(size(delta));
for i=1:N
    e3=A*exp(1i*(2*pi/N)*(delta-1)*(i-1));
    Rk3=Rk3+e3;
end

error2=(abs(Rk3)-abs(Rk1))./(4*abs(Rk2)-2*abs(Rk1)-2*abs(Rk3));

error3=-real((Rk3-Rk1)./(2*Rk2-Rk1-Rk3));

P=1.22;
error4=P*(abs(Rk3)-abs(Rk1))./(abs(Rk1)+abs(Rk2)+abs(Rk3));

Q=0.60;
error5=Q*(Rk1-Rk3)./(2*Rk2+Rk1+Rk3);

hold on;
set(gca, 'YScale', 'log');ylim([0.0001,1]);
plot(delta,abs(error2-delta),'o',delta,abs(error3-delta),'x',delta,abs(error4-delta),'^',delta,abs(error5-delta),'v');
legend({'Error-2','Error-3','Error-4','Error-5'},'Location','northeast');
title('Bias from Different Methods given in Reference Paper');xlabel('Delta');ylabel('Bias');grid on;
hold off;
