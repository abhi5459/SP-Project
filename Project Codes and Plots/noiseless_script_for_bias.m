delta=0:0.005:0.5;
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

error_parabolic=(abs(Rk3)-abs(Rk1))./(4*abs(Rk2)-2*abs(Rk1)-2*abs(Rk3));

error_quinn=zeros(size(delta));
a1=real(Rk1./Rk2);
a2=real(Rk3./Rk2);
b1=a1./(1-a1);
b2=a2./(1-a2);
for i=1:length(delta)
    if b1(i)>0 && b2(i)>0
        error_quinn(i)=b2(i);
    else
        error_quinn(i)=b1(i);
    end
end

d=real(Rk1.*conj(Rk2)-Rk3.*conj(Rk2))./real(2*(abs(Rk2).^2)+Rk1.*conj(Rk2)+Rk3.*conj(Rk2));
error_macleod=(sqrt(1+8*(d.*d))-1)./(4*d);

error_j=real((Rk1-Rk3)./(2*Rk2-Rk1-Rk3));

error_jwbc=(tan(pi/N)/(pi/N))*real((Rk1-Rk3)./(2*Rk2-Rk1-Rk3));
figure(1);
hold on;
set(gca, 'YScale', 'log');ylim([0.0000001, 1]);
plot(delta,abs(error_parabolic-delta),'o',delta,abs(error_quinn-delta),'x',delta,abs(error_macleod-delta),'^',delta,abs(error_j-delta),'v',delta,abs(error_jwbc-delta),'d');
legend({'Parabolic','Quinn','Macleod','Jacobsen','Jacobsen with bias correction'},'Location','northeast');
title('Bias Plot against Delta');xlabel('Delta');ylabel('Bias');grid on;
hold off;
figure(2);
set(gca, 'YScale', 'log');ylim([0.0001, 1]);
hold on;
plot(delta,abs(error_parabolic-delta)-abs(error_jwbc-delta),'o',delta,abs(error_quinn-delta)-abs(error_jwbc-delta),'x',delta,abs(error_macleod-delta)-abs(error_jwbc-delta),'^',delta,abs(error_j-delta)-abs(error_jwbc-delta),'v');
legend({'Parabolic','Quinn','Macleod','Jacobsen'},'Location','northeast');
title('Relative Bias to "Jacobsen with Bias Correction" Plot against Delta');xlabel('Delta');ylabel('Relative Bias to "Jacobsen with Bias Correction"');grid on;
hold off;