delta=0.25;
snr=0.1:0.3:80;
N=8;
A=10;
Rk1=0;Rk2=0;Rk3=0;
for i=1:N
    e1=A*exp(1i*(2*pi/N)*(delta+1)*(i-1));
    Rk1=Rk1+e1;
end

for i=1:N
    e2=A*exp(1i*(2*pi/N)*(delta)*(i-1));
    Rk2=Rk2+e2;
end

for i=1:N
    e3=A*exp(1i*(2*pi/N)*(delta-1)*(i-1));
    Rk3=Rk3+e3;
end

error_parabolic=zeros(size(snr));
error_quinn=zeros(size(snr));
error_macleod=zeros(size(snr));
error_j=zeros(size(snr));
error_jwbc=zeros(size(snr));
c=1000;
for i=1:length(snr)
    error_parabolic_temp=zeros(size(c));
    error_quinn_temp=zeros(size(c));
    error_macleod_temp=zeros(size(c));
    error_j_temp=zeros(size(c));
    error_jwbc_temp=zeros(size(c));
    for j=1:c
        sigma=sqrt(A^2 / snr(i));
        f1=Rk1+sigma/3*randn(1);
        f2=Rk2+sigma/3*randn(1);
        f3=Rk3+sigma/3*randn(1);

        error_parabolic_temp(j)=abs((abs(f3)-abs(f1))/(4*abs(f2)-2*abs(f1)-2*abs(f3))-delta);

        a1=real(f1/f2);
        a2=real(f3/f2);
        b1=a1/(1-a1);
        b2=a2/(1-a2);
        if b1>0 && b2>0
            error_quinn_temp(j)=abs(b2-delta);
        else
            error_quinn_temp(j)=abs(b1-delta);
        end

        d=real(f1*conj(f2)-f3*conj(f2))/real(2*(abs(f2)^2)+f1*conj(f2)+f3*conj(f2));
        error_macleod_temp(j)=abs((sqrt(1+8*(d*d))-1)/(4*d)-delta);

        error_j_temp(j)=abs(real((f1-f3)/(2*f2-f1-f3))-delta);

        error_jwbc_temp(j)=abs((tan(pi/N)/(pi/N))*real((f1-f3)/(2*f2-f1-f3))-delta);
    end
    error_parabolic(i)=sqrt(mean(error_parabolic_temp.^2));
    error_quinn(i)=sqrt(mean(error_quinn_temp.^2));
    error_macleod(i)=sqrt(mean(error_macleod_temp.^2));
    error_j(i)=sqrt(mean(error_j_temp.^2));
    error_jwbc(i)=sqrt(mean(error_jwbc_temp.^2));
    
end

hold on;
set(gca, 'YScale', 'log');ylim([0.001,1]);
plot(snr,error_parabolic,'Marker','o','MarkerFaceColor','red');
plot(snr,error_quinn,'Marker','x','MarkerFaceColor','green');
plot(snr,error_macleod,'Marker','v','MarkerFaceColor','blue');
plot(snr,error_j,'Marker','^','MarkerFaceColor','cyan');
plot(snr,error_jwbc,'Marker','d','MarkerFaceColor','magenta');
legend({'Parabolic','Quinn','Macleod','Jacobsen','Jacobsen with bias correction'},'Location','northeast');
title('RMSE Plot against SNR');xlabel('Input SNR');ylabel('RMSE');grid on;
hold off;