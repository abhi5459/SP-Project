delta=0.25;
snr=0.1:0.05:30;
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
c=10000;
for i=1:length(snr)
    error_parabolic_temp=zeros(length(c));
    error_quinn_temp=zeros(length(c));
    error_macleod_temp=zeros(length(c));
    error_j_temp=zeros(length(c));
    error_jwbc_temp=zeros(length(c));
    for j=1:c
        sigma=sqrt(A^2 / snr(i));
        f1=Rk1+sigma/50 *randn(1);
        f2=Rk2+sigma/50 *randn(1);
        f3=Rk3+sigma/50 *randn(1);

        error_parabolic_temp(j)=(abs(f3)-abs(f1))/(4*abs(f2)-2*abs(f1)-2*abs(f3));

        a1=real(f1/f2);
        a2=real(f3/f2);
        b1=a1/(1-a1);
        b2=a2/(1-a2);
        if b1>0 && b2>0
            error_quinn_temp(j)=b2;
        else
            error_quinn_temp(j)=b1;
        end

        d=real(f1*conj(f2)-f3*conj(f2))/real(2*(abs(f2)^2)+f1*conj(f2)+f3*conj(f2));
        error_macleod_temp(j)=(sqrt(1+8*(d*d))-1)/(4*d);

        error_j_temp(j)=real((f1-f3)/(2*f2-f1-f3));

        error_jwbc_temp(j)=(tan(pi/N)/(pi/N))*real((f1-f3)/(2*f2-f1-f3));
    end
    error_parabolic(i)=abs(mean(error_parabolic_temp-delta));
    error_quinn(i)=abs(mean(error_quinn_temp-delta));
    error_macleod(i)=abs(mean(error_macleod_temp-delta));
    error_j(i)=abs(mean(error_j_temp-delta));
    error_jwbc(i)=abs(mean(error_jwbc_temp-delta));   
end

hold on;
set(gca, 'YScale', 'log');
plot(snr,(error_parabolic),'Marker','o','MarkerFaceColor','red');
plot(snr,(error_quinn),'Marker','x','MarkerFaceColor','green');
plot(snr,(error_macleod),'Marker','v','MarkerFaceColor','blue');
plot(snr,(error_j),'Marker','^','MarkerFaceColor','cyan');
plot(snr,(error_jwbc),'Marker','d','MarkerFaceColor','magenta');
legend({'Parabolic','Quinn','Macleod','Jacobsen','Jacobsen with bias correction'},'Location','northeast');
title('Bias Plot against SNR');xlabel('Input SNR');ylabel('Bias');grid on;
hold off;