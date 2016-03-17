%This program is designed to form broadband beam using LS method
%Using ULA


clear all;
N=4;                      %the number of the sensors
L=20;                     %the number of reference filter laps
d0=0.04;                   %the spacing the array 
c=340;                    %the speed of sound
f_l=200;                  %the lowest frequency
f_u=8000;                 %the highest frequency
step1=10;
f=[f_l:step1:f_u];           %the signal frequency band
fs=2*f_u;                 %the sampling frequency
Ts=1/fs;                  %the sampling period
omega=2*pi*f/fs;             %the signal anglular frequency
lambda_u=c/f_u;           %the wavelength of the highest frequency
p=zeros(1,N);             %the position of each sensor
step2=180;
theta=[0:pi/step2:pi];    %the angle
thetaT=pi/2;              %the steered signal angle
%d2lambda=;                %d to lambda
Theta_p1=70/180*pi;
Theta_p2=110/180*pi;
%Omega_p=
Theta_s1=pi/3;
Theta_s2=2*pi/3;
%Omega_s=
alpha=1;

for(l=1:L)
    e(l,:)=exp(-j*(l-1)*omega);
end

for(n=1:N)
    d(n)=d0*(n-1);
    tau(n,:)=d(n)*cos(theta)*fs/c;
end




%compute Q_e and a
M=L*N;
 step3=1000;
 theta=[Theta_p1:pi/step3:Theta_p2];
 theta1=[0:pi/step3:Theta_s1];
 theta2=[Theta_s2:pi/step3:pi];
for(p=1:M)
    disp(p);
    k=mod(p-1,L);
    n=floor((p-1)/L)+1;
    a=k;
    b=d(n)*fs/c;
    gamma=0;
%      syms theta;
%       omega1=omega((f_u-f_l)/step1+1);
%       omega2=omega(1);
%      g1=subs('omega1*sinc(omega((f_u-f_l)/step1+1)*(a+b*cos(theta))/pi)');
%      f1=inline(g1);
%      g2=subs('omega2*sinc(omega(1)*(a+b*cos(theta))/pi)');
%      f2=inline(g2);
%      A(p)=quad(f1,Theta_p1,Theta_p2)-quad(f2,Theta_p1,Theta_p2);
%       f=inline(''omega1*sinc(omega1*(a+b*cos(theta))/pi)'');
%       A(p)=quad(f,Theta_p1,Theta_p2,omega1,a,b);
    A1(p)=trapz(theta,omega((f_u-f_l)/step1+1)*sinc(omega((f_u-f_l)/step1+1)*(a+b*cos(theta))/pi))...
         -trapz(theta,omega(1)*sinc(omega(1)*(a+b*cos(theta))/pi));
    if b==0
        if a==0
            A(p)=cos(gamma)*(Theta_p2-Theta_p1)*(omega((f_u-f_l)/step1+1)-omega(1));
        else
            A(p)=(sin(omega((f_u-f_l)/step1+1)*a+gamma)-sin(omega(1)*a+gamma))/a*(Theta_p2-Theta_p1);
        end
    elseif abs(b)>abs(a)
        thetan=acos(-a/b);
        if thetan>=Theta_p1 & thetan<=Theta_p2
            A(p)=trapz(theta,(sin(omega((f_u-f_l)/step1+1)*(a+b*cos(theta)))*b*sqrt(1-a^2/b^2)...
                .*(theta-thetan)+sin(gamma)*(a+b*cos(theta)))...
                 ./((a+b*cos(theta))*b*sqrt(1-a^2/b^2).*(theta-thetan)))...
                 -trapz(theta,(sin(omega(1)*(a+b*cos(theta)))*b*sqrt(1-a^2/b^2)...
                 .*(theta-thetan)+sin(gamma)*(a+b*cos(theta)))...
                 ./((a+b*cos(theta))*b*sqrt(1-a^2/b^2).*(theta-thetan)));
        else
            A(p)=trapz(theta,sin(omega((f_u-f_l)/step1+1)*(a+b*cos(theta)))./(a+b*cos(theta)))...
                -trapz(theta,sin(omega(1)*(a+b*cos(theta)))./(a+b*cos(theta)));
        end
    else
        A(p)=trapz(theta,sin(omega((f_u-f_l)/step1+1)*(a+b*cos(theta)))./(a+b*cos(theta)))...
             -trapz(theta,sin(omega(1)*(a+b*cos(theta)))./(a+b*cos(theta)));
    end
 
    for(q=1:M)
        %disp(q);
        l=mod(q-1,L);
        m=floor((q-1)/L)+1;
        a=k-l;
        b=(d(n)-d(m))*fs/c;
        gamma=0;
%         syms theta1 theta2;
%         g1=subs('sinc(omega((f_u-f_l)/step1+1)*(a+b*cos(theta))/pi)');
%         f1=inline('g1');
%         g2=subs('sinc(omega(1)*(a+b*cos(theta))/pi)');
%         f2=inline('g2');
%         Q_ep(p,q)=quad(f1,Theta_p1,Theta_p2)-quad(f2,Theta_p1,Theta_p2);
%         g3=subs('sinc(omega((f_u-f_l)/step1+1)*(a+b*cos(theta1))/pi)');
%         f3=inline('g3');
%         g4=subs('sinc(omega(1)*(a+b*cos(theta1))/pi)');
%         f4=inline('g4');
%         g5=subs('sinc(omega((f_u-f_l)/step1+1)*(a+b*cos(theta2))/pi)');
%         f5=inline('g5');
%         g6=subs('sinc(omega(1)*(a+b*cos(theta2))/pi)');
%         f6=inline('g6');
%         Q_es(p,q)=quad(f3,0,Theta_s1)-quad(f4,0,Theta_s1)+quad(f5,Theta_s2,pi)-quad(f6,Theta_s2,0); 
        Q_ep1(p,q)=trapz(theta,omega((f_u-f_l)/step1+1)*sinc(omega((f_u-f_l)/step1+1)*(a+b*cos(theta))/pi))...
                 -trapz(theta,omega(1)*sinc(omega(1)*(a+b*cos(theta))/pi));
        Q_es1(p,q)=trapz(theta1,omega((f_u-f_l)/step1+1)*sinc(omega((f_u-f_l)/step1+1)*(a+b*cos(theta1))/pi))...
                 -trapz(theta1,omega(1)*sinc(omega(1)*(a+b*cos(theta1))/pi))...
                 +trapz(theta2,omega((f_u-f_l)/step1+1)*sinc(omega((f_u-f_l)/step1+1)*(a+b*cos(theta2))/pi))...
                 -trapz(theta2,omega(1)*sinc(omega(1)*(a+b*cos(theta2))/pi)); 
        if b==0
            if a==0
                Q_ep(p,q)=cos(gamma)*(Theta_p2-Theta_p1)*(omega((f_u-f_l)/step1+1)-omega(1));
                Q_es(p,q)=cos(gamma)*(Theta_s1+pi-Theta_s2)*(omega((f_u-f_l)/step1+1)-omega(1));
            else
                Q_ep(p,q)=(sin(omega((f_u-f_l)/step1+1)*a+gamma)-sin(omega(1)*a+gamma))/a*(Theta_p2-Theta_p1);
                Q_es(p,q)=(sin(omega((f_u-f_l)/step1+1)*a+gamma)-sin(omega(1)*a+gamma))/a*(Theta_s1+pi-Theta_s2);
            end
        elseif abs(b)>abs(a)
            thetan=acos(-a/b);
            if thetan>=Theta_p1 & thetan<=Theta_p2
                Q_ep(p,q)=trapz(theta,(sin(omega((f_u-f_l)/step1+1)*(a+b*cos(theta)))*b*sqrt(1-a^2/b^2)...
                        .*(theta-thetan)+sin(gamma)*(a+b*cos(theta)))...
                        ./((a+b*cos(theta))*b*sqrt(1-a^2/b^2).*(theta-thetan)))...
                        -trapz(theta,(sin(omega(1)*(a+b*cos(theta)))*b*sqrt(1-a^2/b^2)...
                        .*(theta-thetan)+sin(gamma)*(a+b*cos(theta)))...
                        ./((a+b*cos(theta))*b*sqrt(1-a^2/b^2).*(theta-thetan)));
                Q_es(p,q)=trapz(theta1,sin(omega((f_u-f_l)/step1+1)*(a+b*cos(theta1)))./(a+b*cos(theta1)))...
                          -trapz(theta1,sin(omega(1)*(a+b*cos(theta1)))./(a+b*cos(theta1)))...
                          +trapz(theta2,sin(omega((f_u-f_l)/step1+1)*(a+b*cos(theta2)))./(a+b*cos(theta2)))...
                          -trapz(theta2,sin(omega(1)*(a+b*cos(theta2)))./(a+b*cos(theta2)));       
             elseif thetan<=Theta_s1 | thetan>=Theta_s2
                Q_ep(p,q)=trapz(theta,sin(omega((f_u-f_l)/step1+1)*(a+b*cos(theta)))./(a+b*cos(theta)))...
                          -trapz(theta,sin(omega(1)*(a+b*cos(theta)))./(a+b*cos(theta))); 
                Q_es(p,q)=trapz(theta1,(sin(omega((f_u-f_l)/step1+1)*(a+b*cos(theta1)))*b*sqrt(1-a^2/b^2).*(theta1-thetan)+sin(gamma)*(a+b*cos(theta1)))...
                          ./((a+b*cos(theta1))*b*sqrt(1-a^2/b^2).*(theta1-thetan)))...
                        -trapz(theta1,(sin(omega(1)*(a+b*cos(theta1)))*b*sqrt(1-a^2/b^2).*(theta1-thetan)+sin(gamma)*(a+b*cos(theta1)))...
                        ./((a+b*cos(theta1))*b*sqrt(1-a^2/b^2).*(theta1-thetan)))...
                        +trapz(theta2,(sin(omega((f_u-f_l)/step1+1)*(a+b*cos(theta2)))*b*sqrt(1-a^2/b^2)...
                        .*(theta2-thetan)+sin(gamma)*(a+b*cos(theta2)))./((a+b*cos(theta2))*b*sqrt(1-a^2/b^2).*(theta2-thetan)))...
                        -trapz(theta2,(sin(omega(1)*(a+b*cos(theta2)))*b*sqrt(1-a^2/b^2).*(theta2-thetan)+sin(gamma)*(a+b*cos(theta2)))...
                        ./((a+b*cos(theta2))*b*sqrt(1-a^2/b^2).*(theta2-thetan)));
            end
            Q_ep(p,q)=trapz(theta,sin(omega((f_u-f_l)/step1+1)*(a+b*cos(theta)))./(a+b*cos(theta)))...
                      -trapz(theta,sin(omega(1)*(a+b*cos(theta)))./(a+b*cos(theta)));
            Q_es(p,q)=trapz(theta1,sin(omega((f_u-f_l)/step1+1)*(a+b*cos(theta1)))./(a+b*cos(theta1)))...
                      -trapz(theta1,sin(omega(1)*(a+b*cos(theta1)))./(a+b*cos(theta1)))...
                      +trapz(theta2,sin(omega((f_u-f_l)/step1+1)*(a+b*cos(theta2)))./(a+b*cos(theta2)))...
                      -trapz(theta2,sin(omega(1)*(a+b*cos(theta2)))./(a+b*cos(theta2))); 
        else
            Q_ep(p,q)=trapz(theta,sin(omega((f_u-f_l)/step1+1)*(a+b*cos(theta)))./(a+b*cos(theta)))...
                      -trapz(theta,sin(omega(1)*(a+b*cos(theta)))./(a+b*cos(theta)));
            Q_es(p,q)=trapz(theta1,sin(omega((f_u-f_l)/step1+1)*(a+b*cos(theta1)))./(a+b*cos(theta1)))...
                      -trapz(theta1,sin(omega(1)*(a+b*cos(theta1)))./(a+b*cos(theta1)))...
                      +trapz(theta2,sin(omega((f_u-f_l)/step1+1)*(a+b*cos(theta2)))./(a+b*cos(theta2)))...
                      -trapz(theta2,sin(omega(1)*(a+b*cos(theta2)))./(a+b*cos(theta2)));    
       end
    end
end
for n=1:N
    C(:,n*L-L+1:n*L)=eye(L);
end
b=zeros(L,1);
b(1)=1;
Q_LS=Q_ep+alpha*Q_es;%+alpha*Q_e;
% Q_LS1=Q_ep1+alpha*Q_es1+eye(M,M);
w=inv(Q_LS)*C'*inv(C*inv(Q_LS)*C')*(b-C*inv(Q_LS)*A')+inv(Q_LS)*A';
%w=inv((Q_LS+Q_LS')/2)*A';
w1=inv((Q_LS+Q_LS')/2)*A1';

JLS=0;
for(k=1:(f_u-f_l)/step1+1)%omega
    %disp(k);
    for(i=1:step2+1)%theta
        %disp(i);
        for(n=1:N)
            v(n*L-L+1:n*L,i)=e(:,k)*exp(-j*omega(k)*tau(n,i));                       %the array manifold
        end
        B(k,i)=w'*v(:,i);
        B1(k,i)=w1'*v(:,i);
        if abs(B(k,i))<0.1
            B(k,i)=0.1;
        end
        if abs(B1(k,i))<0.1
            B1(k,i)=0.1;
        end
        if i>=71 & i<=111
            JLS=JLS+abs(B(k ,i)-1)^2;
        else
            JLS=JLS+alpha*abs(B(k,i))^2;
        end
    end
end

theta=[0:180];
[theta,f]=meshgrid(theta,f);

figure(1);
mesh(theta,f,20*log10(abs(B)));
axis([0,180,f_l,f_u,-20,0]);
% plot(theta,20*log10(abs(B_d)));
% polardb(abs(B(1,:)));
figure(2);
mesh(theta,f,20*log10(abs(B1)));

axis([0,180,f_l,f_u,-20,0]);
% mesh(theta,d2lambda,B_amp);
% theta=theta*pi/180;
% figure(3);
% polardb(theta(1,:),20*log10(abs(B(151,:))),-40,'b');
% figure(4);
% polardb(theta(1,:),20*log10(abs(B(141,:))),-40,'b');
% figure(5);
% polardb(theta(1,:),20*log10(abs(B(131,:))),-40,'b');
% figure(6);
% polardb(theta(1,:),20*log10(abs(B(121,:))),-40,'b');
% figure(7);
% polardb(theta(1,:),20*log10(abs(B(111,:))),-40,'b');
% figure(8);
% polardb(theta(1,:),20*log10(abs(B(101,:))),-40,'b');
% figure(9);
% polardb(theta(1,:),20*log10(abs(B(91,:))),-40,'b');
% figure(10);
% polardb(theta(1,:),20*log10(abs(B(81,:))),-40,'b');