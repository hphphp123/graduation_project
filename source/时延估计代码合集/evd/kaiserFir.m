win=kaiser(257,3.4);
f=[0 200/20000-0.0001 200/20000 3400/20000 3400/20000+0.0001 1];
m=[0 0 1 1 0 0];
b = fir2(256,f,m,win);
[h,w] = freqz(b,1,512);
plot(f,m,w/pi,abs(h));

%%ºÏ≤‚”Ô“Ùø™ ºµ„%%
x=ones(500,1);
m = length(x);
L = 30;
for n=1:(L*(m-1)+1)
    X=(n-1)/L+1;
    if ( X == fix(X) )
       v(n) = x(X);
    else
        v(n) = 0;
    end
end
% Xfft = abs(fft(x,length(x)));
Vfft = abs(fft(v,length(v)));
for(
if 
H = 

