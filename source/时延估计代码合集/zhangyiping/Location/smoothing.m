function y=smoothing (x,N);

len=length(x);
win=hanning(2*N+1);
win1=win(1:N+1);
win2=win(N+2:2*N+1);

y1=filter(flipud(win1),[1],x);

x2=zeros(len,1);
x2(1:len-N)=x(N+1:len);

y2=filter(flipud(win2),[1],x2);

y=(y1+y2)/norm(win,2);