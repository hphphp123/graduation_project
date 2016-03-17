function y=smoothing (x,N)

len=length(x);
len1=2*N+1;
n=1:1:len1;
swCal= (n/(len1+1))-0.5;
win = 0.62-0.48 * abs(swCal)+0.38*cos(2*pi*swCal); 
% win=hanning(2*N+1);
win1=win(1:N+1);
win2=win(N+2:2*N+1);

pfWin1=flipud(win1');

for swCnt=1:len
	pfY1(swCnt)=0;
    if swCnt==1
       pfY1(swCnt)=pfWin1(1)*x(1); 
    elseif swCnt<=(N+1)
	   for swCntI=1:swCnt
		   pfY1(swCnt) = pfY1(swCnt) + pfWin1(swCntI) * x(swCnt-swCntI+1);
       end

    else
		for swCntI=1:(N+1)
			pfY1(swCnt) = pfY1(swCnt) + pfWin1(swCntI) *  x(swCnt-swCntI+1);
        end
    end
end

y1=filter(flipud(win1'),[1],x);
pfY1=pfY1';

x2=zeros(len,1);
x2(1:len-N)=x(N+1:len);

pfWin2=flipud(win2');
for swCnt=1:1:len
	pfY2(swCnt)=0;
	if swCnt==1
	   pfY2(swCnt)=pfWin2(swCnt)*x2(swCnt);
    elseif swCnt<=N 
		   for swCntI=1:1:swCnt
				pfY2(swCnt) = pfY2(swCnt) + pfWin2(swCntI) *  x2(swCnt-swCntI+1);
           end
    else
		for swCntI=1:1:N
			pfY2(swCnt) = pfY2(swCnt) + pfWin2(swCntI) *  x2(swCnt-swCntI+1);
        end
    end
end
y2=filter(flipud(win2'),[1],x2);
pfY2=pfY2';
norN= norm(win,2);

%y = (pfY1 +pfY2) / norN;
 y=(y1+y2)/norN;