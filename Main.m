clc;clear;
in=Input();
w=0:50:2001;
f=in.circuit.simulator(2,w);

%%d
magnitude=abs(f.volres);
phase=angle(f.volres);
subplot(2,1,1);
plot(w,magnitude);
ylabel('magnitude');
subplot(2,1,2); 
plot(w,phase);
ylabel('phase');
suptitle('voltage frequency response bode diagram');

