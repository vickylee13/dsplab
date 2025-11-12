# dsplab
1.	DFT:
Without using built in fft functions
clc; clear;
x = input('Enter x[n]: ');
x = x(:).';                % ensure row vector
N = length(x);
n = 0:N-1; k = 0:N-1;
X = x * exp(-1j*2*pi/N).^(n'*k);   % DFT computation
disp('X[k] = ');
disp(X);
subplot(3,1,1); stem(n,x);
subplot(3,1,2); stem(k,abs(X));
subplot(3,1,3); stem(k,angle(X));

Using fft built in functions:
clc; clear;
x = input('x[n]= ');
X = fft(x);
disp(X);
n = 0:length(x)-1;
subplot(3,1,1); stem(n,x);
subplot(3,1,2); stem(n,abs(X));
subplot(3,1,3); stem(n,angle(X));



DFT USING MATRIX METHOD:
clc; clear;
x = input('x[n]= ');
N = length(x);
n = 0:N-1; k = n';
X = exp(-1j*2*pi/N).^(k*n) * x(:);
disp(X.');
subplot(3,1,1); stem(n,x);
subplot(3,1,2); stem(n,abs(X));
subplot(3,1,3); stem(n,angle(X));

2.	LINEAR AND CIRCULAR CONVOLUTION USING DFT:
LINEAR:
clc; clear;
x=input('x[n]= '); lx=input('lx= ');
h=input('h[n]= '); lh=input('lh= ');
y=conv(x,h); disp(y);
n1=lx:lx+length(x)-1;
 n2=lh:lh+length(h)-1; 
n3=lx+lh:lx+lh+length(y)-1;
subplot(3,1,1); stem(n1,x);
subplot(3,1,2); stem(n2,h);
subplot(3,1,3); stem(n3,y);


CIRCULAR:
x1=input('x1=[ ]: '); x2=input('x2=[ ]: ');
N=max(length(x1),length(x2));
y=real(ifft(fft(x1,N).*fft(x2,N))); disp(y);
n=0:N-1;
subplot(3,1,1); stem(n,[x1 zeros(1,N-length(x1))]);
subplot(3,1,2); stem(n,[x2 zeros(1,N-length(x2))]);
subplot(3,1,3); stem(n,y);

3.	DFT USING DIT 
clc; clear;
x=input('x= '); N=input('N= ');
x=[x zeros(1,N-length(x))]; x=bitrevorder(x);
for s=1:log2(N)
 m=2^s; W=exp(-1j*2*pi/m);
 for k=1:m:N
  for j=0:m/2-1
   t=W^j*x(k+j+m/2);
   x([k+j k+j+m/2])=[x(k+j)+t x(k+j)-t];
  end
 end
end
disp(x); disp(fft(x,N));

DFT USING DIF:
Above code






4.	FIR FLITER 
RECTANGULAR WINDOW:
clc; clear;
N=input('Enter filter length N: ');
fc=input('Enter cutoff frequency (0–0.5): ');
n=0:N-1; wc=2*pi*fc;
hd=sin(wc*(n-(N-1)/2))./(pi*(n-(N-1)/2)); hd((N+1)/2)=2*fc;
h=hd.*rectwin(N)'; freqz(h,1);
title('FIR Filter using Rectangular Window');

HAMMING WINDOW:
clc; clear;
N=input('Enter filter length N: ');
fc=input('Enter cutoff frequency (0–0.5): ');
n=0:N-1; wc=2*pi*fc;
hd=sin(wc*(n-(N-1)/2))./(pi*(n-(N-1)/2)); hd((N+1)/2)=2*fc;
h=hd.*hamming(N)'; freqz(h,1);
title('FIR Filter using Hamming Window');

HANNING WINDOW:
clc; clear;
N=input('Enter filter length N: ');
fc=input('Enter cutoff frequency (0–0.5): ');
n=0:N-1; wc=2*pi*fc;
hd=sin(wc*(n-(N-1)/2))./(pi*(n-(N-1)/2)); hd((N+1)/2)=2*fc;
h=hd.*hanning(N)'; freqz(h,1);
title('FIR Filter using Hanning Window');


5.	IIR FILTER USING IIT AND BT
BUTTERWORTH LOWPASS FILTER
clc; clear; close all;
wp=input('wp='); ws=input('ws='); Rp=input('Rp='); As=input('As=');
[N,wn]=buttord(wp/pi,ws/pi,Rp,As); [b,a]=butter(N,wn);
[H,w]=freqz(b,a,512);
subplot(2,2,1),plot(w/pi,abs(H)),title('Mag');
subplot(2,2,2),plot(w/pi,20*log10(abs(H))),title('dB');
subplot(2,2,3),plot(w/pi,angle(H)),title('Phase');
subplot(2,2,4),grpdelay(b,a,512),title('Grp Delay');

CHEBYSHEV FILTER:
clc; clear; close all;
wp=input('wp='); ws=input('ws='); Rp=input('Rp='); As=input('As=');
if(wp>pi||ws>pi), error('wp, ws < pi'); end
[N,wn]=cheb1ord(wp/pi,ws/pi,Rp,As);
[b,a]=cheby1(N,Rp,wn,'low');
[H,w]=freqz(b,a,512);
subplot(2,2,1),plot(w/pi,abs(H)),title('Mag');
subplot(2,2,2),plot(w/pi,20*log10(abs(H))),title('dB');
subplot(2,2,3),plot(w/pi,angle(H)),title('Phase');
subplot(2,2,4),grpdelay(b,a,512),title('Grp Delay');
6.	QUANTIZATION EFFECTS IN DIGITAL FILTER DESIGN
CO-EFFICIENT QUANTIZATION EFFECTS ON THE
FREQUENCY REPSONSE OF A DIRECT FORM IIR FILTER
clc; clear; close all;
N=input('N= '); wc=input('wc= '); Ap=input('Ap= '); As=input('As= ');
[b,a]=ellip(N,Ap,As,wc); [H,w]=freqz(b,a,512);
bq=round(b*32)/32; aq=round(a*32)/32;
[Hq,~]=freqz(bq,aq,512);
plot(w/pi,20*log10(abs(H)),'b',w/pi,20*log10(abs(Hq)),'r--'); grid on;
xlabel('\omega/\pi'); ylabel('Gain(dB)');
title('Blue: Original | Red: Quantized');
figure; zplane(b,a);

CO-EFFICIENT QUANTIZATION EFFECTS ON THE
FREQUENCY REPSONSE OF A CASCADE IIR FILTER
clc; clear; close all;
N=input('N= '); wc=input('wc= '); Ap=input('Ap= '); As=input('As= ');
[b,a]=ellip(N,Ap,As,wc);
[H,w]=freqz(b,a,512);
bq=round(b*32)/32; aq=round(a*32)/32;
[Hq,~]=freqz(bq,aq,512);
plot(w/pi,20*log10(abs(H)),'b',w/pi,20*log10(abs(Hq)),'r--'); grid on;
xlabel('\omega/\pi'); ylabel('Gain(dB)');
title('Blue: Original | Red: Quantized');
figure; zplane(b,a);

7.	DECIMATION:
clc; clear; close all;
x=input('x= '); ni=input('Start index= '); M=input('Downsample factor M= ');
n=ni:ni+length(x)-1;
m=ni:M:n(end);
y=x(1:M:end);
subplot(2,1,1); stem(n,x); title('x[n]'); xlabel('n');
subplot(2,1,2); stem(m,y); title('y[m]=x[Mn]'); xlabel('m');

 INTERPOLATION:
clc; clear; close all;
x=input('x= '); ni=input('Start index= '); L=input('Upsample factor L= ');
n=ni:ni+length(x)-1;
y=zeros(1,length(x)*L); y(1:L:end)=x;
m=ni*L:ni*L+length(y)-1;
subplot(2,1,1); stem(n,x); title('x[n]'); xlabel('n');
subplot(2,1,2); stem(m,y); title('y[m]=x[n/L]'); xlabel('m');




8.	ECG SIGNAL:
clc; clear; close all;
fs = input('Enter fs: ');
x = input('Enter ECG samples [ ]: ');
t = (0:length(x)-1)/fs;
[p,l] = findpeaks(x,'MinPeakHeight',0.5);
plot(t,x,'b',t(l),p,'ro'); grid on;
xlabel('Time(s)'); ylabel('Amplitude'); title('ECG Signal');
if length(l)>1
    disp(['Heart Rate = ', num2str(60/mean(diff(t(l)))), ' bpm']);
else
    disp('Not enough peaks');
end

9.	EEG SIGNAL:
clc; clear; close all;
fs = input('Enter sampling frequency (Hz): ');
eeg = input('Enter EEG signal samples [ ]: ');
t = (0:length(eeg)-1)/fs;

delta = bandpass(eeg,[0.5 4],fs);
theta = bandpass(eeg,[4 8],fs);
alpha = bandpass(eeg,[8 13],fs);
beta  = bandpass(eeg,[13 30],fs);

subplot(5,1,1); plot(t,eeg); title('Original EEG Signal');
subplot(5,1,2); plot(t,delta); title('Delta (0.5–4 Hz)');
subplot(5,1,3); plot(t,theta); title('Theta (4–8 Hz)');
subplot(5,1,4); plot(t,alpha); title('Alpha (8–13 Hz)');
subplot(5,1,5); plot(t,beta);  title('Beta (13–30 Hz)');
xlabel('Time (s)');


10.Arduino FFT
#include <arduinoFFT.h>
arduinoFFT FFT;
const uint16_t N = 128;
const double Fs = 1000, F = 100;
double R[N], I[N];
void setup() { Serial.begin(115200); }
void loop() {
  for (int i = 0; i < N; i++) {
    R[i] = 100 * sin(2 * PI * F * i / Fs);
    I[i] = 0;
  }
  FFT.Windowing(R, N, FFT_WIN_TYP_HAMMING, FFT_FORWARD);
  FFT.Compute(R, I, N, FFT_FORWARD);
  FFT.ComplexToMagnitude(R, I, N); 	
  for (int i = 0; i < N / 2; i++) {
    Serial.print((i * Fs) / N);
    Serial.print(" ");
    Serial.println(R[i]);
  }
  Serial.println();
  delay(1000);
}

Linear Arduino
int x[] = {1,2,3,4}, h[] = {1,1,1}, y[20];
int L = sizeof(x)/sizeof(x[0]), M = sizeof(h)/sizeof(h[0]), N;

void setup() {
  Serial.begin(115200);
  N = L + M - 1;
  for (int n = 0; n < N; n++) {
    y[n] = 0;
    for (int k = 0; k < L; k++)
      if (n - k >= 0 && n - k < M) y[n] += x[k] * h[n - k];
    Serial.print(n); Serial.print(" "); Serial.println(y[n]);
  }
}

void loop() {
  for (int i = 0; i < N; i++) Serial.println(y[i]);
  Serial.println();
  delay(1000);
}
FIR Arduino
const int N=128,M=5;
float x[N],y[N],h[M]={0.2,0.2,0.2,0.2,0.2};

void setup(){
  Serial.begin(115200);
  for(int n=0;n<N;n++){ 
    x[n]=50*sin(2*PI*5*n*0.01)+random(-10,10); 
    y[n]=0;
    for(int k=0;k<M;k++) if(n-k>=0) y[n]+=h[k]*x[n-k];
  }
  for(int i=0;i<N;i++){ Serial.print(x[i]); Serial.print(" "); Serial.println(y[i]); }
}

void loop(){
  for(int i=0;i<N;i++) Serial.println(y[i]);
  Serial.println();
  delay(1000);
}
