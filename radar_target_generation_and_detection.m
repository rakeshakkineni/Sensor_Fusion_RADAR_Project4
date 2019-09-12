clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
init_range  = 110;
init_vel = -15;


%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
c = 3e8;                % speed of light 
Rmax = 200; 
range_res = 1;
Vmax = 70;
vel_res = 3;

B = c/2*range_res;
Tchirp = 5.5*2*Rmax/c;
slope  = B/Tchirp;

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq

                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    r = norm(init_range+init_vel*t(i));
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    tau = 2*r/c;
    Tx(i) = cos(2*pi*(fc*t(i)+(slope*t(i)^2/2)));
    Rx(i)  =cos(2*pi*(fc*(t(i)-tau)+(slope*(t(i)-tau)^2/2)));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i)*Rx(i);
    
end

%% RANGE MEASUREMENT


 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix_RFFT = reshape(Mix,Nr,Nd);
 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
RFFT = fft(Mix_RFFT(1:Nr,1),Nr);
 % *%TODO* :
% Take the absolute value of FFT output
RFFT = abs(RFFT);
 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
RFFT = RFFT(1:Nr/2)/Nr;

%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

 % *%TODO* :
 % plot FFT output 
f = linspace(0,400,Nr/2)*((Nr/2)/400);
plot(f,RFFT) 

 
axis ([0 200 0 1]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr = 2;
Td = 4;
% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 1;
Gd = 2;
% *%TODO* :
% offset the threshold by SNR value in dB
thold_SNR = 6;
% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
%noise_level = zeros(2*Tr+2*Gr+1,2*Td+2*Gd+1);
noise_level = zeros(Nr/2,Nd);
cfar_thold = zeros(Nr/2,Nd);
RDM_Thold = zeros(Nr/2,Nd);
% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.
for i = 1:(Nr/2-Tr-Gr-1)
    for j = 1:(Nd-Td-Gd-1)        
        rmin = (i-Tr-Gr-1);
        rmax = (i+Tr+Gr+1);        
        dmin = (j-Td-Gd-1);
        dmax = (j+Td+Gd+1);
        rmin_g = (i-Gr-1);
        rmax_g = (i+Gr+1);        
        dmin_g = (j-Gd-1);
        dmax_g = (j+Gd+1);
        
        if(rmin<1) rmin =1; end
        if(dmin<1) dmin =1; end
        if(rmin_g<1) rmin_g=1; end
        if(dmin_g<1) dmin_g=1; end
        
        if(rmax>Nr/2) rmax = Nr/2; end
        if(rmax_g>Nr/2) rmax_g = Nr/2; end
        if(dmax>Nd) dmax = Nd; end
        if(dmax_g>Nd) dmax_g = Nd; end
        sig_fft2_abs = abs(sig_fft2);   
        
        %z=sum(  reshape( a( 1:mid-0f-row-,1:allcolums,:) [],1 )  );
        windowed_T = sig_fft2_abs(rmin:rmax,dmin:dmax);
        noise_level(i,j) = sum(windowed_T(:));
        windowed_G =sig_fft2_abs(rmin_g:rmax_g,dmin_g:dmax_g);
        noise_level(i,j) = noise_level(i,j)- sum(windowed_G(:));
        noise_level(i,j) = noise_level(i,j)/((2*Tr+2*Gr+1)*(2*Td+2*Gd+1) - (2*Gr+1)*(2*Gd+1));
        cfar_thold(i,j) = pow2db(noise_level(i,j))+thold_SNR;
        if(RDM(i,j)>cfar_thold(i,j))
            RDM_Thold(i,j) = 1;
        else
            RDM_Thold(i,j)=0;
        end
    end
end


   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR





% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 
RDM_Thold(1,:) =0;
RDM_Thold(:,1) =0;
RDM_Thold(Nr/2,:) =0;
RDM_Thold(:,Nd) =0;

% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,RDM_Thold);
colorbar;


 
 
