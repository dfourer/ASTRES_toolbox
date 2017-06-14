% spz_delaunay_example_modular.m
%
% P. Flandrin & Ph. Depalle
% 2015, July 7th
%
% This script runs the whole process: analysis, domain selection,
% reconstruction and comparison.
%
% As a first trial, you might simply launch the script. It will run the
% whole process by default on signal 'Damped Tone - Cosine attack'
%
% To perform the process on another signal:
%   1) Go to the Signal Section, which starts at line 85 of this script
%   2) Comment the 'Damped Tone - Cosine attack section'
%   2) Uncomment the code section from the line after the title (the one
%      starting with %%) to the line x0 = x0' included.
%   2) Tune the parameters of the analysis of the uncommented section
%   3) Runs the script
%
% To learn about the role of the several flags (Variables ending with '_Flag')
% used in the script, please consult the comments from lines 50 to line 82.
%
%
% Calls:
% plot_energy
% plot_signal
% roundgauss
% sigmerge
% spz_delaunay_compare
% spz_delaunay_dom4a
% spz_delaunay_domain_image
% spz_delaunay_domain_image_energy
% spz_delaunay_domain_select
% spz_delaunay_domains_image_save
% spz_delaunay_reconstruct
% spz_delaunay_save
% spz_delaunay_signal_read
% spz_delaunay_signal_save
% Shuffle_Figures

%% Init

tic
Dir_Name = 'Tmp'; % Most of the data are stored in files saved in Folder 'Tmp'. 
if(isdir(Dir_Name) ==0)
    fprintf('\n%s: Folder not found. Please create it before running this script.\n', Dir_Name);
    break;
end
att_time = 0.0; % Legacy to be removed in a future version

Time_in_Seconds_Flag  = true;   % true -> Axes in seconds & Hz; false -> in samples & normalized frequency

Size_Calibration_Flag = true;   % Zero pad to make sure the non-integer part of block-size samples is used; false -> does not
                                % Rmk: Size_Calibration_Flag is usually set
                                % at true by deault except for very short
                                % signals. See the signal section below

Noise_Flag            = true;  % true -> Compute a realisation; false -> Reads it from a file.

% Rmk: The following flags should be all set at true in order to perform the whole
% process. After the first run, some parts might be desactivated (set at
% false).
% Example: experimenting various selections of triangles according to
% different values of ec1, and ec2, does not require to re-compute STFT and
% re-extract all triangles at each time. Then Extraction_Flag can be set at
% false after the first trial. 
% Note: changing type of signal, size of signal, block size, noise
% realisation, and Size_Calibration_Flag requires to re-run Extraction, and
% then to set Extraction_Flag at true.

Extraction_Flag       = true; % true -> Computes STFT, and performs Delaunay triangulation
Selection_Flag        = true; % true -> Selects triangles according to ec1, ec2, and Pruning_Threshold
                              % (see settings of these 3 parameters in each signal section below)
Image_Domain_Flag     = true; % true -> Computes TF zones represented by all domains found; false -> does not.
Energy_Domain_Flag    = true; % true -> Computes energy of each domain and sets a Dom_Select (see line below)
    Dom_Select = [12 13 14 10 4]; % List of indices of selected domains. Used only when Energy_Domain_Flag is at false.
    Energy_Index_Limit  = 0;  % Used only when Energy_Domain_Flag is at true; 0 means all domains used, 
                              %      otherwise chose the Energy_Index_Limit domains of highest energy

Reconstruction_Flag   = true; % true -> Reconstruct signal from selected domains; false -> does not
Mask_Flag             = true; % true -> Use selected domains, false -> no masking (transparent transform)
Difference_Flag       = true; % true -> Compute difference between original signal and recontructed one; false -> does not


% ---------------------------------------------
% Signal Examples - 1/3 - The .mat file section
% ---------------------------------------------
%%   Bat sound
%      Example_Name = 'batrich';
%      load([Example_Name '.mat']);
%      Sample_Rate = 11025;
%      Nx  = 512 ;
%      x0 = zeros(Nx+Nx/64,1) ;
%      x0((Nx/128)+50:(Nx/128)+50+length(x)-1) = (x-mean(x));
%      
%      x0 = x0/max(abs(x0)); % to make sure range lies within [-1, 1]
%      Size_Calibration_Flag = false;
%     
%      Signal_Length = length(x0);
% %     Signal_Length = 2000;
%     
%      % Nx  = 512 ; Defined right after 'load' for pre-conditioning of x0
%      SNR = 60;      % BATCH
%      ec1 = 1.8;     % BATCH
%      ec2 = 0.1;     % ec1/8;  % BATCH
%      Pruning_Threshold   = 2;
%      Energy_Index_Limit  = 0; % 0 means all domains
% 
%      x0 = x0'; 
    
  
%%   MCS 3 Modes - Mix of three modulated components
%       Example_Name = 'mcs_3modes';
%       load([Example_Name '.mat']);
%       Sample_Rate = 11025;
%       x0 = mcs;
%       
%       x0 = x0/max(abs(x0)); % to make sure range lies within [-1, 1]
%       Size_Calibration_Flag = true;
%      
%       Signal_Length = length(x0);
%     
%       Nx  = 1024;
%       SNR = 60;      % BATCH
%       ec1 = 2.0;     % BATCH
%       ec2 = 0.25;     % ec1/8;  % BATCH
%       Pruning_Threshold   = 2;
%       Energy_Index_Limit  = 0; % 0 means all domains
%       
%       x0 = x0';

%%   Damped Tone - Sharp attack
%       Example_Name = 'Damped_Tones_Noise';
%       load([Example_Name '.mat']);
%       Sample_Rate = 11025;
%       x0 = Damped_Tone_Step';
%       
%       x0 = x0/max(abs(x0)); % to make sure range lies within [-1, 1]
%       Size_Calibration_Flag = true;
%      
%       Signal_Length = length(x0);
%     
%       Nx  = 1024;
%       SNR = 60;      % BATCH
%       ec1 = 1.9;     % BATCH
%       ec2 = 0.25;     % ec1/8;  % BATCH
%       Pruning_Threshold   = 2;
%       Energy_Index_Limit  = 0; % 0 means all domains
%       
%       x0 = x0';


%%   Damped Tone - Cosine attack
      Example_Name = 'Damped_Tones_Noise';
      load([Example_Name '.mat']);
      Sample_Rate = 11025;
      x0 = Damped_Tone_Cos';
      
      x0 = x0/max(abs(x0)); % to make sure range lies within [-1, 1]
      Size_Calibration_Flag = true;
     
      Signal_Length = length(x0);
    
      Nx  = 1024;
      SNR = 60;      % BATCH
      ec1 = 2.0;     % BATCH
      ec2 = 0.2;     % ec1/8;  % BATCH
      Pruning_Threshold   = 4;
      Energy_Index_Limit  = 0; % 0 means all domains
      
      x0 = x0';


%%   Damped Tone - Sine attack
%       Example_Name = 'Damped_Tones_Noise';
%       load([Example_Name '.mat']);
%       Sample_Rate = 11025;
%       x0 = Damped_Tone_Sin';
%       
%       x0 = x0/max(abs(x0)); % to make sure range lies within [-1, 1]
%       Size_Calibration_Flag = true;
%      
%       Signal_Length = length(x0);
%     
%       Nx  = 1024;
%       SNR = 60;      % BATCH
%       ec1 = 2.0;     % BATCH
%       ec2 = 0.25;     % ec1/8;  % BATCH
%       Pruning_Threshold   = 2;
%       Energy_Index_Limit  = 0; % 0 means all domains
%       
%       x0 = x0';


% ---------------------------------------------
% Signal Examples - 2/3 - The .wav file section
% ---------------------------------------------
%%   Recorded cello excerpt (properly sampled at 11025 kHz) - All partials
%        Example_Name      = 'Cello_1';
%        [x0, Sample_Rate] = audioread([Example_Name '.wav']);
%        x0 = x0(:, 1);     % Makes the signal monophonic
%      
%        x0 = x0/max(abs(x0)); % to make sure range lies within [-1, 1]
%        Size_Calibration_Flag = false;
%      
%   %    Signal_Length = length(x0);
%        Signal_Length = 4000;
%        x0 = x0(1:Signal_Length);
%      
%        Nx  = 1024;
%        SNR = 60;      % BATCH
%        ec1 = 1.9;     % BATCH
%        ec2 = 0.3;     % ec1/8;  % BATCH
%        Pruning_Threshold   = 2;
%        Energy_Index_Limit  = 0; % 0 means all domains
%      
%        x0 = x0';


%%   Recorded cello excerpt (properly sampled at 11025 kHz) - Partial no. 3
%        Example_Name = 'Cello_1_P3';
%        [x0, Sample_Rate] = audioread([Example_Name '.wav']);
%        x0 = x0(:, 1);     % Makes the signal monophonic
%      
%        x0 = x0/max(abs(x0)); % to make sure range lies within [-1, 1]
%        Size_Calibration_Flag = false;
%      
% %      Signal_Length = length(x0);
%        Signal_Length = 2000;
%        x0 = x0(1:Signal_Length);
%      
%        Nx  = 1024;
%        SNR = 60;      % BATCH
%        ec1 = 1.6;     % BATCH
%        ec2 = 0.2;     % ec1/8;  % BATCH
%        Pruning_Threshold   = 2;
%        Energy_Index_Limit  = 0; % 0 means all domains
%     
%        x0 = x0';

    
%%   Alto_Giove - Soprano (properly sampled at 8000 Hz)
%        Example_Name = 'AltoGiove7';
%        [x0, Sample_Rate] = audioread([Example_Name '.wav']);
%        x0 = x0(:, 1);     % Makes the signal monophonic
%     
%        x0 = x0/max(abs(x0)); % to make sure range lies within [-1, 1]
%        Size_Calibration_Flag = false;
% 
% %       Signal_Length = length(x0);
%        Signal_Length = 2000;
%        x0 = x0(1:Signal_Length);
%      
%        Nx  = 1024;
%        SNR = 60;      % BATCH
%        ec1 = 1.6;     % BATCH
%        ec2 = 0.2;     % ec1/8;  % BATCH
%        Pruning_Threshold   = 2;
%        Energy_Index_Limit  = 0; % 0 means all domains
% 
%        x0 = x0';


% ---------------------------------------------
% Signal Examples - 3/3 - The synthetic section
% ---------------------------------------------
%%   Smooth Damped Sine Wave
%         Example_Name  = 'SDSW';
%         Sample_Rate   = 11025;
%         Signal_Length = 8000; 
%         
%         % Synthesis       
%         t  = 1:Signal_Length ;
%     
%         Time_0      = 1000; %3900;  % BATCH % round(Signal_Length/13); 
%         Freq_0      = 0.15;
%         Freq_Slope  = 0.00000;
%         Bandwidth   = 0.0005; % 0.00005
%     
%         Damped_Coeff= pi * Bandwidth;
%         Inst_Amp    = exp(- Damped_Coeff * (t-1));
%   
%         Arch_Cos    = ones(size(Inst_Amp));
%         att_time    = 1000; % 0 si attaque Heaviside et, par exemple, 2000 si attaque cosinus
%         Arch_Cos(1:att_time) = 0.5*(1 - cos(2*pi*(0:att_time-1)/(2*att_time)));
%         Inst_Amp    = Inst_Amp .* Arch_Cos;
%     
%         Inst_Freq   = Freq_0 + Freq_Slope*(t-1);
%         Inst_Phase  = 2*pi * (Freq_0*(t-1) + 0.5*Freq_Slope*((t-1).^2));
% 
%         x = [zeros(1, Time_0), Inst_Amp .* cos(Inst_Phase)];
%         x0 = x(1:length(t));
%         % End of Synthesis
%         
%         x0 = x0/max(abs(x0)); % to make sure range lies within [-1, 1]
%         Size_Calibration_Flag = false;
%         
%         Nx = 1024;  % Block Size (should be <= Signal_Length)
%         SNR = 60;      % BATCH
%         ec1 = 1.6;     % BATCH
%         ec2 = 0.2;     % ec1/8;  % BATCH
%         Pruning_Threshold   = 2;
%         Energy_Index_Limit  = 0; % 0 means all domains

      
Root_File_Name = [Dir_Name '/' Example_Name '_AT_' int2str(att_time)...
                           '_Nx_' int2str(Nx) '_dB_' int2str(SNR)];

%   Compute L for Gaussian window
Nfft            = 2*Nx ; %2*Nx ; % number of FFT bins
prec            = 10^(-6) ;
[H, L]          = roundgauss(Nfft, prec);
%% Modify x's length so that whole signal is being analysed
if(Size_Calibration_Flag)
    fprintf('\nSignal Size Calibration.\n\tBefore: %d', length(x0));
    
   
    %   Extends signal x0 according to the analysis parameters
    t_border        = 2*ceil(L);
    Step_Incr       = Nx - 2*t_border;
    Begin_Time_Max  = Signal_Length - Nx/64 - Nx;
    Complement_Length = Step_Incr - mod(Begin_Time_Max, Step_Incr);
    x0              = [x0 zeros(1, Complement_Length)];
    Signal_Length   = length(x0);
    Number_Blocks   = idivide(int32(Signal_Length - Nx/64), int32(Step_Incr), 'floor');
    %Effect_Length = Nx/64 + Step_Incr*(idivide(int32(Signal_Length - Nx/64), int32(Step_Incr), 'floor')-1) + Nx;
    Effect_Length   = Nx/64 + (Step_Incr * (Number_Blocks -1)) + Nx;
    Effect_Length   = cast(Effect_Length, 'double');
    
    fprintf('\tAfter:  %d\t Effective: %d\n', length(x0), Effect_Length);
else
    Effect_Length   = Signal_Length;
end



%% Add Noise
    if(Noise_Flag)
        fprintf('\nCompute Noise Realization.\n')
        b = randn(size(x0));
        spz_delaunay_signal_save(b, [Root_File_Name '_noise.spz']);
        spz_delaunay_signal_save(x0,[Root_File_Name '_signal.spz']);
    else
        b  = spz_delaunay_signal_read([Root_File_Name '_noise.spz']);
        b  = b';
%         x0 = spz_delaunay_signal_read([Root_File_Name '_signal.spz']);
%         x0 = x0';
    end
        
    x = sigmerge(x0', b',SNR) ;
    max_x0   = max(abs(x0));

if(Time_in_Seconds_Flag)
    SR_Plot = Sample_Rate;
else
    SR_Plot = 0;
end

% Plot signal to be analyzed
Fig_Num = 1;   
figure(Fig_Num); 
plot_signal(x, [1 Effect_Length], [-max_x0 max_x0], [Example_Name ' - Original signal'], SR_Plot);

soundsc(x, Sample_Rate); % Optional for hearing signal as a sound


%% Extraction of triangles via Delaunay's approach
if(Extraction_Flag)
    Fig_Num = Fig_Num + 1;
    fprintf('\nExtraction.\n');
    [tfr, triZ_Final, xZ_Final, yZ_Final] = spz_delaunay_dom4a(x, Nx, Example_Name, Fig_Num, SR_Plot);
    
    Root_File_Name_All = [Root_File_Name '_All'];
    spz_delaunay_save(triZ_Final, xZ_Final, yZ_Final, {}, Nx, Root_File_Name_All);
end

%% Selection of triangles via Delaunay's approach
if(Selection_Flag)
    [Domains, triZLM] = spz_delaunay_domain_select(triZ_Final, xZ_Final, yZ_Final, Nx, ec1, ec2, Pruning_Threshold);
  
     Root_File_Name_Sel = [Root_File_Name '_Sel' '_ec1_' num2str(ec1) '_ec2_' num2str(ec2)];
     spz_delaunay_save(triZLM, xZ_Final, yZ_Final, Domains, Nx, Root_File_Name_Sel);
end

%% Computes Image of Delaunay domains
if(Image_Domain_Flag)
    Fig_Num = Fig_Num + 1;
    Image_All = spz_delaunay_domain_image(Example_Name,triZLM, xZ_Final, yZ_Final, Domains, Nx+1, Signal_Length, Fig_Num, SR_Plot);
    
    Root_File_Name_Image = [Root_File_Name_Sel '_img.spz'];
    spz_delaunay_domains_image_save(Image_All, Nx+1, Signal_Length, Root_File_Name_Image);
end

%% Computes Energy of Delaunay domains
if(Energy_Domain_Flag)
    Fig_Num = Fig_Num + 1;
    Energy_Domains = spz_delaunay_domain_image_energy(Example_Name, tfr, triZLM, xZ_Final, yZ_Final, Domains, Fig_Num, SR_Plot);
    [Sorted_Energy_Domains, Sorted_Domains_Indices] = sort(Energy_Domains, 'descend');

    if(Energy_Index_Limit == 0)
        Energy_Index_Limit = length(Domains);
    end

    Dom_Select = Sorted_Domains_Indices(1:Energy_Index_Limit);
    Dom_Select = Dom_Select';
    Dom_Select = length(Domains) + 1 - Dom_Select;
    
    fprintf('\n Highest energy domains: %d out of %d selected.\n', Energy_Index_Limit, length(Domains));
    
    % Temporary for checking the energy decrease
    max_Energy = max(10*log10(Sorted_Energy_Domains));
    
    Fig_Num = Fig_Num + 1;
    figure(Fig_Num); 
    plot_energy(10*log10(Sorted_Energy_Domains), [1 length(Domains)], [(max_Energy-100) max_Energy], [Example_Name ' - Energy of domains in decreasing order']);
    
% To be modified   
%     Root_File_Name_Image = [Root_File_Name_Sel '_img.spz'];
%     spz_delaunay_domains_image_save(Image_All, Nx+1, Signal_Length, Root_File_Name_Image);
end

%% Reconstruction
if(Reconstruction_Flag)
    Fig_Num = Fig_Num + 1;
    x_out = spz_delaunay_reconstruct(x0, Example_Name, tfr, triZLM, xZ_Final, yZ_Final, Domains, Dom_Select, Mask_Flag, Fig_Num, SR_Plot);
    x_out = x_out/max(abs(x_out)); % Normalisation to prevent clipping when writing audio file
    
    soundsc(x_out, Sample_Rate);
    
    Root_File_Name_x_out = [Root_File_Name_Sel '_out.wav'];
    audiowrite(Root_File_Name_x_out, x_out, Sample_Rate)
end

%% Error (diff between reconstructed and noise-free signal)
if(Difference_Flag)
    Fig_Num = Fig_Num + 2;
    [diff_x, SNR_Out] = spz_delaunay_compare(x0, Example_Name, 2*L, SNR, x_out, Fig_Num, SR_Plot); 
end

% Energy_Domains_All  = sum(Energy_Domains);
% Energy_tfr_All      = sum(sum(abs(tfr.^2)))/2; % /2 as the energy was computed in the half frequency range 
% Energy_Diff_All     = Energy_tfr_All - Energy_Domains_All;
% SNR_Diff            = -10*log10(Energy_Diff_All/Energy_Domains_All);
% fprintf('\nEnergy_Domains_All: \t%f \nEnergy_tfr_All: \t%f \nEnergy_Diff_All: \t%f \nSNR_Diff: \t\t%f\n',...
%          Energy_Domains_All,  Energy_tfr_All, Energy_Diff_All, SNR_Diff);

     
% Shuffle_Figures(24, 24, true); % Might be used to have figures slightly non overlapping 
time_Elapsed = toc;
fprintf('\n\nTook : %d time\n', time_Elapsed)