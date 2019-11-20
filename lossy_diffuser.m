clear all
close all
tic
%% Parameters
sigma_0=0.7*10^-3; % Pump beam waist at the crystal plane
N=128; % grid resolution NXN
Nx = 64; % macro pixels of the emulated diffusers NxXNx
Ncorr =32; % macro pixels on the SLM for the imperfect correction NcorrXNcorr
t1=[1:-0.01:0]; % loss strength
iter_num=200; % for disorder averaging

%% Grid
D = sigma_0*10; dx = D/N; %m
x = (-N/2:N/2-1)*dx; % following the convention x(N/2+1)=0
[X,Y] = meshgrid(x);

ROI = [round(N/2+1-2*(D/(2*pi*sigma_0))),round(N/2+1+2*(D/(2*pi*sigma_0)))]; % define the target area around N/2+1 (around q=0)

%% Pump and SPDC parameters
G_pump = exp(-(X.^2+Y.^2)/sigma_0^2); % Pump beam, W(\rho) in the paper
G_pump = G_pump/(sum(sum(G_pump))*dx^2); % normalization
G_spdc = G_pump; % W(\rho) in the paper
w_p = 404*10^(-9);w_spdc = w_p*2; %wavelength in [m]
k_pump = 2*pi/w_p; k_spdc = 2*pi/w_spdc; %2pi/m wavenumber

%%
for j=1:length(t1) %looping over loss strength
    t = t1(j);
    j
    for i=1:iter_num % disorder averaging
        %% Generating the amplitude and phase diffusers
        
        trans = 1-rand(Nx,Nx)*(1-t); %uniform distribution between t and 1, transmission coefficients for the pump
        rand_trans = kron(trans,ones(round([N,N]./([Nx,Nx])))); %moving from macropixels to pixels
        phase = 2*pi.*rand(Nx,Nx); % give random phase in each efective pixel
        rand_phase = kron(phase,ones(round([N,N]./([Nx,Nx])))); %moving from macropixels to pixels
        diffuser_pump = rand_trans.*exp(1i*rand_phase);
        diffuser_spdc = ((rand_trans).^2).*exp(1i*rand_phase); %same phase, squared amplitude
        
        %% Focusing the pump at q=0, given NcorrXNcorr macro pixels on the SLM
        part_phase = find_optimal_phase(N,Ncorr,diffuser_pump,G_pump); %see the function for details
        imperfect_correct_pump = exp(1i*part_phase);
        %% perfect correction of the pump
        perfect_correct_pump = exp(-1i*rand_phase); %perfect phase only correction
        
        %% Calculating the output after the diffuser for the different cases
        amp_and_phase_pump  = G_pump.*diffuser_pump; %no correction
        amp_and_phase_spdc = G_spdc.*diffuser_spdc; %no correction
        phase_corrected_pump = G_pump.*diffuser_pump.*perfect_correct_pump; %perfect phase correction
        phase_corrected_spdc = G_spdc.*diffuser_spdc.*perfect_correct_pump; %using the pump's phase correction!
        phase_part_corrected_pump = G_pump.*diffuser_pump.*imperfect_correct_pump; %imperfect phase correction
        phase_part_corrected_spdc = G_spdc.*diffuser_spdc.*imperfect_correct_pump; %using the pump's phase correction!

        
        %% Moving to the far-field plane via Fourier transform
        nodiff = abs(fftshift(fft2(ifftshift(G_pump)))*dx^2).^2; % Intensity/coincidences without a diffuser, note that G_pump=G_spdc
        amp_and_phase_pump = abs(fftshift(fft2(ifftshift(amp_and_phase_pump)))*dx^2).^2; %intensity without correction
        amp_and_phase_spdc = abs(fftshift(fft2(ifftshift(amp_and_phase_spdc)))*dx^2).^2; %coincidences without correction
        phase_corrected_pump = abs(fftshift(fft2(ifftshift(phase_corrected_pump)))*dx^2).^2; %perfect phase correction
        phase_corrected_spdc = abs(fftshift(fft2(ifftshift(phase_corrected_spdc)))*dx^2).^2; %perfect phase correction
        phase_part_corrected_pump = abs(fftshift(fft2(ifftshift(phase_part_corrected_pump)))*dx^2).^2; %imperfect phase correction
        phase_part_corrected_spdc = abs(fftshift(fft2(ifftshift(phase_part_corrected_spdc)))*dx^2).^2; %imperfect phase correction
        
        %% Calculating the total signal at the target area
        sum_nodiff = sum(sum(nodiff(ROI(1):ROI(2),ROI(1):ROI(2)))); % for normalization
        eta_spdc(i,j) = sum(sum(amp_and_phase_spdc(ROI(1):ROI(2),ROI(1):ROI(2))))./sum_nodiff;
        eta_pump(i,j) = sum(sum(amp_and_phase_pump(ROI(1):ROI(2),ROI(1):ROI(2))))./sum_nodiff;
        eta_pump_corrected(i,j) = sum(sum(phase_corrected_pump(ROI(1):ROI(2),ROI(1):ROI(2))))./sum_nodiff;
        eta_spdc_corrected(i,j) = sum(sum(phase_corrected_spdc(ROI(1):ROI(2),ROI(1):ROI(2))))./sum_nodiff;
        eta_spdc_part_corrected(i,j) = sum(sum(phase_part_corrected_spdc(ROI(1):ROI(2),ROI(1):ROI(2))))./sum_nodiff;
        eta_pump_part_corrected(i,j) = sum(sum(phase_part_corrected_pump(ROI(1):ROI(2),ROI(1):ROI(2))))./sum_nodiff;

    end
end

%% plot results
figure('OuterPosition',[481 249.8 568.8 584.8]);
semilogy(1-t1,mean(eta_spdc),'LineWidth',2,'Color',[0 0 0])
hold on
semilogy(1-t1,mean(eta_spdc_part_corrected),'LineWidth',2,'Color',[0 0.447058826684952 0.74117648601532])
semilogy(1-t1,mean(eta_spdc_corrected),'LineWidth',2,'Color',[0 0.498039215803146 0])

xlabel('Loss strength, s')
ylabel('\beta^(^2^)_D_C')
set(gca,'FontSize',16)
legend1=legend('Without correction','Imperfect phase-only correction','Perfect phase-only correction');
set(legend1,...
    'Position',[0.161215806999408 0.156939910351124 0.572510833967299 0.160339260447969]);
legend boxoff
ax1=gca;
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');

xlabel('Mean transmission, \mu_t')
set(ax2,'Color','none','FontSize',16,'XAxisLocation','top','XTick',...
    [0 0.2 0.4 0.6 0.8 1],'XTickLabel',{'1','0.9','0.8','0.7','0.6','0.5'},...
    'YAxisLocation','right','YTick',zeros(1,0));
set(ax2,'Position'...
    ,[0.160937954309871 0.782007077467503 0.744047615675699 0.021904761904763]);
set(ax1,'Position',[0.160952384324301 0.156349209618947 0.744047615675699 0.648888885619148]);

figure('OuterPosition',[481 249.8 568.8 584.8]);
plot(1-t1,mean(eta_spdc_part_corrected)./mean(eta_spdc),'LineWidth',2,'Color',[0 0.447058826684952 0.74117648601532])
hold on
plot(1-t1,mean(eta_spdc_corrected)./mean(eta_spdc),'LineWidth',2,'Color',[0 0.498039215803146 0])
hold on

xlabel('Loss strength, s')
ylabel('Enhancment, \eta')
set(gca,'FontSize',16)
legend('Imperfect phase-only correction','Perfect phase-only correction')
legend boxoff
ax1=gca;
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');

xlabel('Mean transmission, \mu_t')
set(ax2,'Color','none','FontSize',16,'XAxisLocation','top','XTick',...
    [0 0.2 0.4 0.6 0.8 1],'XTickLabel',{'1','0.9','0.8','0.7','0.6','0.5'},...
    'YAxisLocation','right','YTick',zeros(1,0));
set(ax2,'Position'...
    ,[0.160937954309871 0.782007077467503 0.744047615675699 0.021904761904763]);
set(ax1,'Position',[0.160952384324301 0.156349209618947 0.744047615675699 0.648888885619148]);

toc

function [opt_phase] = find_optimal_phase(N,Ncorr,diffuser,G_pump)
% Optimizing the pump's phase using the sequential algorithm (see citation
% [29])
nphases=6; %number of phases to check
phi_vec =2*pi/nphases:2*pi/nphases:2*pi; % Phase vector with equal spacing
exp_phi=exp(1i*phi_vec);
opt_phase1 = zeros(Ncorr,Ncorr);
for i=1:Ncorr %i,j loop over macro pixels on the SLM
    for j=1:Ncorr
        phase=0;
        for k=1:nphases
            phase=phase+2*pi/nphases;
            opt_phase_temp = opt_phase1;
            opt_phase_temp(i,j) = mod(opt_phase_temp(i,j)+phase,2*pi);
            opt_phase_temp = kron(opt_phase_temp,ones(round([N,N]./([Ncorr,Ncorr])))); %macro pixels to pixels
            I_pix(k)=abs(sum(sum(G_pump.*diffuser.*exp(1i*opt_phase_temp)))).^2; %the cost function is the intensity at the center of the far-field(DC term)

        end
        A_mod=exp_phi*I_pix.'; % Dot product
        phi=mod(angle(A_mod)+2*pi,2*pi); % Extracting the phase which yields the maximum value
        opt_phase1(i,j) = phi;
    end
end
opt_phase = kron(opt_phase1,ones(round([N,N]./([Ncorr,Ncorr]))));% macro pixels to pixels

end





