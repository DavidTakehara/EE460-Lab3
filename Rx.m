%
%   David Takehara
%   5/20/19
%   EE 460 - Lab 6
%   Digital Receiver Project (work on indexing for decoding)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [decoded_msg y] = Rx(r,rolloff,desired_user)
    beta = rolloff;                             %   beta
    user = desired_user;                        %   user to decode
    DEBUG = 0;                                  %   debug (Y == 1, N == 0)
    f_if = 2000000;                             %   intermediate mod. frequency
    Fs = 850000;                                %   sampling frequency
    Ts = 1/Fs;                                  %   sampling period
    T_sym = 0.0000064;                          %   nominal sym. period
    preamble = 'A0Oh well whatever Nevermind';  %   decoded preamble
    pre_lng = 245;                              %   length of preabmle
    data_lng = 875;                             %   length of data
    frame_lng = 2870;                           %   length of whole frame

    %%  Frequency and Phase Tracking
    %   PLL Preprocess (squaring Rx signal)
    q=r.^2;                           % square nonlinearity  
    fl=230; ff=[0 235000 247500 252500 265000 Fs/2]/(Fs/2); % BPF center frequency at 250 kHz
    fa=[0 0 1 1 0 0];                 % which is twice f_0
    b=firpm(fl,ff,fa);                % BPF tor freq recovery
    rp=filter(b,1,q);                 % filter gives preprocessed r

    if DEBUG == 1
        figure (1);
        plotspec(rp,Ts);
        title('Output of SQing');
    else
    end

    %   Determine phase of filter
    fftrBPF=fft(rp);                     % spectrum of rBPF
    [m,imax]=max(abs(fftrBPF(1:end/2))); % find freq of max peak
    ssf=(0:length(rp))/(Ts*length(rp));  % frequency vector
    freqS=ssf(imax);                     % freq at the peak
    [IR,f]=freqz(b,1,length(rp),1/Ts);   % freq response of filter (change back to h if IIR doesnt work)
    [mi,im]=min(abs(f-freqS));           % freq where peak occurs
    bpfPhase=angle(IR(im));              % < of BPF at peak freq

    time=length(r)/Fs; t=0:Ts:time-Ts;   % time vector
    mu1 = 0.005; mu2 = .00005;            % algorithm stepsizes

    f0 = f_if - 2*Fs;                    % assumed freq at receiver(aliased to 300k)

    lng_t=length(t); th1=zeros(1,lng_t);    % initialize estimates
    th2=zeros(1,lng_t); carest=zeros(1,lng_t);

    for k=1:lng_t-1                          % Dual PLL Here
      th1(k+1)=th1(k)-mu1*rp(k)*sin(4*pi*f0*t(k)+2*th1(k)+bpfPhase);           
      th2(k+1)=th2(k)-mu2*rp(k)*sin(4*pi*f0*t(k)+2*th1(k)+2*th2(k)+bpfPhase);  
    end
    carest = cos(2*pi*f0*t+th1+th2);  %   estimated Rx oscillator
    if DEBUG == 1                           %   plot if DEBUGGING
        figure (2);
        subplot(2,1,1), plot(t,th1);        % plot first theta
        title('output of first PLL');
        ylabel('\theta_1');
        subplot(2,1,2), plot(t,th2);        % plot second theta
        title('output of second PLL');
        ylabel('\theta_2');
    else
    end

    %%  Demod w/ Recovered Carrier
    r_demod = r .* carest';                 %   demodulate from IF (column)

    if DEBUG == 1
        figure (3);                         %   plot Rx signal
        plotspec(r,Ts);
        title('Plotspec of recieved signal');
        figure (4);                         %   plot demod Rx signal
        plotspec(r_demod,Ts);
        title('Plotspec of demod recieved signal');
    else
    end

    %%  Matched Filtering
    M = T_sym/Ts;                    % downsample factor
    ps = srrc(4,beta,M,0);
    r_mfilter = conv(ps,r_demod');   % es now a row vector
    if DEBUG == 1
        figure (5);
        plotspec(r_mfilter,Ts);
        title('Plotspec of matched filter output');
    else
    end

    %%  Clock Recovery & Downsampling (fix delt and mu values, es broken)
    n = floor(lng_t/M);
    l = 4;                                  % half length of pulse shape
    tnow=l*M+1; tau=0; rs=zeros(1,n);       % initialize variables
    tausave=zeros(1,n); tausave(1)=tau; i=0;
    mu = 9000000;                            % algorithm stepsize
    delta = 1e-8;                           % time for derivative
    while tnow<length(r_mfilter)-l*M        % MAX Output Power - run iterations, rs is recovered received sig.
      i=i+1;
      rs(i)=interpsinc(r_mfilter,tnow+tau,l);           % interp at tnow+tau
      x_deltap=interpsinc(r_mfilter,tnow+tau+delta,l);  % value to right
      x_deltam=interpsinc(r_mfilter,tnow+tau-delta,l);  % value to left
      dx=x_deltap-x_deltam;                             % numerical derivative
      tau=tau+mu*dx*rs(i);                              % alg update (energy)
      tnow=tnow+M; tausave(i)=tau;                      % save for plotting
    end

    if DEBUG == 1
        figure (6);
        subplot(2,1,1), plot(rs(1:i-2),'b.')    % plot constellation diagram of rs
        title('constellation diagram (timing recovery)');
        ylabel('estimated symbol values')
        subplot(2,1,2), plot(tausave(1:i-2))    % plot trajectory of tau
    else
    end

    %%  Equilization             need to finish (alg. broken 
    s_train = letters2pam2(preamble);   % known encoded marker sequence
    n=10;
    f = zeros(n,1);                     % initialize equalizers at 0
    mu1=.01;  mu2 = 0.0025;          % stepsize 
    delta=4;                            % delay delta

    find_mrkr = xcorr(s_train,rs);      % locate preamble
    [~,strt_data] = max(abs(find_mrkr));      % locate preamble
    strt_data=mod(strt_data,frame_lng):frame_lng:length(rs)-frame_lng;
    if strt_data(1)<200
        strt_data=strt_data(2:end);
    end
        
    if min(find_mrkr) > max(find_mrkr)               % inversion check
        find_mrkr = -find_mrkr;
        rs = -rs;
    end

    if DEBUG == 1
        figure(7);
        plot(find_mrkr);
        title('Correlation of known Marker and Received');
    end

    r_eq = 0;                               % initialize w/ offset of +1
    k = 1;
    err = zeros(1,length(rs));              % error tracking term
    pre_eq = zeros(1,pre_lng);

    % start of equilization, frame by frame
    for cnt = 1:length(strt_data)
        strt = strt_data(1,cnt);            % copy starting index of data
        j = 1;                            % reset for eq. alg.

        % LMS
        for i=strt:strt+pre_lng-1         % iterate through preamble
          rr=rs(i+delta:-1:i-n+1+delta)'; % vector of received signal
          e=s_train(j)-rr'*f;             % calculate error
          err(k) = e;                       % save e for error checking 
          f=f+mu1*e*rr;                     % update equalizer coefficients
          k = k+1;
          pre_eq(j)=rr'*f;
          j = j+1;
        end
        r_eq = horzcat(r_eq,pre_eq);        % combine equalized preambles

        % DD-LMS
        for i=strt+pre_lng:strt+pre_lng+3*data_lng-1      % iterate through user 1,2 and 3
          rr=rs(i:-1:i-n+1)';               % vector of received signal
          user_eq(i) = quantalph(rr'*f,[-3 -1 1 3]);   % quantize given index for error calc
          e=user_eq(i)-rr'*f;               % calculate error
          err(k) = e;                       % save e for error checking
          f=f+mu2*e*rr;                     % update equalizer coefficients
          k = k+1;                          % index error term
        end
        r_eq = horzcat(r_eq,user_eq(strt+pre_lng:strt+pre_lng+3*data_lng-1));               % combine to equalization vector
    end

    if DEBUG == 1
        figure(8);
        subplot(2,1,1), plot(r_eq(2:end),'b.')      % plot constellation diagram of rs
        title('constellation diagram (equilization)');
        subplot(2,1,2), plot(err(2:end))            % plot trajectory of error
    end

    %% Quantization (implement frame by frame quantizing and preamlbe locating)    
    find_mrkr = xcorr(s_train,r_eq);           % locate preamble
    if DEBUG == 1
       figure (9);
       plot(find_mrkr);
       title('Correlation of EQ data and Preamble');
    end
 
    y = quantalph(r_eq(2:end),[-3 -1 1 3]);    % quantize vector
    y = y';

    %% Decodeing
    user1_msg = 0;
    user2_msg = 0;
    user3_msg = 0;

    find_mrkr = xcorr(s_train,y);       % locate preamble
    if DEBUG == 1
        figure (11);
        plot(find_mrkr);
        title('Correlation of Preamble and Quantized sig');
    end
    
    find_mrkr = xcorr(s_train,rs);            % correlate
    if min(find_mrkr) > max(find_mrkr)        % inversion check
        find_mrkr = -find_mrkr;
        y = -y;
    end
    [~,strt_data] = max(abs(find_mrkr));      % locate preamble
    strt_data=mod(strt_data,frame_lng):frame_lng:length(y)-frame_lng;
    strt = strt_data(1,1) + pre_lng;          % set to start of U1 data for frame
       
    for i = 1:length(strt_data)                  % index through frame by frame (broken on last iteration)
        user1 =  pam2letters2(y(strt:strt+data_lng-1));
        strt = strt+data_lng;             % set to start of U2 data for frame
        user2 =  pam2letters2(y(strt:strt+data_lng-1));
        strt = strt+data_lng;             % set to start of U3 data for frame
        user3 =  pam2letters2(y(strt:strt+data_lng-1));
        strt = strt+data_lng+pre_lng;     % point at next U1 data frame

        user1_msg = horzcat(user1_msg,user1);
        user2_msg = horzcat(user2_msg,user2);
        user3_msg = horzcat(user3_msg,user3);
    end

    if user == 1
        disp(user1_msg);
        decoded_msg=user1_msg(2:end);

    elseif user == 2
        disp(user2_msg);
        decoded_msg=user2_msg(2:end);

    elseif user == 3
        disp(user3_msg);
        decoded_msg=user3_msg(2:end);

    end

end