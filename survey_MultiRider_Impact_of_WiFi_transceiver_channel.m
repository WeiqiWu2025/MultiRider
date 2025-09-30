clc;
clear;
close all;
addpath(genpath(pwd));
rng(1993); % For repeatable results

%%%%%*** Waveform Configuration ***%%%%%
% Create a format configuration object for a 1-by-1 HT transmission
cfgHT = wlanHTConfig;
cfgHT.ChannelBandwidth = 'CBW20'; % 20 MHz channel bandwidth
cfgHT.NumTransmitAntennas = 1; % 1 transmit antennas
cfgHT.NumSpaceTimeStreams = 1; % 1 space-time streams
cfgHT.PSDULength = 2000; % PSDU length in bytes % 64
cfgHT.MCS = 0; % 1 spatial streams, BPSK rate-1/2
cfgHT.ChannelCoding = 'BCC'; % BCC channel coding

fs = wlanSampleRate(cfgHT); % Get the baseband sampling rate
ofdmInfo = wlanHTOFDMInfo('HT-Data',cfgHT); % Get the OFDM info
ind = wlanFieldIndices(cfgHT); % Indices for accessing each field within the time-domain packet

%%%%%*** Simulation Parameters ***%%%%%
snr = 40;
global numTags;
numTags = 3;

global seqLenForEstChannel;
seqLenForEstChannel = 20;
refSyms = survey_MultiRider_funcGeneratePreamble(seqLenForEstChannel,numTags);
refSyms(refSyms == 0) = 2;
refSyms(refSyms == 1) = 0;

global len_refSyms;
len_refSyms = size(refSyms,1);

maxNumPackets = 30000; % The maximum number of packets at an SNR point

S = numel(snr); % 返回数组snr中元素的个数
numBitErrs = zeros(S,numTags);
berEst = zeros(S,numTags);

numDataSubcarrier = 52;
numPilotSubcarrier = 4;
estH = ones(numDataSubcarrier+numPilotSubcarrier,numTags);
global interval;
interval = ceil(numDataSubcarrier/numTags);

WiFi_transceiver_link = 'high quality';
packetLoss = zeros(S,1);

for i = 1:S
    disp(['SNR: ',num2str(snr(i)),' dB...']);
    % Set random substream index per iteration to ensure that each
    % iteration uses a repeatable set of random numbers
    stream = RandStream('combRecursive','Seed',0);
    stream.Substream = i;
    RandStream.setGlobalStream(stream);
    
    % Loop to simulate multiple packets
    n = 1; % Index of packet transmitted
    while  n<=maxNumPackets
        disp(['snr: ',num2str(snr(i)),' dB -> ','n: ',num2str(n),'-th packet']);
        %%%%%*** TX side ***%%%%%
        % Generate a packet waveform
        txPSDU = randi([0 1],cfgHT.PSDULength*8,1); % PSDULength in bytes
        tx = wlanWaveformGenerator(txPSDU,cfgHT);
        tx = [tx; zeros(15,cfgHT.NumTransmitAntennas)]; % Add trailing zeros to allow for channel filter delay
        
        exSig = [];
        H_TX_Tags = [];
        %%%%%*** TX-Tags backscatter channel & AWGN
        for chan_tx_tag_idx1 = 1:numTags
            bxCoeffForTxTag_real = -1+(1+1)*rand(1,1);
            bxCoeffForTxTag_imag = -1+(1+1)*rand(1,1);
            bxCoeffForTxTag_real = bxCoeffForTxTag_real*0.1;
            bxCoeffForTxTag_imag = bxCoeffForTxTag_imag*0.1;
            bxCoeffForTxTag = bxCoeffForTxTag_real+1i*bxCoeffForTxTag_imag;
            tmp_exSig = tx.*bxCoeffForTxTag;
            exSig = [exSig,tmp_exSig];
            H_TX_Tags = [H_TX_Tags,bxCoeffForTxTag];
        end
        
        
        %%%%%*** Tags side ***%%%%%
        % Backscatter at the tag
        temp = ceil((cfgHT.PSDULength*8+16+6)/26);
        numSymForPsdu = 0;
        numSymForTailPad = 0;
        if mod(temp,2) == 1
            numSymForPsdu = (numel(tx)-720-15-80-80-80)/80;
            numSymForTailPad = 2;
        else
            numSymForPsdu = (numel(tx)-720-15-80-80)/80;
            numSymForTailPad = 1;
        end
        numTagData = numSymForPsdu; % modulate one tag data per one symbol
        
        % Initial tags data
        tagData = zeros(numTagData,numTags);
        numPayload = numTagData-len_refSyms;
        actualPayloadBits = zeros(numPayload,numTags);
        for tag_idx1 = 1:numTags
            payload = randi([0,1],numPayload,1);
            actualPayloadBits(:,tag_idx1) = payload;
            tagData(:,tag_idx1) = [refSyms(:,tag_idx1);payload];
        end
        
        % backscatter modulation
        for tag_idx2 = 1:numTags
            bxSig{tag_idx2} = survey_MultiRider_funcBackscatter(exSig(:,tag_idx2),tagData(:,tag_idx2),1);
        end
        
        %%%%%***** Backscatter channel & AWGN ***%%%%%
        H_Tags_RX = [];
        for chan_tag_rx_idx1 = 1:numTags
            bxCoeffForTagRx_real = -1+(1+1)*rand(1,1);
            bxCoeffForTagRx_imag = -1+(1+1)*rand(1,1);
            bxCoeffForTagRx_real = bxCoeffForTagRx_real*0.01;
            bxCoeffForTagRx_imag = bxCoeffForTagRx_imag*0.01;
            bxCoeffForTagRx = bxCoeffForTagRx_real+1i*bxCoeffForTagRx_imag;
            bxSig{chan_tag_rx_idx1} = bxSig{chan_tag_rx_idx1}.*bxCoeffForTagRx; % backscatter channel
            bxSig{chan_tag_rx_idx1} = awgn(bxSig{chan_tag_rx_idx1},snr(i),'measured');
            H_Tags_RX = [H_Tags_RX,bxCoeffForTagRx];
        end
        actualH = H_TX_Tags.*H_Tags_RX;
        
        %%%%%*** WiFi TX to WiFi RX channel & AWGN ***%%%%%
        coeffForTxRx_real = -1+(1+1)*rand(1,1);
        coeffForTxRx_imag = -1+(1+1)*rand(1,1);
        coeffForTxRx_real = coeffForTxRx_real*0.1;
        coeffForTxRx_imag = coeffForTxRx_imag*0.1;
        coeffForTxRx = coeffForTxRx_real+1i*coeffForTxRx_imag;
        rxSig_from_WiFi_TX = coeffForTxRx.*tx;
        
        WiFi_transceiver_link_snr = survey_MultiRider_funcChannelQuality2SNR(WiFi_transceiver_link);
      
        
        for rx_idx1 = 1:numTags
            ofdmDataSymbols{rx_idx1} = survey_MultiRider_funcGetFrequencyDomainSymbols(bxSig{rx_idx1}(ind.HTData(1):ind.HTData(2)),cfgHT);
            channelEstimation{rx_idx1} = ofdmDataSymbols{rx_idx1}(:,1+(1+(rx_idx1-1)*seqLenForEstChannel:rx_idx1*seqLenForEstChannel));
            ofdmCarriedData{rx_idx1} = ofdmDataSymbols{rx_idx1}(:,(1+len_refSyms)+(1:numPayload));
        end
        
        % 信号错位叠加
        rxSig = zeros(size(ofdmCarriedData{1},1),size(ofdmCarriedData{1},2));
        for rx_idx2 = 1:numTags
            tmp_rx_ref_l1 = (rx_idx2-1)*interval+1;
            tmp_rx_ref_r1 = numDataSubcarrier;
            tmp_rx_l1 = 1;
            tmp_rx_r1 = numDataSubcarrier-tmp_rx_ref_l1+1;
            rxSig(tmp_rx_ref_l1:tmp_rx_ref_r1,:) = rxSig(tmp_rx_ref_l1:tmp_rx_ref_r1,:) + ofdmCarriedData{rx_idx2}(tmp_rx_l1:tmp_rx_r1,:);
        end
        
        [rx_From_WiFi_TX,~,~] = func_awgn(rxSig_from_WiFi_TX,WiFi_transceiver_link_snr,'measured');
        [decodedPSDU,flag] = survey_Mecha_funcWiFiRX(rx_From_WiFi_TX,cfgHT);
        if flag == 0
            packetLoss(i,1) = packetLoss(i,1) + 1;
            for rx_idx3 = 1:numTags
                numBitErrs(i,rx_idx3) = numBitErrs(i,rx_idx3) + numPayload;
            end
            n = n+1;
            continue;
        end
        
        [~,ofdmSymDerived] = survey_MultiRider_funcOFDMSymDerived(decodedPSDU,cfgHT);
        [cfgOFDM,dataInd,pilotInd] = wlan.internal.wlanGetOFDMConfig(cfgHT.ChannelBandwidth, cfgHT.GuardInterval, 'HT', cfgHT.NumSpaceTimeStreams);
        ofdmDataDerived = ofdmSymDerived(cfgOFDM.DataIndices,:,:);
        ofdmPilotsDerived = ofdmSymDerived(cfgOFDM.PilotIndices,:,:);
        
        ofdm_ref_data = ofdmDataDerived(:,(1+len_refSyms)+(1:numPayload));
        ofdm_ref_channelestimation = ofdmDataDerived(:,1+(1:len_refSyms));
        
        % Calculate channel coefficients
        for rx_idx4 = 1:numTags
            A = channelEstimation{rx_idx4};
            B = ofdm_ref_channelestimation(:,1+(rx_idx4-1)*seqLenForEstChannel:rx_idx4*seqLenForEstChannel);
            tmp_LL = size(A,1);
            for rx_idx5 = 1:tmp_LL
                tmp_H_real = funcLSEstimator(B(rx_idx5,:)',real(A(rx_idx5,:))');
                tmp_H_imag = funcLSEstimator(B(rx_idx5,:)',imag(A(rx_idx5,:))');
                tmp_H = tmp_H_real + 1i*tmp_H_imag;
                chanEst{rx_idx4}(rx_idx5,1) = tmp_H;
            end
        end
%         chanEst = ones(52,1);
        
        for rx_idx6 = 1:numTags
            tmp_rx_l2 = (rx_idx6-1)*interval+1;
            tmp_rx_r2 = tmp_rx_l2+interval-1;
            tmp_rx_r2 = min(tmp_rx_r2,numDataSubcarrier);
            tmp_distance = tmp_rx_r2 - tmp_rx_l2 + 1;
            subcarriers_ref_tag = ofdm_ref_data(1:tmp_distance,:);
            subcarriers_received_tag = rxSig(tmp_rx_l2:tmp_rx_r2,:);
            demodPayloadBits(:,rx_idx6) = survey_MultiRider_funcDemodulation(chanEst{rx_idx6}(1:tmp_distance),subcarriers_ref_tag,subcarriers_received_tag);
            % 重构干扰信号
            re_construct_interference_signal = survey_MultiRider_funcReConstructInterferenceSignal(demodPayloadBits(:,rx_idx6),ofdm_ref_data,chanEst{rx_idx6});
            % 消去干扰信号
            tmp_ref_l = (rx_idx6-1)*interval+1;
            tmp_ref_r = numDataSubcarrier;
            tmp_l = 1;
            tmp_r = numDataSubcarrier-tmp_ref_l+1;
            rxSig(tmp_ref_l:tmp_ref_r,:) = rxSig(tmp_ref_l:tmp_ref_r,:) - re_construct_interference_signal(tmp_l:tmp_r,:);
        end
        
        
        % calculate the number of bits
        for comm_idx1 = 1:numTags
            numBitErrs(i,comm_idx1) = numBitErrs(i,comm_idx1) + biterr(actualPayloadBits(:,comm_idx1),demodPayloadBits(:,comm_idx1));
        end
        n = n+1;
        
    end
    % calculate bit error rate
    for comm_idx2 = 1:numTags
        berEst(i,comm_idx2) = numBitErrs(i,comm_idx2)/(numPayload*maxNumPackets);
    end
    
end

aaa = 1;


