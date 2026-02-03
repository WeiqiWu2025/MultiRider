function [WiFi_transceiver_link_snr] = survey_MultiRider_funcChannelQuality2SNR(WiFi_transceiver_link)


low_l = -15;
low_r = 0;
medium_l = 0;
meduim_r = 15;
high_l = 15;
high_r = 30;
extremely_l = 30;
extremely_h = 45;


switch WiFi_transceiver_link
    case 'low quality'
        WiFi_transceiver_link_snr = randi([low_l,low_r]);
    case 'medium quality'
        WiFi_transceiver_link_snr = randi([medium_l,meduim_r]);
    case 'high quality'
        WiFi_transceiver_link_snr = randi([high_l,high_r]);
    case 'extremely high quality'
        WiFi_transceiver_link_snr = randi([extremely_l,extremely_h]);
    otherwise
        error('channel error!');
end

end

