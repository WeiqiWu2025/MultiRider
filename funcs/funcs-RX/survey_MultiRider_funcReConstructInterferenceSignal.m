function re_construct_interference_signal = survey_MultiRider_funcReConstructInterferenceSignal(tagData,ofdm_ref_data,chanEst)

len = length(tagData);
re_construct_interference_signal = ofdm_ref_data;
for idx1 = 1:len
    if tagData(idx1) == 1
        re_construct_interference_signal(:,idx1) = re_construct_interference_signal(:,idx1).*(-1);
    end
end
re_construct_interference_signal = re_construct_interference_signal.*chanEst;

end

