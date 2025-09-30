function demd_tag_data = survey_MultiRider_funcDemodulation(chanEst,subcarriers_ref_tag,subcarriers_received_tag)

len = length(subcarriers_ref_tag);
demd_tag_data = zeros(len,1);

for idx1 = 1:len
    y = subcarriers_received_tag(:,idx1);
    x = subcarriers_ref_tag(:,idx1);
    value1 = norm(y-(-1.*chanEst).'.*x);
    value0 = norm(y-(1.*chanEst).'.*x);
    if value1>value0
        demd_tag_data(idx1,1) = 0;
    else
        demd_tag_data(idx1,1) = 1;
    end
end

end

