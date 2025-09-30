function [refSyms] = survey_MultiRider_funcGenerateReferenceSymbols(numTags)
% 0 and 1 --> reflecting; 2 --> absorbing
% 0 -> not flip
% 1 -> flip \pi
% 2 -> absorbing

refSyms = [];

for idx1 = 1:numTags
    start_2_count = 2*(idx1-1);
    remaining_2_count = 2*(numTags-idx1);
    if start_2_count > 0
        start_part = repmat(2,start_2_count,1);
    else
        start_part = [];
    end
    
    middle_part = [0;1];
    
    if remaining_2_count > 0
        end_part = repmat(2,remaining_2_count,1);
    else
        end_part = [];
    end
    
    tmp_refSyms = [start_part;middle_part;end_part];
    refSyms = [refSyms,tmp_refSyms];
    
end


end

