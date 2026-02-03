function [bxSig] = survey_MultiRider_funcBackscatter(exSig,tagData,isLeader)
% 0 and 1 --> reflecting; 2 --> absorbing
% 0 -> not flip
% 1 -> flip \pi
% 2 -> absorbing

global len_refSyms;

bxSig = exSig;
numTagData = length(tagData);
pulseLen = 80;
pfo = comm.PhaseFrequencyOffset('PhaseOffset',180);
for idx1 = 1:numTagData
    if idx1 <= len_refSyms
        switch tagData(idx1)
            case 1
                bxSig(801+(idx1-1)*pulseLen:800+idx1*pulseLen) = pfo(bxSig(801+(idx1-1)*pulseLen:800+idx1*pulseLen));
            case 2
                bxSig(801+(idx1-1)*pulseLen:800+idx1*pulseLen) = 0.*bxSig(801+(idx1-1)*pulseLen:800+idx1*pulseLen);
            case 0
                bxSig(801+(idx1-1)*pulseLen:800+idx1*pulseLen) = bxSig(801+(idx1-1)*pulseLen:800+idx1*pulseLen);
        end
    else
        if tagData(idx1) == 1
            bxSig(801+(idx1-1)*pulseLen:800+idx1*pulseLen) = pfo(bxSig(801+(idx1-1)*pulseLen:800+idx1*pulseLen));
        end
    end
end



end

