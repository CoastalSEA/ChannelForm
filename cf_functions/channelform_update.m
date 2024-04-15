function channelform_update(obj,oldV,newV)
%
%-------header-------------------------------------------------------------
% NAME
%   channelform_update.m 
% PURPOSE
%   update saved models to newer versions of Asmita
% USAGE
%   channelform_update(oldV,newV) 
% INPUTS
%   obj - instance of model
%   oldV - old version number as a character string
%   newV - new version number as a character string
% RESULTS
%   saved model updated to new version. If this is called from Asmita this
%   will be the version that is being run at the time.
% NOTES
%   Called in muiModelUI.loadModel when old and new version numbers do not
%   match.
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2023 
%--------------------------------------------------------------------------
%
    if strcmp(newV,'3.3')
        if strcmp(oldV,'2.1')
            update_v21_to_v31(obj);
        end
        update_v32_to_v33(obj);
    else
        warndlg(sprintf('No update for version %s to version %s', oldV,newV))
    end
end

