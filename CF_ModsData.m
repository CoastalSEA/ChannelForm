classdef CF_ModsData < muiPropertyUI               
%
%-------class help---------------------------------------------------------
% NAME
%   CF_ModsData.m
% PURPOSE
%   Class for morphological modification parameters used in ChannelForm model
% USAGE
%   obj = CF_ModsData.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'Distance from mouth to start of modification (m)',...
                          'Distance from mouth to end of modification (m)',...
                          'Distance from centre-line to left side (m +/-0)',...
                          'Distance from centre-line to right side (m +/-0)',...                          
                          'Elevation of modification (mOD)'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        ModStart      %Distance from mouth to start of modification (m)
        ModEnd        %Distance from mouth to end of modification (m)
        ModLeft       %Distance from centre-line to right side (m +/-0)
        ModRight      %Distance from centre line to left side (m +/-0)
        ModElev       %Elevation of modification (mOD) 
    end    

%%   
    methods (Access=protected)
        function obj = CF_ModsData(mobj)       
            %constructor code:            
            %TabDisplay values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);  %muiPropertyUI function
            
            %to use non-numeric entries then one can either pre-assign 
            %the values in the class properties defintion, above, or 
            %specify the PropertyType as a cell array here in the class 
            %constructor, e.g.:
            % obj.PropertyType = [{'datetime','string','logical'},...
            %                                       repmat({'double'},1,8)];
        end 
    end
%%  
    methods (Static)  
        function obj = setInput(mobj,editflag)
            %gui for user to set Parameter Input values
            classname = 'CF_ModsData';          
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = CF_ModsData(mobj);       
            end
            
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                %add nrec to limit length of props UI (default=12)
                obj = editProperties(obj);  
                nrec = length(obj.ModStart);
                if length(obj.ModEnd)~=nrec || ...
                        length(obj.ModRight)~=nrec || ...
                        length(obj.ModLeft)~=nrec || ...
                        length(obj.ModElev)~=nrec
                        warndlg('Check Input. Number of entries must be the same for all variables',...
                            'Morphological Modifications');
                        return;
                end
                checkInput(obj,obj.ModStart,obj.ModEnd);
                checkInput(obj,obj.ModLeft,obj.ModRight); 
            end
            setClassObj(mobj,'Inputs',classname,obj);
        end     
%%
        function addMorphMods(mobj)
            %add a prismatic dredge channel to a channel form
            if isempty(mobj.Cases), return; end
            
            idf = contains(mobj.Cases.CaseType,'model');            
            if isempty(idf)
                return;
            elseif sum(idf)>1
                [useCase,~,ok] = ScenarioList(mobj.Cases,'_model',...
                  'PromptText','Select Channel Form','ListSize',[200,140]);
                if ok<1, return; end
            else
                useCase = find(idf);
            end
            robj = mobj.Cases;
            caseid = mobj.Cases.CaseID(useCase); 
            [h_f,id_f,prop,id_p] = getCaseRecord(robj,mobj,caseid);    
            handle = mobj.(h_f)(id_f);
            seldata = handle.(prop{1}){id_p};
            x = seldata.Properties.UserData.XYZ{1};
            y = seldata.Properties.UserData.XYZ{2};  
            z = squeeze(seldata.zLevel);          
            if min(x)<0, x = max(x)-x; end            
            [X,Y] = ndgrid(x,y);
            
            %get channel dimensions
            inp = mobj.(getClassHandle(mobj,'MorphMods'));
            for i=1:length(inp.ModStart) 
                idx = X>inp.ModStart(i) & X<inp.ModEnd(i);
                idy = Y>inp.ModLeft(i) & Y<inp.ModRight(i);
                z(idx & idy) = inp.ModElev(i);
            end
            %overwrite exisitng form data set with new form   
            [m,n] = size(z);
            new_z = reshape(z,1,m,n);            
            handle.(prop{1}){id_p}.zLevel = new_z;            
            ModelUI.myDialog('Data modified');
        end    
    end
%%  
    methods
        function ok = checkInput(~,var1,var2)
            %check that inputs are correctly matched
            BC = [var1',var2'];
            [~,idx] = max(BC,[],2);
            if any(idx<2)
                warndlg('Check Input. Requires End>Start',...
                    'Morphological Modifications');
                ok = 0;
            else
                ok = 1;
            end
        end
    end  
end