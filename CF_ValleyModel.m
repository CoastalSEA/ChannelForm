classdef CF_ValleyModel < GDinterface  
%
%-------class help---------------------------------------------------------
% NAME
%   CF_ValleyModel.m
% PURPOSE
%   Class for valley form model used in ChannelForm App
% NOTES
%   Co-ordinate convention is that the x-axis is seaward from the head or
%   tidal limit. The y-axis is the cross-valley axis mirrored about the
%   centre-line.
% SEE ALSO
%   muiDataSet and GDinterface
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%     
    properties
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties:
    end
    
    methods (Access = private)
        function obj = CF_ValleyModel()                   
            %class constructor
        end
    end      
%%
    methods (Static)        
%--------------------------------------------------------------------------
% Model implementation
%--------------------------------------------------------------------------         
        function obj = runModel(mobj)
            %function to run a simple 2D diffusion model
            obj = CF_ValleyModel;             
            %now check that the input data has been entered
            %isValidModel checks the InputHandles defined in ModelUI
            if ~isValidModel(mobj, metaclass(obj).Name)  
                warndlg('Use Setup to define model input parameters');
                return;
            end          
%--------------------------------------------------------------------------
% Model code
%--------------------------------------------------------------------------
%             wlvobj = getClassObj(mobj,'Inputs','WaterLevels');
            hydobj = getClassObj(mobj,'Inputs','CF_HydroData');
            setTransHydroProps(hydobj,mobj); %initialise transient properties
            
            %assign the run parameters to the model instance           
            setRunParam(obj,mobj); 
            %add water level definition to the run parameters
            setWaterLevels(obj);

            %call valley model
            [xi,yi,zi,yz,Lv,Ls] = cf_valley_model(obj);
            meta.data = sprintf('Valley with Lv=%d and Ls=%d',Lv,Ls);
            if isempty(xi), return; end
            
            griddata = reshape(zi,1,length(xi),length(yi));  
            %x and y are 1xn and yz is {[1xn]x3}
            if size(yz,2)~=3, yz = yz'; end
            results = [{griddata},yz];
%             xydata = {xi,yi};     %assign xyz coordinates as row vectors
            %now assign results to object properties  
            gridyear = years(0);  %durataion data for rows 
            dims = struct('x',xi,'y',yi,'t',gridyear);
%--------------------------------------------------------------------------
% Assign model output to dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------                   
            %assign metadata about model and save grid
            meta.source = metaclass(obj).Name;
            dst = setGridOuput(obj,results,dims,meta);
%--------------------------------------------------------------------------
% Add property dstables in function GDinterface.setFormProperties
%--------------------------------------------------------------------------  
            dst = setFormProps(obj,dst,1,meta);            
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------             
            setDataSetRecord(obj,mobj.Cases,dst,'form_model');
            getdialog('Run complete');
            DrawMap(mobj);
        end
%%
        function addForm2Valley(mobj)
            %merge a selected form with the selected valley form 
            obj = getClassObj(mobj,'Cases','CF_ValleyModel');
            if isempty(obj) 
                warndlg('No valley model available');
                return; 
            end
            ftxt = 'Select Form Model to use:';
            [fgrid,frmobj] = getGrid(obj,mobj,1,{'CF_FormModel'},ftxt);
            xf = fgrid.x;            
            if min(xf)<0, xf = max(xf)-xf; end
            vtxt = 'Select Valley Model to use:';
            [vgrid,valobj] = getGrid(obj,mobj,1,{'CF_ValleyModel'},vtxt);
            
            [X,Y] = ndgrid(xf,fgrid.y);
            zv = griddata(vgrid.x,vgrid.y,vgrid.z',X,Y);    %valley elevations
            new_z = max(fgrid.z,zv);
            [m,n] = size(new_z);
            new_z = reshape(new_z,1,m,n);
            
%             
%             [f_useCase,v_useCase] = selectForm(obj,mobj,promptxt);
%             if isempty(f_useCase) || isempty(v_useCase), return; end
%             [new_z,frmobj] = Form2Valley(obj,mobj,v_useCase,f_useCase);
            
            
            answer = questdlg('Add or update existing?','Add Valley','Add','Update','Add');
            
            if strcmp(answer,'Add')
                %create new record
%                 lobj = mobj.(h_f)(id_c);
%                 irec = length(lobj.(fprop{1}))+1;
%                 newds = lobj.(fprop{1}){id_rec};
%                 newds.zLevel = new_z;
%                 
%                 lobj.ModelData{irec} = newds;                      %assign form to class property
%                 lobj.GrossProps{irec} = lobj.GrossProps{id_rec} ;  %asign gross property results
%                 lobj.Hypsommetry{irec} = lobj.Hypsommetry{id_rec}; %assign hypsommetry results
%                 lobj.PrismArea{irec} = lobj.PrismArea{id_rec};     %assign area-prism results
%                 if isa(lobj,'HydroFormModel') || isa(lobj,'CKFAmodel')
%                     lobj.xWaterLevels{irec}  = lobj.xWaterLevels{id_rec};%assign water levels
%                 end
%                 lobj.cidValley(irec) = mobj.Cases.CaseID(v_useCase);
%                 Results.saveResults(mobj,h_f,lobj);
            else            
                %overwrite exisitng form data set with new form             
                frmobj.Data.Form.Z(1,:,:) = new_z;  
                frmobj.MetaData.valleyID = caseID(mobj.Cases,valobj.CaseIndex);
            end            
            getdialog('Data modified');
        end 
%%
        function new_z = updateValley(F,F0,zhw,trp,trg,incFP)
            %add the valley form to the channel form 
%             fp_offset = 0.1; 

            [X,Y] = ndgrid(F.xi,F.yi);
            zv = griddata(F.xv,F.yv,F.zv',X,Y);    %valley elevations
            
            new_z = max(F.zi,zv);

%             idx = F.zi>trp.Zhw & zv<=trp.Zhw;  %areas that are undefined
%             F.zi(F.zi>trp.Zhw) = 0;
%             zv(zv<=trp.Zhw) = 0;
%             new_z = F.zi+zv;
%             new_z(idx) = trp.newHW(idx);
            
            if ~trg.isConstrained && ~incFP
                fact = 4; %needs to be consistent with value in TransgressionModel.applyConstraints
                idz = Y<(F.Yhw+fact*trp.dy) & Y>(F.Yhw);
                Zhw = repmat(zhw,1,length(F.yi));
                new_z(idz) = Zhw(idz)+trp.offset;  %add offset from high water to flood plain surface
                %reset area outside to original levels
                idz = Y>(F.Yhw+fact*trp.dy);
                new_z(idz) = F0.zi(idz);
            elseif trg.isConstrained
                idz = Y<(F.Yhw+trp.dy) & Y>(F.Yhw);
                new_z(idz) = trp.wcrest;
                %reset area outside to original levels
                idz = Y>(F.Yhw+trp.dy);
                new_z(idz) = F0.zi(idz);
            end       
        end
%%
        function restoreForm(mobj)
            %remove valley data from form to restore base form
            if isempty(mobj.Cases), return; end

            %idf = strcmp(mobj.Cases.CaseType,'_model');  
            idf = contains(mobj.Cases.CaseType,'_model');
            if isempty(idf)
                return;
            elseif sum(idf)>1
                [f_useCase,~,ok] = ScenarioList(mobj.Cases,'_model',...
                                                    'ListSize',[200,140]);
                if ok<1, return; end
            else
                f_useCase = find(idf);
            end
            
            robj = mobj.Cases;
            f_caseid = mobj.Cases.CaseID(f_useCase);             
            [h_f,id_c,fprop,id_rec,aprop] = getCaseRecord(robj,mobj,f_caseid);  
            localObj = mobj.(h_f)(id_c);
            zf = localObj.(fprop{1}){id_rec}.zLevel;            
            rnp = mobj.(getClassHandle(mobj,'CFRunProps'));
            
            if isa(localObj,'ChannelFormModel') || isa(localObj,'PRFormModel')
                hyd = mobj.(getClassHandle(mobj,'HydroFormData'));
                zhw = hyd.MTLatMouth+hyd.TidalAmplitude;
                zf(zf>zhw) = zhw+2*rnp.histint;
            else
                cst_res = localObj.(aprop{4}){id_rec};
                zhw = cst_res{1}+cst_res{2};
                zf = squeeze(zf);
                [m,n] = size(zf);
                for i=1:m
                    idx = zf(i,:)>zhw(i); %areas that are undefined
                    zf(i,idx) = zhw(i)+2*rnp.histint;
                end
                zf = reshape(zf,1,m,n);
            end     

            %overwrite exisitng form data set with new form
            localObj.(fprop{1}){id_rec}.zLevel = zf;
            localObj.cidValley = [];
            ModelUI.myDialog('Data modified');
        end 
%%
        function componentsPlot(mobj)
            %generate a plot of the channel, valley and combined form
            obj = getClassObj(mobj,'Cases','CF_ValleyModel');
            if isempty(obj) 
                warndlg('No valley model available');
                return; 
            end
            ftxt = 'Select Form Model to use:';
            [fgrid,frmobj] = getGrid(obj,mobj,1,{'CF_FormModel'},ftxt);
            xf = fgrid.x/1000;  %change x-y axes to km
            if min(xf)<0, xf = max(xf)-xf; end
            yf = fgrid.y/1000;
            zf = grid.z;
            
            vtxt = 'Select Valley Model to use:';
            [vgrid,valobj] = getGrid(obj,mobj,1,{'CF_ValleyModel'},vtxt);
            zv = vgrid.z;
            
            [X,Y] = ndgrid(xf,fgrid.y);
            zv = griddata(vgrid.x,vgrid.y,vgrid.z',X,Y);    %valley elevations
            new_z = max(fgrid.z,zv);
            [m,n] = size(new_z);
            zc = reshape(new_z,1,m,n);
            
            
            
            
            
            defaultvals = {num2str(floor(min(zf,[],'All'))),'12'};
            prompt = {'Lower z-limit','Upper z-limit'};
            inp = inputdlg(prompt,'Zlimits',1,defaultvals);
            zlimits = [str2double(inp{1}),str2double(inp{2})];
                        
            hf = figure('Name','Form Components','Tag','PlotFig');
            ax = axes('Parent',hf,'Tag','SurfacePlot');
            
            s1 = subplot(1,3,1,ax);
            contouredSurfacePlot(obj,s1,xf,yf,zf',zlimits,'(a) Channel Form');
            colorbar('off');
            zlim(zlimits)
            s1.Position = [0.05,0.1,0.25,0.8]; 

            s2 = subplot(1,3,2);
            contouredSurfacePlot(obj,s2,xf,yf,zv',zlimits,'(b) Valley Form');
            colorbar('off');
            zlim(zlimits)
            s2.Position = [0.35,0.1,0.25,0.8];
            
            s3 = subplot(1,3,3);
            contouredSurfacePlot(obj,s3,xf,yf,zc',zlimits,'(c) Combined Form');
            zlim(zlimits)
            s3.Position = [0.65,0.1,0.25,0.8];
            
            channel = frmobj.Data.Form.Description;
            valley = valobj.Data.Form.Description;
            ht = sgtitle(sprintf('Form using %s and %s',channel,valley));
            ht.FontSize = 12;        
        end        
        
    end
%%
    methods
        function tabPlot(obj,src) %abstract class method for muiDataSet
            %generate plot for display on Q-Plot tab
            cf_model_tabs(obj,src);
        end
    end 
%%    
    methods (Access = private) 

       
%%
        function setWaterLevels(obj)
            %set the water levels at the mouth as constant values along
            %channel
            grdobj = obj.RunParam.GD_GridProps;
            hydobj = obj.RunParam.CF_HydroData;
            Lt = diff(grdobj.XaxisLimits);   %length of model domain (m)
            nintx = grdobj.Xint;             %no of intervals in the x direction
            xi = 0:Lt/nintx:Lt;
            obj.RunParam.CF_HydroData.zhw = ones(size(xi))*hydobj.zhw; %high water
            obj.RunParam.CF_HydroData.zmt = ones(size(xi))*hydobj.zmt; %mean tide level
            obj.RunParam.CF_HydroData.zlw = ones(size(xi))*hydobj.zlw; %low water
            obj.RunParam.CF_HydroData.cstres = [];
        end
    end
end