function cfm_open_manual()
%find the location of the asmita app and open the manual
appinfo = matlab.apputil.getInstalledAppInfo;
idx = find(strcmp({appinfo.name},'ChannelForm'));
fpath = [appinfo(idx(1)).location,[filesep,'ChannelForm',filesep,'doc',...
                                    filesep,'ChannelForm manual.pdf']];
open(fpath)
