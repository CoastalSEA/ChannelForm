function open_manual()
%find the location of the asmita app and open the manual
appinfo = matlab.apputil.getInstalledAppInfo;
idx = find(strcmp({appinfo.name},'ChannelForm'));
fpath = [toolboxes(idx(1)).location,'/doc/ChannelForm manual.pdf'];
open(fpath)
