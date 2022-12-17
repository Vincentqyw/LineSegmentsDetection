function[sfilenames] = check_save_requirements(save_flags,outputdir,etype);

    name1       = ['edge',etype];
    sfilenames  = {name1,'blur.mat','dark.mat','light.mat',...
                   'g1mag.mat','g1dir.mat','g1scale.mat','g2mag.mat','g2scale.mat',...
                   'xzero1.mat','yzero1.mat','xzero2.mat','yzero2.mat'};
    
	if any(save_flags)
        
        % Output directory does not exist:
        if ((exist(outputdir) ~= 7) & (length(outputdir) ~= 0))
            mkdir(outputdir);        % Create this directory    
            mess1   = 'Requested output directory has been created.';
            uiwait(msgbox(mess1,'Output Directory Creation','none','modal'));
        end;
        
        if (length(outputdir) == 0)
            outputdir   = 'tmp_output';
            opmessage   = ['Output directory name not specified.  ',...
                           'Files will be saved in directory "tmp_output".'];
            uiwait(msgbox(opmessage,'Output Directory Warning','warn','modal'));
            if (exist(outputdir) ~= 7)
                mkdir(outputdir);
            end;
        end;
            
        filestosave = sfilenames(find(save_flags));
        
        x   = dir(outputdir);
        xx  = {x.name}';
        commonfiles = intersect(filestosave,xx);
    
        if ~isempty(commonfiles)
            maxnumx = 0;
            for k = 1:1:length(xx)
                y = xx{k};
                numx = 0;
                if (length(y)>4)
                    if (y(end-4)=='x')
                        numx = sum(y=='x');
                    end;
                end;
                if (numx > maxnumx)
                    maxnumx = numx;
                end;
            end;
            assignin('base','maxnumx',maxnumx+1);
            assignin('base','commonfiles',commonfiles);
  
%             uiwait(outputdir_erase_info);  % Call GUI
%         
%             changename_flag = evalin('base','changename_flag;');
%         
%             if changename_flag
%                repx = repmat('x',1,maxnumx+1);
%                for k = 1:1:length(sfilenames)
%                    s    = sfilenames{k};
%                    s    = [s(1:end-4),repx,s(end-3:end)];
%                    sfilenames{k} = s;
%                end;
%             end;
           
            evalin('base','clear changename_flag maxnumx commonfiles;');    
        end;    
    end;
    
    for k = 1:1:length(sfilenames)
        s   = sfilenames{k};
        s   = [outputdir,'/',s];
        sfilenames{k} = s;
    end;
    
return;