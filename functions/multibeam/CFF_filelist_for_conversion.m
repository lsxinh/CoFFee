function [IN2,OUT2] = CFF_filelist_for_conversion(IN,varargin)
% function [IN2,OUT2] = CFF_filelist_for_conversion(IN,varargin)
%
% DESCRIPTION
%
% function to sort out filenames, prior to .all to .mat conversion. 
%
% INPUT VARIABLES
%
% In all following cases, a specified folder MUST finish with a filesep.
%
% 1. only IN as a string filename.
% -> turn IN as a cell of one string. Write OUT with same foldername and same filename as IN
% ex: [in1,out1] = CFF_filelist_for_conversion('D:\Alex\test\test1.all')
% -> in1 = {'D:\Alex\test\test1.all'}
% -> out1 = {'D:\Alex\test\test1.mat'}
%
% 2. IN as a string filename & varargin{1} is OUT as a string filename
% -> turn IN and OUT as cells of one string
% ex: [in2,out2] = CFF_filelist_for_conversion('D:\Alex\test\test1.all','D:\Alex\test\tagada\test.mat')
% -> in2 = {'D:\Alex\test\test1.all'}
% -> out2 = {'D:\Alex\test\tagada\test.mat'}
%
% 3. only IN as a cell of string filenames
% -> Keep IN as is. Write OUT with same foldername and same filename as IN
% ex: [in3,out3] = CFF_filelist_for_conversion({'D:\Alex\test\test1.all','D:\Alex\test\test2.all'})
% -> in3 = {'D:\Alex\test\test1.all','D:\Alex\test\test2.all'}
% -> out3 = {'D:\Alex\test\test1.mat','D:\Alex\test\test2.mat'}
%
% 4. IN as cell of string filenames & varargin{1} is OUT as a cell of string filenames
% -> Keep IN and OUT as is
% ex: [in4,out4] = CFF_filelist_for_conversion({'D:\Alex\test\test1.all','D:\Alex\test\test2.all'},{'D:\Alex\test\boum\test1.mat','D:\Alex\test\boum\test2.mat'})
% -> in4 = {'D:\Alex\test\test1.all','D:\Alex\test\test2.all'}
% -> out4 = {'D:\Alex\test\boum\test1.mat','D:\Alex\test\boum\test2.mat'}
%
% 5. IN as a cell of string filenames & varargin{1} is OUT as a string folder name
% -> Keep IN as is. Write OUT with requested foldername and same filenames as IN
% ex: [in5,out5] = CFF_filelist_for_conversion({'D:\Alex\test\test1.all','D:\Alex\test\test2.all'},'D:\Alex\test\tsoin\')
% -> in5 = {'D:\Alex\test\test1.all','D:\Alex\test\test2.all'}
% -> out5 = {'D:\Alex\test\tsoin\test1.mat','D:\Alex\test\tsoin\test2.mat'}
%
% 6. IN as a string folder name & varargin{1} is the extension type (starting with a period) as string
% -> parse for files in IN and turn to cells. Write OUT in same foldername and same filenames as IN
% ex: [in6,out6] = CFF_filelist_for_conversion('D:\Alex\test\','.all')
%
% 7. IN as a string folder name & varargin{1} is the extension type as string & varargin(2) is OUT as string folder name
% -> parse for files in IN and turn to cells. Write OUT with requested foldername and same filenames as IN
% ex: [in7,out7] = CFF_filelist_for_conversion('D:\Alex\test\','.all','D:\Alex\test\paf\')
%
% RESEARCH NOTES
%
% EXAMPLES
%
% % all following examples need the said input files to exist:
% [in1,out1] = CFF_filelist_for_conversion('D:\Alex\test\test1.all')
% [in2,out2] = CFF_filelist_for_conversion('D:\Alex\test\test1.all','D:\Alex\test\tagada\test.mat')
% [in3,out3] = CFF_filelist_for_conversion({'D:\Alex\test\test1.all','D:\Alex\test\test2.all'})
% [in4,out4] = CFF_filelist_for_conversion({'D:\Alex\test\test1.all','D:\Alex\test\test2.all'},{'D:\Alex\test\boum\test1.mat','D:\Alex\test\boum\test2.mat'})
% [in5,out5] = CFF_filelist_for_conversion({'D:\Alex\test\test1.all','D:\Alex\test\test2.all'},'D:\Alex\test\tsoin\')
% [in6,out6] = CFF_filelist_for_conversion('D:\Alex\test\','.all')
% [in7,out7] = CFF_filelist_for_conversion('D:\Alex\test\','.all','D:\Alex\test\paf\')
%
% NEW FEATURES
%%% 
% Alex Schimel, Deakin University
% Version 1 (20-02-2014)
%%%




if ischar(IN) && ~exist(IN,'dir') && exist(IN,'file')
    % IN is a string filename. Turn IN as a cell of one string.
    IN2 = {IN};
    if nargin==1
        % Write OUT with same foldername and same filename as IN
        [p,n,e] = fileparts(IN);
        OUT2 = {[p filesep n '_' e(2:end) '.mat']};
    elseif nargin==2 && ischar(varargin{1})
        % varargin{1} is OUT as a string filename. Turn as cell of one
        % string
        OUT2 = {varargin{1}};
    else
        error('Cannot recognize input');
    end
elseif iscellstr(IN)
    % IN is a cell of string filenames. Keep IN as is.
    IN2 = IN;
    if nargin==1
        % Write OUT with same foldername and same filename as IN
        for ii=1:length(IN)
            [p,n,e] = fileparts(IN{ii});
            OUT2{ii} = [p filesep n '_' e(2:end) '.mat'];
        end
    elseif nargin==2 && iscellstr(varargin{1})
        % varargin{1} is OUT as a cell of string filenames. Keep OUT as is
        OUT2 = varargin{1};
    elseif nargin==2 && ischar(varargin{1}) && strcmp(varargin{1}(end),filesep)
        % varargin{1} is OUT as a string folder name (finishing with a
        % filesep). Write OUT with requested foldername and same filenames
        % as IN
        for ii=1:length(IN)
            [p,n,e] = fileparts(IN{ii});
            OUT2{ii} = [varargin{1} n  '_' e(2:end) '.mat'];
        end
    else
        error('Cannot recognize input');
    end
elseif ischar(IN) && exist(IN,'dir') && nargin>1 && ischar(varargin{1}) && strcmp(varargin{1}(1),'.')
    % IN is a string folder name and varargin{1} is extension. Parse for
    % files in IN and turn to cells.
    listing = dir(IN);
    iout=0;
    for ii = 1:length(listing)
        [p,n,e] = fileparts(listing(ii).name);
        if listing(ii).isdir==0 && length(n)>0 && strcmp(e,varargin{1})
            iout=iout+1;
            IN2{iout}=[IN listing(ii).name];
        end
    end
    if nargin==2
        % Write OUT in same foldername and same filenames as IN
        for ii=1:length(IN2)
            [p,n,e] = fileparts(IN2{ii});
            OUT2{ii} = [p filesep n  '_' e(2:end) '.mat'];
        end
    elseif nargin==3 && ischar(varargin{2}) && strcmp(varargin{2}(end),filesep)
        % varargin{2} is OUT as string folder name (finishing with a
        % filesep). Write OUT with requested foldername and same filenames
        % as IN
        for ii=1:length(IN2)
            [p,n,e] = fileparts(IN2{ii});
            OUT2{ii} = [varargin{2} n  '_' e(2:end) '.mat'];
        end
    else
        error('Cannot recognize input');
    end
else
    error('Cannot recognize input');
end


