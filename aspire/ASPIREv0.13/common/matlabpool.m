%
% F.Poitevin - September 2017
% The following was adapted to recent Matlab versions,
% where matlabpool is replaced with parpool.
function varargout = matlabpool(varargin)
%
%function matlabpool(varargin)
%     varargout={};
     arg1=varargin{1};
     switch arg1
       case 'open'
%         parpool(varargin{2:end});
	 parpool('local', varargin{2:end});  
       case 'close'
         delete(gcp('nocreate'));
       case  'size'
         p = gcp('nocreate'); % If no pool, do not create new one.
         if isempty(p)
            poolsize = 0;
         else
            poolsize = p.NumWorkers;
         end
         varargout={poolsize};
    otherwise %number of workers
         if ~isnumeric(arg1)
           arg1=num2str(arg1);
         end
         parpool(arg1);
end
% - additional code ends here.
