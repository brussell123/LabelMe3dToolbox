function params = ReadObjectHeights(txtfile)
% Inputs:
% txtfile - Text file name
%
% Outputs:
% params.objNames
% params.mu_obj
% params.sig_obj

params.objNames = [];
params.mu_obj = [];
params.sig_obj = [];
fp = fopen(txtfile);
while 1
  tline1 = fgetl(fp); % Object name
  if ~ischar(tline1), break, end
  tline2 = fgetl(fp); % Height distribution
  if ~ischar(tline2), break, end

  stats = regexp(tline2,' ','split');

  if length(stats)>0
    params.objNames{end+1} = tline1;
    params.mu_obj(end+1) = str2num(stats{1});
  end
  if length(stats)>1
    params.sig_obj(end+1) = str2num(stats{2});
  else
    params.sig_obj(end+1) = 0;
  end
end
fclose(fp);
