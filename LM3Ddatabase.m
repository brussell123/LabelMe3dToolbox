function DB = LM3Ddatabase(HOMEANNOTATIONS,HOMEIMAGES,folderlist,filelist)

if nargin > 2
  % Prepend "users/labelme3d" to folderlist:
  for i = 1:length(folderlist)
    folderlist{i} = fullfile('users/labelme3d',folderlist{i});
  end
end

switch nargin
 case 1
  DB = LMdatabase(HOMEANNOTATIONS);
 case 2
  DB = LMdatabase(HOMEANNOTATIONS,HOMEIMAGES);
 case 3
  DB = LMdatabase(HOMEANNOTATIONS,HOMEIMAGES,folderlist);
 case 4
  DB = LMdatabase(HOMEANNOTATIONS,HOMEIMAGES,folderlist,filelist);
 otherwise
  error('Invalid number of arguments.');
end
