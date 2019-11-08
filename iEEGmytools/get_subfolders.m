
function dfolders = get_subfolders(path)
d = dir(path)
% remove all files (isdir property is 0)
dfolders = d([d(:).isdir]==1)
% remove '.' and '..'
dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}))

end