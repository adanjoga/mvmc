function [fun,opt] = cviconfig(cvi)
% Type of clustering criteria
if strcmpi(cvi,'sil')
    fun = @silindex; % Silhoutte index 1
    opt = 'mx'; 
elseif strcmpi(cvi,'wgs')
    fun = @wgsindex;
    opt = 'mn';   
else
    error('Unknown clustering criteria');
end