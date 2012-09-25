function PI = getPlanes(plane)

pix = str2num(char({plane(:).pix}));
piy = str2num(char({plane(:).piy}));
piz = str2num(char({plane(:).piz}));
piw = str2num(char({plane(:).piw}));
PI = [pix(:) piy(:) piz(:) piw(:)]';
