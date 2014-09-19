function testread

A=read_array('out/combined.dat');

B=fftshift(fft2(fftshift(A)));
figure; imagesc(abs(B)); colormap('gray');

end

function A=read_array(fname)

F=fopen(fname,'rb');

fseek(F,0,'bof');
RO=read_readout(F);
N=length(RO);
M=1;
while (~feof(F))
	X=read_readout(F);
	if (size(X,1)>1)
		M=M+1;
	end;
end;

disp([N,M]);

A=zeros(N,M);
fseek(F,0,'bof');
for j=1:M
	X=read_readout(F);
	A(:,j)=X;
end;

fclose(F);


end

function X=read_readout(F)

X=0;
N=fread(F,1,'uint32');
fread(F,31,'uint32');
if (~feof(F))
	X=fread(F,512,'float32');
	X=X(1:2:end)+i*X(2:2:end);
end;

end