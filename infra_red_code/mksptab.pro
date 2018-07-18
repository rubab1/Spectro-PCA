path='/Users/ds/data/SSGSS/'
outpath='R-spect-new/'

m=mrdfits(path+'measured.ssgss.fits',1)

openw,1,outpath+'SSGSS_info-R.txt'
printf,1,'num   ra   dec   z'
for i=0,100 do begin
	printf,1,i+1,m[i].sdss_ra,m[i].sdss_dec,m[i].sdss_z,format='(I8,F12.6,F12.6,F12.4)'
endfor
close,1

readcol,'fitslist',files, format='A'


nspect=n_elements(files)

outnm=strmid(files,7,9)

for i=0,nspect-1 do begin

	d=mrdfits(path+files[i],1)
	openw,1,outpath+'SSGSS_'+strcompress(i+1,/remove_all)+'-R.txt'
	printf,1,'wave flux error channelflag bitflag stitchflag'
	for j=0,n_elements(d.wave)-1 do begin
		printf,1,d.wave[j],d.flux[j],d.error[j],d.channelflag[j],d.bitflag[j],d.stitchflag[j],format='(F9.3,D14.7,D14.7,I12,I12,I12)'
	endfor
	close,1

endfor


end

