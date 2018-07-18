window,0,xsize=800,ysize=400
readcol,'fitslist',files,format='A'

n=n_elements(files)

lamint=findgen(800)*6+3600

darr=fltarr(800,404)

for i=0,n-1 do begin
	d=mrdfits(files[i],0,h)
	c0=fxpar(h,'COEFF0')
	c1=fxpar(h,'COEFF1')
	z=fxpar(h,'Z')
	npix=n_elements(d[*,0])
	lam=10.0^(c0+c1*findgen(npix))
	arr=interpol(smooth(d[*,1],3),lam/(1+z),lamint)
	for j=0,3 do darr[*,(i*4)+j]=arr
endfor
tvscl,(-10) > (darr)< 30
end