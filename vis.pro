pro vis,noerase=noerase

if ~keyword_set(noerase) then begin
device,decomposed=0
window,xs=800,ys=800


nx=0l
ny=0l
nz=0l
openr,10,'data/newmag.dat'
readu,10,nx,ny,nz
close,10


surface,dist(100),xrange=[0,nx],yrange=[0,ny],zrange=[0,nz],/nodata,/save,az=-30,ax=30,xtitle='x',ytitle='y',ztitle='z'

endif

nnulls=0l
openr,10,'nulls.dat',/f77_unformatted
readu,10,nnulls
close,10

iloop=255

for nulls=1,nnulls do begin
;iloop=nulls*127

n=0l
openr,10,'ringidl00'+strtrim(nulls,2)+'.dat',/f77_unformatted
readu,10,n
while n gt 0 do begin
  x=dblarr(n)
  y=dblarr(n)
  z=dblarr(n)
  readu,10,x,y,z
  plots,x,y,z,/t3d,psym=3
  readu,10,n
endwhile
close,10


;iloop=0
;if nulls eq 1 then begin
openr,11,'sep00'+strtrim(nulls,2)+'.dat',/f77_unformatted
readu,11,n
while n gt 0 do begin
 ;   iloop=iloop+85
    print, n
        
	x=dblarr(n)
	y=dblarr(n)
	z=dblarr(n)
	
	readu,11,n
	readu,11,x,y,z
	

	loadct,13,/silent
	plots,x,y,z,/t3d,psym=2,symsize=1,color=iloop
	loadct,0,/silent
	
	readu,11,n
	
endwhile
close,11

;endif

endfor
 
stop

end

