pro map, n, rsphere=rsphere, ball=ball, converge=converge, fan=fan

  if keyword_set(converge) then begin
    spawn, 'make'
    spawn, 'sf_converge n=' + string(n, format='(I4.4)')
  endif
  
  plotstruct, n, rsphere=rsphere, ball=ball, converge=converge, fan=fan
  
end
