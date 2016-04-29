pro map, n, rsphere=rsphere, ball=ball, converge=converge, fan=fan

  spawn, 'make'
  spawn, 'sf_converge n=' + string(n, format='(I4.4)')
  
  plotstruct, n, rsphere=rsphere, ball=ball, converge=converge, fan=fan
  
end
