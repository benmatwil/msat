function convert, gridcoord, x, y, z

    r = gridcoord[0]
    t = gridcoord[1]
    p = gridcoord[2]
    
    r = x[floor(r)] + (r-floor(r))(x[ceil(r)]-x[floor(r)])
    t = y[floor(t)] + (t-floor(t))(y[ceil(t)]-y[floor(t)])
    p = z[floor(p)] + (p-floor(p))(z[ceil(p)]-z[floor(p)])
    
    return, [r*cos(p)*sin(t), r*sin(p)*sin(t), r*cos(t)]

end