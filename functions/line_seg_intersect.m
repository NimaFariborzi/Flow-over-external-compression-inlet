function [ix, iy] = line_seg_intersect(x1,y1,x2,y2,x3,y3,x4,y4)
    dx1 = x1 - x2;
    dx2 = x3 - x4;
    dy1 = y1 - y2;
    dy2 = y3 - y4;
    den = dx1*dy2 - dy1*dx2;
    ix = ((x1*y2 - y1*x2)*dx2 - dx1*(x3*y4 - y3*x4)) / den;
    iy = ((x1*y2 - y1*x2)*dy2 - dy1*(x3*y4 - y3*x4)) / den;
end