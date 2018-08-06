function alphainterp = ThrottleInterp(t,AlphaList,tsearch)

alphainterp = interp1(t,AlphaList,tsearch);
% alphainterp = spline(t,AlphaList,tsearch);
end