global a b epsilon;
a = 1;
b = 0.1;
epsilon = 1;

mxiter = 30;
toler = 1.0e-6;

plow = 0.1;
phigh = 10.0;

niter = mxiter;

for i = 1:mxiter;

  pcur = (plow + phigh)/2;

  yd = demand(pcur);
  ys = supply(pcur);

  excesssupply = ys - yd;

  if (excesssupply > 0); 
     phigh = pcur;
  else;
     plow = pcur;
  end;

  diff = abs(phigh - plow);

  if (diff <= toler);
     niter = i;
     break;
   end;

end;

pclear = (plow + phigh)/2;
yd = demand(pcur);
ys = supply(pcur);
excesssupply = ys - yd;

[niter pclear yd ys excesssupply]
