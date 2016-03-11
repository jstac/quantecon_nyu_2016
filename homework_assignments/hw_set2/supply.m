function ys = supply(price);

   global b;

   ys = exp(b*price) - 1;

end

