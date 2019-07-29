function z=Benchmark(x , iter)
% FDA1
    
%     f1=x(1);
%     
%     nt = 10;
%     toet = 10;
%     
%     t = ( 1 / nt ) * floor(iter / toet);
%     Gt = sin(0.5 * pi * t);
%     g = 1 + sum(( x(2: end) - Gt ) .^ 2);
%     h = 1 - sqrt(f1/g);
%     f2 = g * h;
%     z=[f1 f2]';

    a=0.8;
    
    b=3;
    
    z1=sum(-10*exp(-0.2*sqrt(x(1:end-1).^2+x(2:end).^2)));
    
    z2=sum(abs(x).^a+5*(sin(x)).^b);
    
    z=[z1 z2]';


end