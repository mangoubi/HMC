function G = Hessian(xx,AA,UU)

  G = AA;
       
        Q = xx'*UU;
        
        logist = ((1+(exp(1)).^(-Q)).^(-2));
        
        logist = logist.*((exp(1)).^(-Q));
        
         
        
       G=G+UU*diag(logist)*UU';
        
end
