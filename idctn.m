
function a = idctn(a)

    number_of_dim=ndims(a);

    transpose=0;
    if (number_of_dim==2) && (size(a,1)==1)
        transpose=1; a=a';
    end

    siz=size(a);
    number_of_dim=ndims(a);
    
    for i=1:number_of_dim,
        n=siz(i); clear tmp;
        
        ww{i} = 2*exp((-1i*pi/(2*n))*(0:n-1)')/sqrt(2*n);
        ww{i}(1) = ww{i}(1)/sqrt(2);
        
        tmp(1:2:n)=(1:ceil(n/2)); 
        tmp(2:2:n)=(n:-1:ceil(n/2)+1);
        ind{i}=bsxfun(@plus, tmp', 0:n:n*(prod(siz)/n-1));
        if (siz(i)==1), break; end;
    end
    if number_of_dim==2
        if (min(siz)==1),
            a=idct(a,ww{1},ind{1});if transpose, a=a'; end;  % 1D case
        else
            a = idct(idct(a,ww{1},ind{1}).',ww{2},ind{2}).'; % 2D case
        end       
    else
        
    for i=1:number_of_dim
        a=reshape(idct(reshape(a,siz(1),[]),ww{i},ind{i}), siz); %3D case
        siz=[siz(2:end) siz(1)];                  
        a=shiftdim(a,1);                           
    end
   end


function a = idct(a,ww,ind)

    
    a = bsxfun(@times,  ww, a); 
    a = fft(a);                 % fft
    a=real(a(ind));              
    
    


