function a = dctn(a)
number_of_dim=ndims(a);
   

transpose=0;
if (number_of_dim==2) && (size(a,1)==1)
    transpose=1; a=a';
end

    siz=size(a);
    number_of_dim=ndims(a);
    
    for i=1:number_of_dim
        n=siz(i);
        
        ww{i} = 2*exp(((-1i*pi)/(2*n))*(0:n-1)')/sqrt(2*n);
        ww{i}(1) = ww{i}(1) / sqrt(2);
        ind{i}=bsxfun(@plus, [(1:2:n) fliplr(2:2:n)]', 0:n:n*(prod(siz)/n-1));
        if (siz(i)==1), break; end
    end
    

if number_of_dim==2
    if (min(siz)==1)
        a=dct(a,ww{1},ind{1});if transpose, a=a'; end  % 1D case
    else a = dct(dct(a,ww{1},ind{1}).',ww{2},ind{2}).'; % 2D case
    end       
else
    
    for i=1:number_of_dim
        a=reshape(dct(reshape(a,siz(1),[]),ww{i},ind{i}), siz); %3D case
        siz=[siz(2:end) siz(1)];                   
        a=shiftdim(a,1);                        
    end
end


function a=dct(a,ww,ind)

    

     a=a(ind);     
     a = fft(a);  % ifft
     a = real(bsxfun(@times,  ww, a)); 
   
   

      

