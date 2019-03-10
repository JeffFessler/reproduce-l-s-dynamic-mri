function [Samp] = strucrand(n1,n2,n3,line);

 for frameno = 1:n3,
     
%    
      y=-n1/2:1:n1/2-0.5;  % Create an array of N points between -n1/2 and n1/2
      x=(linspace(-n2/2,n2/2,length(y)));
      i=1;
     
      inc = 0 + (pi /line)*rand;
    
     % Loop to traverse the kspace ; 0 to pi in steps of pi/lines -- 
     % Succesive frames rotated by a small random angle (pi/line)*rand (see
     % above)
   
     for ang=0:pi/line:pi-1e-3;
         klocn=complex(y*cos(ang+ (inc)),x*sin(ang+ (inc)));
         kloc_all(:,i)=klocn;
         i=i+1;
     end
     
     
     % Round the collected data to the nearest cartesian location   
     kcart=round(kloc_all+(0.5+0.5*1i));
    % plot(kcart,'*');title('k locations after nearest neighbor interpolation: Center (0,0)');
%     
%     
    % Next, shift the cartesian locations accordingly such that the center
    % is now at (n1/2,n1/2); {Previously the center in kcart was (0,0)}
    kloc1 = round(kcart)+((n1/2+1)+(n2/2+1)*1i);
    kloc1real = real(kloc1); kloc1real = kloc1real - n1*(kloc1real>n1);
    kloc1imag = imag(kloc1); kloc1imag = kloc1imag - n2*(kloc1imag>n2);
    kloc1real = kloc1real + n1*(kloc1real<1);
    kloc1imag = kloc1imag + n2*(kloc1imag<1);
    kloc1 = kloc1real + 1i*kloc1imag;
    
 
  %
%   Create the sampling pattern
    for i=1:size(kloc1,1)
        for j=1:size(kloc1,2)
     
        Samp(real(kloc1(i,j)),imag(kloc1(i,j)),frameno) = 1;
        
        end
    end

 end
         
end
