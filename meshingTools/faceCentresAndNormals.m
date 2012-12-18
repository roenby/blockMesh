function [Cf,Sf] = faceCentresAndNormals(p,f)

%Initial guess of face centres
x = zeros(size(f,1),3);
for k = 1:size(f,2)
    x = x + p(f(:,k),:);
end
x = x/size(f,2);

%Calculating surface normal and face centre by triangulation with
%initially guessed face centre
Sf = zeros(size(f,1),3);
Cf = zeros(size(f,1),3);
sumAbsSf = zeros(size(f,1),3);
for k = 1:size(f,2)
    x1 = p(f(:,k),:);
    kPlus1 = mod(k,size(f,2))+1;
    x2 = p(f(:,kPlus1),:);
    Sk = 1/2*cross(x2-x1,x-x1);
    Sf = Sf + Sk;
    Cf = Cf + 1/3*norm(Sk).*(x1 + x2 + x);
    sumAbsSf = sumAbsSf + norm(Sk); 
end
Cf = Cf./sumAbsSf;