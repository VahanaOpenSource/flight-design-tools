function t=plotAirfoil(filename)

%Read airfoil file
fid=fopen(filename);
tline=fgetl(fid); tline=fgetl(fid);
j=1; x=[];
while ischar(tline)
    if ~isempty(tline) && sum(isstrprop(tline,'alpha'))<2 && ~any(strfind(tline,'...'))
        x(j,:)=str2num(tline); j=j+1;
    end
    tline=fgetl(fid);
end
fclose(fid);

%Calculate thickness if requested
if nargout==1
    xu=x([-1;diff(x(:,1))]<=0,:);
    xl=x([-1;diff(x(:,1))]>0,:);
    
    t=0;
    for j=1:length(xu)
        t=max(t,xu(j,2)-interp1(xl(:,1),xl(:,2),xu(j,1)));
    end
    t=t/max(x(:,1));
end

%Plot arifoil
plot(x(:,1),x(:,2))
xlabel('x');ylabel('y');
axis equal