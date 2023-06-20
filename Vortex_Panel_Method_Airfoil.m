clc; clear 

 

alpha_value = linspace(-2,14,(17)); 

cl_value = zeros(1,length(alpha_value)); 

cm_c_four_value = zeros(1,length(alpha_value)); 

xcp_value = zeros(1,length(alpha_value)); 

 

for i_final = 1:length(alpha_value) 

alpha = alpha_value(i_final); 

 

V = 1; 

num = 8; 

 

%% Body input section, unit radius 

% panel = [1, 2, 3, 4 ,5 ... 

 

m = 0.02; % max camber 

p=0.4; % psotion of max camber 

t = 0.12; 

k1 = 15.957;%??? 

 

r = 1.1019.*(t^2); %radius of leading edge circle 

x1 = linspace(r/3,p,num); %x coordinates nose cicle to m 

x2 = linspace(p,1,num); %x coordinates m to 1 

x = ([x1 x2(:,[2:num])]); 

 

y_cam_1 = (m/p^2).*(2*p.*x1-(x1.^2)); %camber line forward of max ordinate 

y_cam_2 = (m/((1-p)^2)).*((1-2*p)+(2*p.*x2)-x2.^2); %camber line afte of maximum ordinate 

y_cam = [y_cam_1 y_cam_2(:,[2:num])]; 

 

dy_cam_1 = (m/p^2).*(2*p-(2.*x1)); %derivative forward 

dy_cam_2 = (m/((1-p)^2)).*(2*p-(2.*x2)); %derivative of aft 

dy_cam = [dy_cam_1 dy_cam_2(:,[2:num])]; %merged derivative of camber line 

theta = atan(dy_cam); %slope of camber line 

 

y_t = 5.*t.*((0.29690.*sqrt(x))-(0.12600.*x)-(0.35160.*(x.^2)) +(0.28430.*(x.^3))-(0.10150.*(x.^4))); %thickness equation 

 

x_upper = x-(y_t.*sin(theta)); %x coordinates of upper surface 

x_lower = x+(y_t.*sin(theta)); %x coordinates of lower surface 

y_upper = y_cam+(y_t.*cos(theta)); %y coordinates of upper surface 

y_lower = y_cam-(y_t.*cos(theta)); %y coordinates of lower surface 

 

x_lower = flip(x_lower); 

y_lower = flip (y_lower); 

 

X = [x_lower x_upper]; 

Y = [y_lower y_upper]; 

X(1) = 1; 

Y(1) = 0; 

X(end) = 1; 

Y(end) = 0; 

 

 

% vortex panel length, mid-point and orientation angle 

n = length(X)- 1; 

 

Sj = zeros(n,1); phi = zeros(n,1); 

xmp = zeros(n,1); ymp = zeros(n,1); 

for i = 1:n 

    % Length of each panel 

    Sj(i) = sqrt((X(i+1) - X(i))^2 + (Y(i+1) - Y(i))^2); 

    % mid-point of each panel 

    xmp(i) = 0.5*(X(i+1) + X(i)); 

    ymp(i) = 0.5*(Y(i+1) + Y(i)); 

    phi(i) = atan2d( (Y(i+1) - Y(i)),((X(i+1) - X(i))) ); 

end 

 

%  

% subplot(2,1,1);hold on; 

% plot(X,Y,'o-'); plot(xmp,ymp,'*'); label_1 = string(1:n); 

% text(xmp,ymp,label_1,'VerticalAlignment','bottom');  

% grid on; hold off; 

 

 

%% build the Cn(i,j) matrix   

Cn2 = eye(n,n)*1;  Cn1 = eye(n,n)*-1;Ct1 = eye(n,n)*pi/2; 

Ct2 = eye(n,n)*pi/2; 

for i = 1:n 

    for j = 1:n 

        if i ~= j 

     

            A = -(xmp(i) - X(j))*cosd(phi(j))-(ymp(i) - Y(j))*sind(phi(j));

            C = sind(phi(i)-phi(j)); 

            B = (xmp(i) - X(j))^2 + (ymp(i) - Y(j))^2; 

            D = cosd(phi(i) - phi(j)); 

            E = (xmp(i)-X(j))*sind(phi(j))-(ymp(i)-Y(j))*cosd(phi(j)); 

            F = log(1 + ((Sj(j))^2 + (2*A*Sj(j))) / B); 

            G = atan2((E*Sj(j)) , (B + A*Sj(j))); 

            P = (xmp(i)-X(j))*sind(phi(i)-2*phi(j)) + (ymp(i)-Y(j))*cosd(phi(i)-2*phi(j)); 

            Q = ((xmp(i) - X(j)) * cosd(phi(i) - 2*phi(j))) - ((ymp(i) - Y(j)) * sind(phi(i) - 2*phi(j))); 

            % normal component 

            Cn2(i,j) = D + ((0.5*Q*F)/Sj(j))-((A*C + D*E)*(G/Sj(j))); 

            Cn1(i,j) = 0.5*D*F + C*G - Cn2(i,j); 

            Ct2(i,j) = C + ((0.5*P*F)/Sj(j))+((A*D - C*E)*(G/Sj(j))); 

            Ct1(i,j) = 0.5*C*F - D*G - Ct2(i,j); 

 

        end 

    end 

end 

 

An = zeros(n+1,n+1); RHS = zeros(n+1,1); 

for i= 1:n 

    An(i,1) = Cn1(i,1); 

    An(i,n+1) = Cn2(i,n); 

    RHS(i) = sind(phi(i)-alpha); 

    for j = 2:n 

        An(i,j) = Cn1(i,j) + Cn2(i,j-1); 

    end 

end 

An(n+1,1) = 1; 

An(n+1,n+1) = 1;  

 

for j = 2:n 

        An(n+1,j) = 0; 

end 

RHS(n+1) = 0; 

 

At = zeros(n+1,n+1); 

for i=1:n 

    At(i,1) = Ct1(i,1); 

    At(i,n+1) = Ct2(i,n); 

    for j =2:n 

        At(i,j) = Ct1(i,j) + Ct2(i,j-1); 

    end 

   

 

end  

gamma_prime = An\RHS; 

 

V = zeros(n+1,1); 

cp = zeros(n,1); 

 

for i = 1:n  

    V(i) = cosd(phi(i)-alpha); 

    for j = 1:n+1 

    V(i) = V(i) + At(i,j)*gamma_prime(j); 

    cp(i) = 1 - (V(i))^2; 

    end 

   

end 

% subplot(2,1,2); hold on; 

% plot(xmp,cp); 

% set(gca,'Ydir','reverse'); 

% %text(theta_p1,Cp_p1,label_1,'VerticalAlignment','bottom'); 

% grid on; hold off 

 

 

xmp_lower = xmp(1:(length(xmp))/2); 

cp_lower = flip(cp(1:(length(cp))/2)); 

ymp_lower = flip (ymp(1:(length(xmp))/2)); 

 

xmp_upper = xmp(((length(xmp))/2)+1:end); 

ymp_upper = ymp(((length(xmp))/2)+1:end); 

 

cp_upper = cp(((length(xmp))/2)+1:end); 

 

dcp = cp_lower-cp_upper; 

 

dy_upper = zeros(length(ymp_upper)-1,1); 

dy_lower = zeros(length(ymp_upper)-1,1); 

dx_upper = zeros(length(ymp_upper)-1,1); 

dx_lower = zeros(length(ymp_upper)-1,1); 

 

 

 

for i=1:(length(ymp_upper)-1)  

    dy_upper(i) = ymp_upper(i+1)-ymp_upper(i); 

    dy_lower(i) = ymp_lower(i+1)-ymp_lower(i); 

    dx_upper(i) = xmp_upper(i+1)-xmp_upper(i); 

    dx_lower(i) = xmp_lower(i+1)-xmp_lower(i);  

end 

dyu_dx = dy_upper./dx_upper; 

dyl_dx = dy_lower./dx_upper; 

 

cd_cp_upper = cp_upper; 

cd_cp_lower = cp_lower; 

cd_xmp_upper = xmp_upper; 

 

cd_cp_upper(end)=[]; 

cd_cp_lower(end) = [];  

cd_xmp_upper (end) =[]; 

 

 

cd_integral_value = (cd_cp_upper.*dyu_dx) - (cd_cp_lower.*dyl_dx); 

cd = trapz(cd_xmp_upper,cd_integral_value); 

cl = trapz(xmp_upper,dcp); 

cmle = trapz(xmp_upper,-dcp.*xmp_upper); 

xcp = -cmle/cl; 

cm_c_four = cmle +(cl/4); 

 

cl_value(i_final) = cl;  

cm_c_four_value(i_final) = cm_c_four;  

xcp_value(i_final) = xcp; 

end 

alpha_value = transpose(alpha_value); 

cl_value = transpose(cl_value); 

cm_c_four_value= transpose(cm_c_four_value); 

 

 

 

 

subplot(3,1,1); hold on; 

plot(alpha_value,cl_value); 

xlabel("angle of attack") 

ylabel("Section Lift Coeffecient") 

 

grid on; hold off 

 

 

 

 

 

subplot(3,1,2); hold on; 

plot(alpha_value,cm_c_four_value); 

xlabel("angle of attack") 

ylabel("Moment Coeffecient") 

%text(theta_p1,Cp_p1,label_1,'VerticalAlignment','bottom');grid on; hold off 

grid on; hold off 

 

subplot(3,1,3); hold on; 

plot(alpha_value,xcp_value); 

xlabel("angle of attack") 

ylabel("xcp") 

%text(theta_p1,Cp_p1,label_1,'VerticalAlignment','bottom'); 

grid on; hold off 

 

a_zero_fit = fit(alpha_value,cl_value,'poly1'); 

p=coeffvalues(a_zero_fit); 

a_zero = p(1); 

cm_slope = fit(alpha_value,cm_c_four_value,'poly1'); 

p=coeffvalues(cm_slope); 

m_zero = p(1); 

 

x_ac = (-m_zero/a_zero)+0.25; 

cm_ac = cl_value(3)*(x_ac-0.25); 