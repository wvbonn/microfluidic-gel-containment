clear all
close all

% All units in c/g/s
% The ambient pressure is set at zero
% The injection is assumed to occur at a low enough velocity that the
% contact angle of the moving pre-gel is its static contact angle

% some geometry parameters
w_gap = 40e-4; % post-gap width
w_gel = 0.25e-1; % gel region width
h = 300e-4; % device height (nominal)
n_gap = 25; % number of post gaps
a_post = 100e-4; % wide base of posts
d_post = 40e-4; % depth of posts (along y)

% equivalent hydraulic radius
equi_r = 0.65*(h*(w_gel-2*d_post))^0.625/(h+w_gel-2*d_post)^0.25;

% physical constants
gamma = 60; % surface tension coefficient collagen/air (static, Kezwon & Wojciechowski (2014) Colloids & Surfaces A)
theta_c = 2*pi/3; % static contact angle water/hydrophobic PDMS/air (assumed same as water/hydrophobic glass/air...)
mu_inf = 11e-2;
sigma_y = 120e-2; % Casson fit for 2 mg/ml collagen precursor at 4°C (Gudapati et al. (2020) Table SI 2)

% some code parameters
n_front = ceil(n_gap/2); % post-gap number at which the front of the pre-gel is assumed to be
v = 0.1; % injection velocity

% vector of post-gap positions (in x)
x_gap = 0:(a_post+w_gap):(a_post+w_gap)*n_gap;

% Pressure at the front of the pre-gel during injection 
r_front(1) = (w_gel-2*d_post)/2/cos(pi/2-theta_c);
r_front(2) = h/2/cos(pi/2-theta_c);
p_front = gamma*(1/r_front(1)+1/r_front(2));

% viscosity
shr = v*2/h; % max shear rate
mu = (sqrt(mu_inf)+sqrt(sigma_y/shr))^2; % Casson viscosity (assumed that the viscosity is uniform in the device, not true for gel near air interfaces in y) 
p_gel = zeros(1,n_gap);
p_gel(n_front+1:n_gap) = 0;
p_gel(n_front) = p_front;
p_gel(1:n_front-1) = p_front + 8*mu*v*(n_front-1:-1:1)*(a_post+w_gap)*h*(w_gel-2*d_post)/pi/equi_r^4;

% at the post-gap interfaces, radii of curvature and contact angles with
% fluid channel surfaces
r_gap_z = NaN*ones(1,n_gap); % yz plane
r_gap_x = NaN*ones(1,n_gap); % xy plane
d_gap = NaN*ones(1,n_gap); % max extension of interface into fluid channel (i.e where the radii are calculated)
theta_xy = NaN*ones(1,n_gap); % contact angle of interface in xy plane
theta_yz = NaN*ones(1,n_gap); % same in yz plane
for i = 1:n_front
    r_gap_z(i) = 1/2*(lsqnonlin(@(x)2*gamma*(1/x+1/(x/2-sqrt(x^2-h^2)/2+w_gap^2/(2*x-2*sqrt(x^2-h^2))))-p_gel(i),5*h,h,[]));
    d_gap(i) = r_gap_z(i) - sqrt(r_gap_z(i)^2-h^2/4);
    r_gap_x(i) = d_gap(i)/2 + w_gap^2/8/d_gap(i);
    theta_xy(i) = pi-acos(sign(d_gap(i)-w_gap/2)*w_gap/2/r_gap_x(i));
    theta_yz(i) = pi-acos(h/2/r_gap_z(i));
end

% plot the Young-Laplace pressure difference
% z = h:0.0001:20*h;
% figure;plot(z,2*gamma*(1./z+1./(z/2-sqrt(z.^2-h^2)/2+w_gap^2./(2*z-2*sqrt(z.^2-h^2))))-p_gel(1));


figure('position',[50 0 800 800],'color','w');
subplot(2,2,1);hold all;box on;grid on;
set(gca,'fontsize',13);
xlabel('post number from inlet'); xlim([0 n_gap+1]);
ylabel('radius of curvature (z, \mum)');
plot(1e4*r_gap_z,'o','markers',10,'linew',1.5);
plot(0);
line([0 n_gap+1],1e4*[h/2 h/2],'linew',1.5);
Y = get(gca,'YLim');
text(n_gap-3,1e4*h/2+0.05*diff(Y),'half-height','fontsize',13);

subplot(2,2,2);hold all;box on;grid on;
set(gca,'fontsize',13);
xlabel('post number from inlet'); xlim([0 n_gap+1]);
ylabel('radius of curvature (x, \mum)');
plot(1e4*r_gap_x,'+','markers',10,'linew',1.5);
plot(0);
line([0 n_gap+1],1e4*[w_gap/2 w_gap/2],'linew',1.5);
Y = get(gca,'YLim');
text(n_gap-3,1e4*w_gap/2+0.05*diff(Y),'half-width','fontsize',13);

subplot(2,2,3);hold all;box on;grid on;
set(gca,'fontsize',13);
xlabel('post number from inlet'); xlim([0 n_gap+1]);
ylabel('interface extension (y, \mum)');
plot(1e4*d_gap,'o','markers',10,'linew',1.5);

subplot(2,2,4);hold all;box on;grid on;
set(gca,'fontsize',13);
xlabel('post number from inlet'); xlim([0 n_gap+1]);
ylabel('interface meniscus angle (°)'); 
plot(180/pi*theta_xy,'o','markers',10,'linew',1.5);
plot(180/pi*theta_yz,'+','markers',10,'linew',1.5);
line([0 n_gap+1],180/pi*[theta_c theta_c],'linew',1.5);
Y = get(gca,'YLim');
text(1,180/pi*theta_c+0.05*diff(Y),'static contact angle','fontsize',13);
legend('xy plane','yz plane');


