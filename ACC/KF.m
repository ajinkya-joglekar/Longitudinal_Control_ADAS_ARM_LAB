
% Kalman Variables Initiation

%Previous states initialization
yk_1 = 0;
vk_1 = 0;
% System state matrix intialization
sigma_x = 0;
sigma_v = 0;
%Process variables
Pk = 0.0;    %Initial process covariance
G = 0.0;     %Initial Kalman Gain
zk = 0.0;    %Initial measured value of system states
Yk = [0;0];  %Initial corrected value of system states

function Yk = kalman(dist,vel, measurement_noise)

dt = 0.1; % TImestep increment between each measurement

%Control Paramters and state matrices
ak = 0 ;
A = [1 dt ; 0 1] ;
B = [0 ; dt];
H = [1 0 ; 0 1];

% State predictions
Yk_1 = A*[yk_1 ; vk_1] + B*ak ;

%State measurement
zk = H*Yk_1;

%Process covariance matrix
Pk = A*[(sigma_x*sigma_x) (sigma_x*sigma_v) ; (sigma_v*sigma_x) (sigma_v*sigma_v)] * transpose(A);

%Kalman Gain
K = Pk*transpose(H) * inv(H*Pk*transpose(H));

%Update estimate
Yk = Yk_1 + K(zk - H*Yk_1);

%Update error covariance
Pk = (1-K*H)*Pk_;
end
