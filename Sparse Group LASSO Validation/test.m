% --- beta ---
KL = 0.68096902;
L0 = -0.89643015;
KF = 0.91955807;
tib_KF = -0.70354371;

beta_ = zeros(39,1);
beta_(15) = KL;
beta_(16) = KL*L0;
beta_(26) = KF;
beta_(27) = tib_KF;
% --- X ---
sol_L = X(:, XMap('soleus_r.fiber_length_norm'));
sol_F = X(:, XMap('soleus_r.mtu_force_norm'));
tib_F = X(:, XMap('tib_ant_r.mtu_force_norm'));
% --- sigma ---
sigma = [X, ones(size(X,1),1)] * M * beta_;
% --- plot ---
figure, plot(y), hold on
plot(sigma)
legend('sigma1', 'sigma2')