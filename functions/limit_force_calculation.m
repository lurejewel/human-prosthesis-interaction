function force = limit_force_calculation(ang, vel, angUpTh, angLowTh)
kUp = 2;
kLow = 2;
damp = 0.2;

kUp_ = kUp * ( atan(10*(ang-angUpTh))/pi + 0.5 );
kLow_ = kLow * ( atan(10*(ang-angLowTh))/pi + 0.5 );
fUpLim = kUp_ .* ( angUpTh - ang );
fLowLim = kLow_ .* ( angLowTh - ang );
fDamp = -damp * (kUp_/kUp + kLow_/kLow_) .* vel;
force = fUpLim + fLowLim + fDamp;

end