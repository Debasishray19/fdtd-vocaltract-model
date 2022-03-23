% This function implements the two-mass vocal model algorithm to compute
% the glottal volume velocity.
% Calculate, sourceVelocity = ug/Area

function [vf_flowParam] = vf_TwoMassModel(SRATE, airParam, vf_structuralParam, vf_flowParam)
    mu  = airParam.mu;
    rho =  airParam.rho;
    
    % Vocal fold mass
    m1 = vf_structuralParam.m1;
    m2 = vf_structuralParam.m2;
    
    % Thickness of each mass
    d1 = vf_structuralParam.d1;
    d2 = vf_structuralParam.d2;
    
    % Linear spring stiffness
    k1 = vf_structuralParam.k1;
    k2 = vf_structuralParam.k2;
    kc = vf_structuralParam.kc;
    
    % Non-linear spring stiffness
    etak1 = vf_structuralParam.etak1;
    etak2 = vf_structuralParam.etak2;

    % Linear stiffness during vocal fold collision
    h1 = vf_structuralParam.h1; 
    h2 = vf_structuralParam.h2;
    
    % Non-linear stiffness during vocal fold collision - Don't delete this
    etah1 = vf_structuralParam.etah1;
    etah2 = vf_structuralParam.etah2;

    % Cross sectional area of glottal slit at rest
    Ag01 = vf_structuralParam.Ag01;
    Ag02 = vf_structuralParam.Ag02;
    
    % Vocal cord effective length
    lg = vf_structuralParam.lg;
    
    % Calculate viscous resitance
    r1_open = vf_structuralParam.r1_open;
    r2_open = vf_structuralParam.r2_open;
    r1_close = vf_structuralParam.r1_close;
    r2_close = vf_structuralParam.r2_close;
    
    % Input area to the vocal tract
    A1 = vf_structuralParam.A1;
    
    % Sub-glottal(Lungs) pressure
    ps = vf_flowParam.ps; % In Pascal(8cmH2O) ==== 1cmH2O = 98.0665 => 784.532

    % Pressure downstream of glottal expansion
    p1 = vf_flowParam.p1;
    
    % Flow velocity change
    ug_old = vf_flowParam.ug_old;
    ug_curr = vf_flowParam.ug_curr;

    % Latteral displacement of vocal fold masses
    x1_old = vf_flowParam.x1_old;
    x1_curr = vf_flowParam.x1_curr;

    x2_old = vf_flowParam.x2_old;
    x2_curr = vf_flowParam.x2_curr;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Two Mass Model Simulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculate the current glottal areas
    Ag1 = Ag01 + 2 * lg * x1_curr;
    Ag2 = Ag02 + 2 * lg * x2_curr;
    
	% Use Equation 18 from IF72 to compute the force
	% Use Equation 17 from IF72 to compute the following parameters: Lg1, Lg2, Rv1, Rv2
	% Implement the conditional table of Equation 14 from IF72 (p-1244)

    % Calculate the parameters (Required to compute f_mass1 & f_mass2 ): Rv1, Lg1, Rv2, Lg2
    Rv1 = (12*mu*lg*lg*d1) / (Ag1*Ag1*Ag1);
    Lg1 = (rho*d1) / (Ag1);
    
    Rv2 = (12*mu*lg*lg*d2) / (Ag2*Ag2*Ag2);
    Lg2 = (rho*d2) / (Ag2);
    
    % COMPUTE THE FORCE F1 FOR THE MASS M1
    if (x1_curr > -Ag01 / (2 * lg) && x2_curr > -Ag02 / (2 * lg))
			p_mass1 = ps - 1.37*(rho / 2)*(ug_curr / Ag1)*(ug_curr / Ag1)...
                      - 0.5*(Rv1*ug_curr + Lg1 * SRATE*(ug_curr - ug_old));
                  
			f_mass1 = p_mass1 * d1*lg; % Force = Pressure* Area;
            
    else
        p_mass1 = ps;
        f_mass1 = p_mass1 * d1*lg; % Force = Pressure* Area
    end
    
    % COMPUTE THE FORCE F2 FOR THE MASS M2
    if (x1_curr > -Ag01 / (2 * lg))   
        if (x2_curr > -Ag02 / (2 * lg))
            p_mass2 = p_mass1 - (0.5*(Rv1 + Rv2)*ug_curr + (Lg1 + Lg2)*(ug_curr - ug_old)*SRATE)...
                      - (rho / 2 * ug_curr*ug_curr)*(1 / (Ag2*Ag2) - 1 / (Ag1*Ag1));
                  
            f_mass2 = p_mass2 * d2*lg; % Forec = Pressure* Area;
            
        else
            p_mass2 = ps;
            f_mass2 = p_mass2 * d2*lg; % Forec = Pressure* Area
        end      
    else
        p_mass2 = 0;
        f_mass2 = p_mass2 * d2*lg; %Forec = Pressure* Area
    end       
    
    % GLOTTAL MASS MOVEMENT
  
    % Resultance force on each mass : Sum of the forces arise due to the below conditions
    %  1. Deflection from the equilibrium position (Linear & Non-Linear stiffness of the vocal cord)
    %  2. Contact force when both vocal cords collide (Deformation of each mass)  
    % Consider these new parameters (h=Linear stiffness during VF collison & r = damping coefficient)
    % Collison Condition: xi + Ag0i/(2*lg)<=0 then h exist otherwise h=0 (p1237 from IF72)
    
    % For mass1
    if (x1_curr + Ag01 / (2 * lg) < 0)        
        h1_upd = h1;
        r1_upd = r1_close;
    else
        h1_upd = 0;
        r1_upd = r1_open;
    end
    
    % For mass2
    if (x2_curr + Ag02 / (2 * lg) < 0)
        h2_upd = h2;
        r2_upd = r2_close;
    else
        h2_upd = 0;
        r2_upd = r2_open;
    end
    
    % To calculate the matrix elements use (A-9) from S87 : a11, a12, a21, a22
	% These matrix elemets will be used for furthur vocal cord displacements: x1_next, x2_next for mass m1 and m2
	% Follow (A-10) solution for the code implementation
	a11 = (k1 + h1_upd + kc) / (SRATE*SRATE) + r1_upd / SRATE + m1;
	a12 = -kc / (SRATE*SRATE);
	a21 = a12; 
	a22 = (k2 + h2_upd + kc) / (SRATE*SRATE) + r2_upd / SRATE + m2;
    
    % Follow (A-14) from S87
    s1_prime = k1*etak1*x1_curr*x1_curr*x1_curr + h1_upd*(Ag01 / (2 * lg) + etah1*(Ag01 / (2* lg) + x1_curr)*(Ag01 / (2 * lg) + x1_curr)*(Ag01 / (2 * lg) + x1_curr)); %KEES
	s2_prime = k2*etak2*x2_curr*x2_curr*x2_curr + h2_upd*(Ag02 / (2 * lg) + etah2*(Ag02 / (2* lg) + x2_curr)*(Ag02 / (2 * lg) + x2_curr)*(Ag02 / (2 * lg) + x2_curr)); %KEES
    
%     if (x1_curr > -Ag01 / (2 * lg))
%         s1_prime = k1 * etak1*x1_curr*x1_curr*x1_curr;
%     else
%         s1_prime = (k1 * etak1*x1_curr*x1_curr*x1_curr) - h1_upd * (Ag01 / (2 * lg) + etak1 * power((x1_curr+ Ag01/(2*lg)), 3));
%     end
%     
%     if (x2_curr > -Ag02 / (2 * lg))
%         s2_prime = k2 * etak2*x2_curr*x2_curr*x2_curr;
%     else
%         s2_prime = (k2 * etak2*x2_curr*x2_curr*x2_curr) - h2_upd * (Ag02 / (2 * lg) + etak2 * power((x2_curr + Ag02 / (2 * lg)), 3));
%     end
           
    % Follow (A-11) from S87 to calculate b1 and b2
	% NOTE: You'll find a mismatch between the code implementation and equation
	% that has been provided in the S87 paper in the last term. To understand 
	% the error please validate the units for each terms (which should be: kg*m -in SI)
	b1 = (2 * m1 + r1_upd / SRATE)*x1_curr - m1 * x1_old - s1_prime / (SRATE*SRATE) + f_mass1/(SRATE*SRATE);
	b2 = (2 * m2 + r2_upd / SRATE)*x2_curr - m2 * x2_old - s2_prime / (SRATE*SRATE) + f_mass2/(SRATE*SRATE);
    
    % Compute the determinant 
    det = a11*a22 - a21*a12;
    if det == 0
        det = 1;
    end
    
    x1_next = (a22*b1 - a12 * b2) / det;
	x2_next = (a11*b2 - a21 * b1) / det;

	% Update x1_curr and x1_old
    vf_flowParam.x1_old = x1_curr;
	vf_flowParam.x1_curr = x1_next;
    vf_flowParam.x1_next = x1_next;
    x1_curr = x1_next;
    
	vf_flowParam.x2_old = x2_curr;
	vf_flowParam.x2_curr = x2_next;
    vf_flowParam.x2_next = x2_next;
    x2_curr = x2_next;
    
	% Calculate new area from x1 nad x2
	Ag1 = Ag01 + 2 * lg * x1_curr;
	Ag2 = Ag02 + 2 * lg * x2_curr;
    
    %*********************GLOTTAL FLOW VELOCITY Ug****************************%
    
    if (Ag1 < 0 || Ag2 < 0)
       ug_next = 0;
    else
       % Using Equation 3-4 from S87(p957)
       Rtot = (rho / 2)*abs(ug_curr)*((0.37 / (Ag1*Ag1)) + (1 - 2 * (Ag2 / A1)*(1 - Ag2 / A1)) / (Ag2*Ag2))  + (12 * mu* lg* lg*(d1 / (Ag1*Ag1*Ag1) + d2 / (Ag2*Ag2*Ag2)));
%        Rtot = (rho / 2)*abs(ug_curr)*((0.37 / (Ag2*Ag2)) + (1 - 2 * (Ag2 / A1)*(1 - Ag2 / A1)) / (Ag2*Ag2))  + (12 * mu* lg* lg*(d1 / (Ag1*Ag1*Ag1) + d2 / (Ag2*Ag2*Ag2)));
	   Ltot = rho * (d1 / Ag1 + d2 / Ag2);
        
       % Equation 2 : Png Noise pressure source has not been considered here
	   ug_next = ((ps - p1) / SRATE + Ltot * ug_curr) / (Rtot / SRATE + Ltot);
    end
    
	vf_flowParam.ug_old = ug_curr;
	vf_flowParam.ug_curr = ug_next;
    vf_flowParam.ug_next = ug_next;
end