% The vocal fold modelling can be divided into two halves:
% 1. Vocal Fold Structural Model
% 2. Vocal Fold Flow Model
% This function initializes and sets the structural and flow parameters for
% a two-mass vocal fold model

function [airParam, vf_structuralParam, vf_flowParam, vf_matParam] = vf_SetVocalFoldParams()

    % DEFINE UNITS AND CONSTANTS
    METER = 1;      % Unit of length
    KILOGRAM = 1;   % Unit of mass
    NEWTON = 1;     % Unit of force
    SECOND = 1;    % Unit of time

    CM     = 1e-2*METER; 
    GRAM   = 1e-3*KILOGRAM; 
    DYN    = 1e-5*NEWTON; 

    airParam.mu = 1.86e-4 * (DYN *SECOND / (CM*CM)); % Air viscocity coefficient
    airParam.rho = 1.14e-3 * (GRAM / (CM*CM*CM)); % Air density
    
    % GLOTTAL MODEL PARAMETERS
    % cord-tension parameter: Determines the mechanical constant of the 
    % oscillator. Read and understand how does it impact the self-oscillation 
    % of the vocal fold.
    vf_structuralParam.q_fact = 1; % Pitch factor(As per S87) =>Yet to verify
    vf_structuralParam.gs     = 1; %Dimensionless damping factor
    
    % Vocal fold mass
    vf_structuralParam.m1 = (0.125*GRAM) / vf_structuralParam.q_fact;
    vf_structuralParam.m2 = (0.025*GRAM) / vf_structuralParam.q_fact;

    % Thickness of each mass
    vf_structuralParam.d1 = (0.25*CM) / vf_structuralParam.q_fact;
    vf_structuralParam.d2 = (0.05*CM) / vf_structuralParam.q_fact;

    % Linear spring stiffness
    vf_structuralParam.k1 = 80000*(DYN/CM)*vf_structuralParam.q_fact;
    vf_structuralParam.k2 = 8000*(DYN/CM)*vf_structuralParam.q_fact;
    vf_structuralParam.kc = 25000 * (DYN / CM)*(vf_structuralParam.q_fact*vf_structuralParam.q_fact); % Coupled spring coefficient
    
    % Non-linear spring stiffness
    vf_structuralParam.etak1 = 100 / (CM*CM);
    vf_structuralParam.etak2 = 100 / (CM*CM);

    % Linear stiffness during vocal fold collision
    vf_structuralParam.h1 = 3 * vf_structuralParam.k1; 
    vf_structuralParam.h2 = 3 * vf_structuralParam.k2; 

    % Non-linear stiffness during vocal fold collision
    vf_structuralParam.etah1 = 500 / (CM*CM);
    vf_structuralParam.etah2 = 500 / (CM*CM);

    % Cross sectional area of glottal slit at rest
    vf_structuralParam.Ag0 = 0.05*CM*CM; %Glottal rest area
    vf_structuralParam.Ag01 = vf_structuralParam.Ag0;
    vf_structuralParam.Ag02 = vf_structuralParam.Ag0;
    
    % Vocal cord effective length
    vf_structuralParam.lg = 1.4*CM;

    % Damping ratio
    vf_structuralParam.zeta1_open = 0.2;
    vf_structuralParam.zeta2_open = 0.6;
    vf_structuralParam.zeta1_close = 1.1;
    vf_structuralParam.zeta2_close = 1.6;
    
    % Calculate viscous resitance
    vf_structuralParam.r1_open = (2 * vf_structuralParam.zeta1_open * ...
        sqrt(vf_structuralParam.k1*vf_structuralParam.m1)) / (vf_structuralParam.gs^2);
    vf_structuralParam.r2_open = (2 * vf_structuralParam.zeta2_open * ...
        sqrt(vf_structuralParam.k2*vf_structuralParam.m2)) / (vf_structuralParam.gs^2);
    vf_structuralParam.r1_close = (2 * vf_structuralParam.zeta1_close * ...
        sqrt(vf_structuralParam.k1*vf_structuralParam.m1)) / (vf_structuralParam.gs^2);
    vf_structuralParam.r2_close = (2 * vf_structuralParam.zeta2_close * ...
        sqrt(vf_structuralParam.k2*vf_structuralParam.m2)) / (vf_structuralParam.gs^2);
    
    % Input area to the vocal tract
    vf_structuralParam.A1 = 3*CM*CM;
    
    % Glottal area changes
    vf_structuralParam.Ag1 = 0;
    vf_structuralParam.Ag2 = 0;
    
    % Sub-glottal(Lungs) pressure
    vf_flowParam.ps = 686.46; % In Pascal(7cm H2O) ==== 1cmH2O => 98.0665 => 784.532

    % Pressure downstream of glottal expansion
    vf_flowParam.p1 = 0;
    
    % DECLARE TIME DEPENDENT VARIABLES
    vf_flowParam.f_mass1 = 0; % Force on mass m1 of vocal cord
    vf_flowParam.f_mass2 = 0; % Force on mass m2 of vocal cord
    vf_flowParam.p_mass1 = 0; % Pressure on mass m1 of vocal cord
    vf_flowParam.p_mass2 = 0; % Pressure on mass m2 of vocal cord

    vf_flowParam.Rv1 = 0; % Glottal resitive parameters
    vf_flowParam.Rv2 = 0; % Glottal resitive parameters
    vf_flowParam.Lg1 = 0; % Glottal inductive parameters
    vf_flowParam.Lg2 = 0; % Glottal inductive parameters
    
    % Flow velocity change
    vf_flowParam.ug_old = 0;
    vf_flowParam.ug_curr = 0;
    vf_flowParam.ug_next = 0;

    % Latteral displacement of vocal fold masses
    vf_flowParam.x1_old = 0;
    vf_flowParam.x1_curr = 0;
    vf_flowParam.x1_next = 0;

    vf_flowParam.x2_old = 0;
    vf_flowParam.x2_curr = 0;
    vf_flowParam.x2_next = 0;
    
    % Matrix Parameters from S87
    vf_matParam.a11 = 0;
    vf_matParam.a12 = 0;
    vf_matParam.a21 = 0;
    vf_matParam.a22 = 0;
    vf_matParam.b1 = 0;
    vf_matParam.b2 = 0;
    
    vf_matParam.s1_prime = 0;
    vf_matParam.s2_prime = 0;
end