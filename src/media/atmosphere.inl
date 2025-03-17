#include "../volume.h"
#include <iostream>

Spectrum get_majorant_op::operator()(const AtmosphereMedium &m) {

    Real beta_r = 0.00519673173; //for 680nm 
    Real beta_g = 0.01214269792; //for 550nm
    Real beta_b = 0.02964525861; //for 440nm
    return Spectrum{ beta_r, beta_g, beta_b };
}

Spectrum get_sigma_s_op::operator()(const AtmosphereMedium &m) {
    Real beta_r = 0.00519673173; //for 680nm 
    Real beta_g = 0.01214269792; //for 550nm
    Real beta_b = 0.02964525861; //for 440nm

    
    Real earth_radius = 6378;
    Real dist = distance(p, Vector3{ 0, 0, 0 });

    
    if (dist > 6438 || dist < 6378) {
        return make_zero_spectrum();
    }

    Real height = (dist - earth_radius);

   
    Real scale_height = 8.5;
    Real roh = exp(-height / scale_height);

    // correction factor for the anisotropic properties of air molecules
    Real anisotropic = 1.06081;

    Real sigma_r = beta_r * roh;
    Real sigma_g = beta_g * roh;
    Real sigma_b = beta_b * roh;

    return Spectrum{ sigma_r, sigma_g, sigma_b };
}

Spectrum get_sigma_a_op::operator()(const AtmosphereMedium &m) {

    return make_zero_spectrum();
    
}