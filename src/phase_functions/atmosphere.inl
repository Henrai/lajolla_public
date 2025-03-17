Spectrum eval_op::operator()(const AtomspherePhase &) const {
    Real mu = -dot(dir_in, dir_out);
    Spectrum rayleigh = make_const_spectrum(c_INVFOURPI * 0.75 * (Real(1) + mu * mu));
    return rayleigh;
}

std::optional<Vector3> sample_phase_function_op::operator()(const AtomspherePhase &) const {
   

    Real term = 2 * rnd_param.x - 1;
    Real term1 = 2 * term;
    Real term2 = sqrt(4 * term * term + 1);
    Real u = -pow(term1 + term2, Real(1) / Real(3));
    Real cos_elevation = u - 1 / u;
    Real sin_elevation = sqrt(max(1 - cos_elevation * cos_elevation, Real(0)));
    Real azimuth = 2 * c_PI * rnd_param.y;
    Frame frame(dir_in);
    return to_world(frame,
                    Vector3{ sin_elevation * cos(azimuth),
                    sin_elevation * sin(azimuth),
                    cos_elevation });
}

Real pdf_sample_phase_op::operator()(const AtomspherePhase &) const {
    Real mu = -dot(dir_in, dir_out);
    Real rayleigh_pdf = c_INVFOURPI * 0.75 * (Real(1) + mu * mu);
    return rayleigh_pdf;
}
