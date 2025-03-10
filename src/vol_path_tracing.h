#pragma once

int update_medium(const PathVertex &vertex, const Ray& ray, int current_medium_id) {
    if(vertex.interior_medium_id != vertex.exterior_medium_id) {
        if(dot(ray.dir, vertex.geometric_normal) > 0) {
            return vertex.exterior_medium_id;
        } else {
            return vertex.interior_medium_id;
        }
    }
    return current_medium_id;
}

Spectrum next_event_estimation(Vector3 p, int meidum_id, const Ray& ray_in, const Scene& scene,  pcg32_state &rng) {
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real shape_w = next_pcg32_real<Real>(rng);
    Real light_w = next_pcg32_real<Real>(rng);
    Real light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];

    PointAndNormal p_prime = sample_point_on_light(light, p, light_uv, shape_w, scene);

    Real transmittance_light = 1.0;
    int shadow_medium_id = meidum_id;
    int shadow_bounces = 0;
    Real p_trans_dir = 1; // for mis

    Vector3 p_i = p;

    while (true) {
       Vector3 shadow_dir = normalize(p_prime.position - p);
       Ray shadow_ray{p, 
                    shadow_dir, 
                    get_shadow_epsilon(scene),
                    (1 - get_shadow_epsilon(scene)) * distance(p_prime.position, p)};
        auto isect = intersect(scene, shadow_ray);
        Real next_t = distance(p, p_prime.position);
        if(isect) {
            next_t = distance(p, isect->position);
        } 
        if( shadow_medium_id != -1) {
            Spectrum sigma_a = get_sigma_a(scene.media[shadow_medium_id], p);
            Spectrum sigma_s = get_sigma_s(scene.media[shadow_medium_id], p);
            Spectrum sigma_t = sigma_a + sigma_s;
            transmittance_light *= exp(-sigma_t[0] * next_t);
            p_trans_dir *= exp(-sigma_t[0] * next_t);
        }
        if(!isect) {
            break;
        } 
        // blocked by a surface
        if(isect->material_id > 0) {
            return make_zero_spectrum();
        }

        shadow_bounces++;
        if(scene.options.max_depth != -1 && shadow_bounces >= scene.options.max_depth) {
            return make_zero_spectrum();
        }
        shadow_medium_id = update_medium(*isect, shadow_ray, shadow_medium_id);
        p = p + shadow_dir * next_t;
       
    }

    if (transmittance_light > 0) {
        Real dist_sqr = distance_squared(p_i, p_prime.position);
        Vector3 dir = normalize(p_i - p_prime.position);
        Real cos_theta = abs(dot(dir, p_prime.normal));
        Real G= cos_theta/ dist_sqr;

        Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, p_prime, p_i, scene);
        PhaseFunction phase_function = get_phase_function(scene.media[meidum_id]);
        Spectrum rho = eval(phase_function, -ray_in.dir, -dir);
        Spectrum L = emission(light, dir, 0, p_prime, scene);
        Real pdf_phase = pdf_sample_phase(phase_function, -ray_in.dir, -dir) * G * p_trans_dir;

        Real weight = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
        Spectrum contrib = transmittance_light * G * rho * L / pdf_nee;

        return weight * contrib;

    }

    return make_zero_spectrum();
}

Spectrum next_event_estimation_bsdf(Vector3 p, int meidum_id, const Ray& ray_in, const Scene& scene,  pcg32_state &rng, const PathVertex& vertex) {
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real shape_w = next_pcg32_real<Real>(rng);
    Real light_w = next_pcg32_real<Real>(rng);
    Real light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];

    PointAndNormal p_prime = sample_point_on_light(light, p, light_uv, shape_w, scene);

    Real transmittance_light = 1.0;
    int shadow_medium_id = meidum_id;
    int shadow_bounces = 0;
    Real p_trans_dir = 1; // for mis

    Vector3 p_i = p;

    while (true) {
       Vector3 shadow_dir = normalize(p_prime.position - p);
       Ray shadow_ray{p, 
                    shadow_dir, 
                    get_shadow_epsilon(scene),
                    (1 - get_shadow_epsilon(scene)) * distance(p_prime.position, p)};
        auto isect = intersect(scene, shadow_ray);
        Real next_t = distance(p, p_prime.position);
        if(isect) {
            next_t = distance(p, isect->position);
        } 
        if( shadow_medium_id != -1) {
            Spectrum sigma_a = get_sigma_a(scene.media[shadow_medium_id], p);
            Spectrum sigma_s = get_sigma_s(scene.media[shadow_medium_id], p);
            Spectrum sigma_t = sigma_a + sigma_s;
            transmittance_light *= exp(-sigma_t[0] * next_t);
            p_trans_dir *= exp(-sigma_t[0] * next_t);
        }
        if(!isect) {
            break;
        } 
        // blocked by a surface
        if(isect->material_id > 0) {
            return make_zero_spectrum();
        }

        shadow_bounces++;
        if(scene.options.max_depth != -1 && shadow_bounces >= scene.options.max_depth) {
            return make_zero_spectrum();
        }
        shadow_medium_id = update_medium(*isect, shadow_ray, shadow_medium_id);
        p = p + shadow_dir * next_t;
       
    }

    if (transmittance_light > 0) {
        Real dist_sqr = distance_squared(p_i, p_prime.position);
        Vector3 dir = normalize(p_i - p_prime.position);
        Real cos_theta = abs(dot(dir, p_prime.normal));
        Real G= cos_theta/ dist_sqr;

        const Material& mat = scene.materials[vertex.material_id];

        Spectrum L = emission(light, dir, 0, p_prime, scene);
        Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, p_prime, p_i, scene);

        Spectrum bsdf = eval(mat, -ray_in.dir, -dir, vertex, scene.texture_pool);
        Real pdf_bsdf = pdf_sample_bsdf(mat, -ray_in.dir, -dir, vertex, scene.texture_pool) * G * p_trans_dir; 
        Real weight = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_bsdf * pdf_bsdf);
        Spectrum contrib = transmittance_light * G * bsdf * L / pdf_nee;

        return weight * contrib;

    }

     return make_zero_spectrum();
}


// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);

    if (vertex_ && is_light(scene.shapes[vertex_->shape_id])) {
        int medium_id = vertex_->exterior_medium_id;
        Spectrum sigma_a = get_sigma_a(scene.media[medium_id], vertex_->position);
        Real t = distance(ray.org, vertex_->position);
        Spectrum transmittance = exp(-sigma_a * t);
        return transmittance * emission(*vertex_, -ray.dir, scene);
    }


    return make_zero_spectrum();
}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!

    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    int medium_id = scene.camera.medium_id;
    Spectrum sigma_a = get_sigma_a(scene.media[medium_id], vertex_->position);
    Spectrum sigma_s = get_sigma_s(scene.media[medium_id], vertex_->position);
    Spectrum sigma_t = sigma_a + sigma_s;

    Real u = next_pcg32_real<Real>(rng);
    Real t = -log(1-u)/ sigma_t[0];

    if(!vertex_ || t < distance(ray.org, vertex_->position)){
        Real trans_pdf = exp(-sigma_t[0] * t) * sigma_t[0];
        Real transmittance = exp(-sigma_t[0] * t);
        Vector3 p = ray.org + ray.dir * t;


        
        Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real light_w = next_pcg32_real<Real>(rng);
        Real shape_w = next_pcg32_real<Real>(rng);
        int light_id = sample_light(scene, light_w);
        const Light &light = scene.lights[light_id];
        PointAndNormal point_on_light =
            sample_point_on_light(light, p, light_uv, shape_w, scene);

        Vector3 dir_light = normalize(point_on_light.position - p);
        Real L_s1_pdf = light_pmf(scene, light_id) * pdf_point_on_light(light, point_on_light, p, scene);
       
        Ray shadow_ray{p, dir_light, 
            get_shadow_epsilon(scene),
            (1 - get_shadow_epsilon(scene)) *
                distance(point_on_light.position, p)};
        
        Spectrum L_s1_estimate = make_zero_spectrum();

        if (!occluded(scene, shadow_ray)) {
            Spectrum Le = emission(light, -dir_light, Real(0), point_on_light, scene);
            Spectrum phase = eval(get_phase_function(scene.media[medium_id]), -ray.dir, dir_light);
            Real cos_theta  = abs(dot(dir_light, point_on_light.normal));
            Real dist = distance(p, point_on_light.position);
           L_s1_estimate = phase * Le  * exp(-sigma_t[0] * dist) * cos_theta / (dist * dist);
        }

        return  (transmittance /trans_pdf) * sigma_s * L_s1_estimate / L_s1_pdf;
    } else {
        Real dist = distance(ray.org, vertex_->position);

        Real trans_pdf = exp(-sigma_t[0] * dist);
        Real transmittance = exp(-sigma_t[0] * dist);

        Spectrum Le = make_zero_spectrum();
        if(vertex_ && is_light(scene.shapes[vertex_->shape_id])){
            Le =  emission(*vertex_, -ray.dir, scene);
        }
        return transmittance / trans_pdf * Le;

    }
    return make_zero_spectrum();
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    int medium_id = scene.camera.medium_id;

    Spectrum current_throughput = make_const_spectrum(1.0);
    Spectrum radiance = make_zero_spectrum();
    Real bounce = 0;

    while(true) {
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray);

        Real transmittance = 1;
        Real trans_pdf = 1;

        if(medium_id != -1) {
            Spectrum sigma_a = get_sigma_a(scene.media[medium_id], vertex_->position);
            Spectrum sigma_s = get_sigma_s(scene.media[medium_id], vertex_->position);
            Spectrum sigma_t = sigma_a + sigma_s;
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1-u)/ sigma_t[0];

            Real t_hit = distance(ray.org, vertex_->position);

            if(!vertex_ || t < t_hit) {
                scatter = true;
                transmittance = exp(-sigma_t[0] * t);
                trans_pdf = exp(-sigma_t[0] * t) * sigma_t[0];
                ray.org = ray.org + ray.dir * t;
            } else {
                transmittance = exp(-sigma_t[0] * t_hit);
                trans_pdf = exp(-sigma_t[0] * t_hit);
                ray.org = ray.org + ray.dir * t_hit;
            }
        }

        current_throughput *= transmittance/trans_pdf;

        if(!scatter && is_light(scene.shapes[vertex_->shape_id])) {
            radiance += current_throughput * emission(*vertex_, -ray.dir, scene);
        }

        if(bounce == scene.options.max_depth - 1) {
            break;
        }

        if(!scatter && vertex_) {
            if(vertex_->material_id == -1) {
                medium_id = update_medium(*vertex_, ray, medium_id);
                // get in to the next media
                ray.org = vertex_->position + ray.dir * get_intersection_epsilon(scene);
                bounce++;
                continue;
            }
        }

        if(scatter) {
            auto medium = scene.media[medium_id];
            auto phase_function = get_phase_function(medium);
            Vector2 randomUV = {next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            std::optional<Vector3> next_dir = sample_phase_function(
                phase_function, -ray.dir, randomUV
            );
            if(next_dir) {
                Real sigma_s = get_sigma_s(medium, vertex_->position)[0];
                current_throughput *= (eval(phase_function, -ray.dir, *next_dir) 
                    / pdf_sample_phase(phase_function, -ray.dir, *next_dir)) * sigma_s;
                ray.dir = *next_dir;
            } 
        } else {
            break;
        }

        if(bounce >= scene.options.rr_depth) {
            Real rr_prob = min(luminance(current_throughput), 0.95);
            if(next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_throughput /= rr_prob;
            }
        }
        bounce++;
    }
    return radiance;
}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!


    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);

    Ray ray = sample_primary(scene.camera, screen_pos);
    
    int medium_id = scene.camera.medium_id;

    Spectrum current_throughput = make_const_spectrum(1.0);
    Spectrum radiance = make_zero_spectrum();
    Real bounces = 0;

    Real dir_pdf = 0;
    Real mis_trans_pdf = 1; 
    Vector3 nee_cache;
    bool never_scatter = true;
 
    while (true) {
        bool scatter = false; 
        std::optional<PathVertex> vertex_ = intersect(scene, ray); 
        Real transmittance = 1;
        Real trans_pdf = 1;

        if(medium_id != -1) {
            Spectrum sigma_a = get_sigma_a(scene.media[medium_id], vertex_->position);
            Spectrum sigma_s = get_sigma_s(scene.media[medium_id], vertex_->position);
            Spectrum sigma_t = sigma_a + sigma_s;
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1-u)/ sigma_t[0];

            Real t_hit = distance(ray.org, vertex_->position);

            if(!vertex_ || t < t_hit) {
                scatter = true;
                transmittance = exp(-sigma_t[0] * t);
                trans_pdf = exp(-sigma_t[0] * t) * sigma_t[0];
                ray.org = ray.org + ray.dir * t;
            } else {
                transmittance = exp(-sigma_t[0] * t_hit);
                trans_pdf = exp(-sigma_t[0] * t_hit);
                ray.org = ray.org + ray.dir * t_hit;
            }
        }    
         current_throughput *= (transmittance / trans_pdf);
 
         if(!scatter && vertex_ && is_light(scene.shapes[vertex_->shape_id])) {
            if(never_scatter)
                radiance += current_throughput * emission(*vertex_, -ray.dir, scene);
            else {
                int light_id = get_area_light_id(scene.shapes[vertex_->shape_id]);
                PointAndNormal light_point = {vertex_->position, vertex_->geometric_normal};
                Real pdf_nee = pdf_point_on_light(scene.lights[light_id], light_point, nee_cache, scene);
                Light light = scene.lights[light_id];
                Vector3 ray_prime = normalize(vertex_->position - nee_cache);
                Real G = abs(dot(ray_prime, light_point.normal)) / distance_squared(vertex_->position, nee_cache);
                Real dir_pdf_ = dir_pdf * mis_trans_pdf * G;
                Real weight = (dir_pdf_ * dir_pdf_) / (pdf_nee * pdf_nee + dir_pdf_ * dir_pdf_);
                radiance += current_throughput * emission(*vertex_, -ray.dir, scene) * weight;
            }
        }

 
        if(bounces == scene.options.max_depth - 1) {
            break;
        }

        if(!scatter && vertex_) {
            if(vertex_->material_id == -1) {
                medium_id = update_medium(*vertex_, ray, medium_id);
                // get in to the next media
                ray.org = vertex_->position + ray.dir * get_intersection_epsilon(scene);
                mis_trans_pdf *= trans_pdf;
                bounces++;
                continue;
            }
        }
         
        if(scatter) {
            never_scatter = false;
            auto medium = scene.media[medium_id];
            auto phase_function = get_phase_function(medium);
            Vector2 randomUV = {next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            std::optional<Vector3> next_dir = sample_phase_function(
                phase_function, -ray.dir, randomUV
            );
            if(next_dir) {
                Real sigma_s = get_sigma_s(medium, vertex_->position)[0];
                Spectrum nee = next_event_estimation(ray.org, medium_id, ray, scene, rng);
                radiance += current_throughput * nee * sigma_s;
                nee_cache = ray.org;
                current_throughput *= (eval(phase_function, -ray.dir, *next_dir) 
                    / pdf_sample_phase(phase_function, -ray.dir, *next_dir)) * sigma_s;
                ray.dir = *next_dir;
                mis_trans_pdf = 1;
                
            } 
        } else {
            break;
        }

        if(bounces >= scene.options.rr_depth) {
            Real rr_prob = min(luminance(current_throughput), 0.95);
            if(next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_throughput /= rr_prob;
            }
        }
        bounces++;
     }
     return radiance;
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);

    Ray ray = sample_primary(scene.camera, screen_pos);
    
    int medium_id = scene.camera.medium_id;

    Spectrum current_throughput = make_const_spectrum(1.0);
    Spectrum radiance = make_zero_spectrum();
    Real bounces = 0;

    Real dir_pdf = 0;
    Real mis_trans_pdf = 1; 
    Vector3 nee_cache;
    bool never_scatter = true;
 
    while (true) {
        bool scatter = false; 
        std::optional<PathVertex> vertex_ = intersect(scene, ray); 
        Real transmittance = 1;
        Real trans_pdf = 1;

        
        if (medium_id == -1) {
            if (vertex_) {
                PathVertex vertex = vertex_.value();
                ray.org = ray.org +(distance(ray.org, vertex.position)  + get_intersection_epsilon(scene) )* ray.dir;
            } else break;
        }

        if(medium_id != -1) {
            Spectrum sigma_a = get_sigma_a(scene.media[medium_id], vertex_->position);
            Spectrum sigma_s = get_sigma_s(scene.media[medium_id], vertex_->position);
            Spectrum sigma_t = sigma_a + sigma_s;
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1-u)/ sigma_t[0];

            Real t_hit = distance(ray.org, vertex_->position);

            if(!vertex_ || t < t_hit) {
                scatter = true;
                transmittance = exp(-sigma_t[0] * t);
                trans_pdf = exp(-sigma_t[0] * t) * sigma_t[0];
                ray.org = ray.org + ray.dir * t;
            } else {
                transmittance = exp(-sigma_t[0] * t_hit);
                trans_pdf = exp(-sigma_t[0] * t_hit);
                ray.org = ray.org + ray.dir * t_hit;
            }
        }    
         current_throughput *= (transmittance / trans_pdf);
 
         if(!scatter && vertex_ && is_light(scene.shapes[vertex_->shape_id])) {
            if(never_scatter)
                radiance += current_throughput * emission(*vertex_, -ray.dir, scene);
            else {
                int light_id = get_area_light_id(scene.shapes[vertex_->shape_id]);
                PointAndNormal light_point = {vertex_->position, vertex_->geometric_normal};
                Real pdf_nee = pdf_point_on_light(scene.lights[light_id], light_point, nee_cache, scene);
                Light light = scene.lights[light_id];
                Vector3 ray_prime = normalize(vertex_->position - nee_cache);
                Real G = abs(dot(ray_prime, light_point.normal)) / distance_squared(vertex_->position, nee_cache);
                Real dir_pdf_ = dir_pdf * mis_trans_pdf * G;
                Real weight = (dir_pdf_ * dir_pdf_) / (pdf_nee * pdf_nee + dir_pdf_ * dir_pdf_);
                radiance += current_throughput * emission(*vertex_, -ray.dir, scene) * weight;
            }
        }

 
        if(bounces == scene.options.max_depth - 1) {
            break;
        }

        if(!scatter && vertex_) {
            if(vertex_->material_id == -1) {
                medium_id = update_medium(*vertex_, ray, medium_id);
                // get in to the next media
                ray.org = vertex_->position + ray.dir * get_intersection_epsilon(scene);
                mis_trans_pdf *= trans_pdf;
                bounces++;
                continue;
            }
        }
         
        if(scatter) {
            never_scatter = false;
            auto medium = scene.media[medium_id];
            auto phase_function = get_phase_function(medium);
            Vector2 randomUV = {next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            std::optional<Vector3> next_dir = sample_phase_function(
                phase_function, -ray.dir, randomUV
            );
            if(next_dir) {
                Real sigma_s = get_sigma_s(medium, vertex_->position)[0];
                Spectrum nee = next_event_estimation(ray.org, medium_id, ray, scene, rng);
                radiance += current_throughput * nee * sigma_s;
                current_throughput *= (eval(phase_function, -ray.dir, *next_dir) 
                    / pdf_sample_phase(phase_function, -ray.dir, *next_dir)) * sigma_s;
                
                nee_cache = ray.org;
                ray.dir = *next_dir;
                mis_trans_pdf = 1;
                
            } 

        } else if(vertex_) {
            // hit a surface
            never_scatter = false;
            Spectrum nee = next_event_estimation_bsdf(ray.org, medium_id, ray, scene, rng, *vertex_);
            radiance += current_throughput * nee;
        

            const Material &mat = scene.materials[vertex_->material_id];
            Vector3 dir_view = -ray.dir;
            Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
            std::optional<BSDFSampleRecord> bsdf_sample_ =
                sample_bsdf(mat,
                            dir_view,
                            *vertex_,
                            scene.texture_pool,
                            bsdf_rnd_param_uv,
                            bsdf_rnd_param_w);

             if (!bsdf_sample_) {
                // BSDF sampling failed. Abort the loop.
                break;
            }

            const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;

            Vector3 dir_bsdf = bsdf_sample.dir_out;
            ray.dir = dir_bsdf;
            medium_id = update_medium(*vertex_, ray, medium_id);

            Spectrum f = eval(mat, dir_view, dir_bsdf, *vertex_, scene.texture_pool);
            Real bsdf_pdf = pdf_sample_bsdf(mat, dir_view, dir_bsdf, *vertex_, scene.texture_pool);
            current_throughput *= f / bsdf_pdf;


            mis_trans_pdf = 1;
            nee_cache = ray.org;
        }

        if(bounces >= scene.options.rr_depth) {
            Real rr_prob = min(luminance(current_throughput), 0.95);
            if(next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_throughput /= rr_prob;
            }
        }
        bounces++;
     }
     return radiance;
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}
