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
    int media_id = scene.camera.medium_id;
    Spectrum sigma_a = get_sigma_a(scene.media[media_id], vertex_->position);
    Spectrum sigma_s = get_sigma_s(scene.media[media_id], vertex_->position);
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
            Spectrum phase = eval(get_phase_function(scene.media[media_id]), -ray.dir, dir_light);
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
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    int media_id = scene.camera.medium_id;

    Spectrum current_throughput = make_const_spectrum(1.0);
    Spectrum radiance = make_zero_spectrum();
    Real bounce = 0;

    while(true) {
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);

        Real transmittance = 1;
        Real trans_pdf = 1;

        if(media_id != -1) {
            Spectrum sigma_a = get_sigma_a(scene.media[media_id], vertex_->position);
            Spectrum sigma_s = get_sigma_s(scene.media[media_id], vertex_->position);
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
                media_id = update_medium(*vertex_, ray, media_id);
                // get in to the next media
                ray.org = vertex_->position + ray.dir * get_intersection_epsilon(scene);

                bounce++;
                continue;
            }
        }

        if(scatter) {
            auto medium = scene.media[media_id];
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
    return make_zero_spectrum();
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
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
