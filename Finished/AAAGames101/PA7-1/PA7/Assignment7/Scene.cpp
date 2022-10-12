//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    Intersection ip = intersect(ray);
    if (!ip.happened)
        return Vector3f(0.0f, 0.0f, 0.0f);
    if (ip.m->hasEmission())
        return ip.m->getEmission();

    Intersection ix;
    float pdf_light = 1.0f;
    sampleLight(ix, pdf_light);
    Vector3f& p = ip.coords;
    Vector3f& x = ix.coords;

    Ray p2x(p, normalize(x - p));
    Intersection ip2x = intersect(p2x);

    Vector3f L_direct(0.0f, 0.0f, 0.0f);
    if ( (p-x).norm() - ip2x.distance < 0.1f) //bounce ray hit light, we assume ix equals ip2x
    {
        L_direct = ix.emit * ip.m->eval(ray.direction, p2x.direction, ip.normal) * 
            dotProduct(ip.normal, p2x.direction) * dotProduct(-p2x.direction, ix.normal) / (p-x).len2() / pdf_light;
    }

    Vector3f L_indirect = {};
    float num = get_random_float();
    if (num > RussianRoulette)
        return L_direct;
    Vector3f other_dir = ip.m->sample(ray.direction, ip.normal).normalized();
    Ray other_ray(p, other_dir);
    Intersection io = intersect(other_ray);
    if (io.happened && !io.m->hasEmission())
    {
        L_indirect = ip.m->eval(ray.direction, other_dir, ip.normal) * castRay(other_ray, depth + 1) *
            dotProduct(other_dir, ip.normal) / ip.m->pdf(ray.direction, other_dir, ip.normal) / RussianRoulette;
    }

    return L_direct + L_indirect;
}