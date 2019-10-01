class SamplerProgram extends WebGL2ProgramContainer {

    constructor(gl) {
        const glUtils = WebGL2Utils.from(gl);

        const program = glUtils.createProgramFromSources(
            SamplerProgram.vertShaderSource.trim(), SamplerProgram.fragShaderSource.trim());

        super(gl, program);

        const cameraBlockBinding = 0;
        const cameraBuffer = new glUtils.UniformBuffer(cameraBlockBinding, gl.DYNAMIC_DRAW);
        gl.uniformBlockBinding(program, gl.getUniformBlockIndex(program, "CameraGeometry"), cameraBlockBinding);

        this.setCamera = (location, focalLength, cameraX, cameraY, cameraZ) => {
            // Only std140 layout supported in uniform blocks
            cameraBuffer.update([
                ...location, focalLength,
                ...cameraX, 0,
                ...cameraY, 0,
                ...cameraZ, 0
            ]);
        };
    }

}

SamplerProgram.vertShaderSource = `
    #version 300 es

    const vec2 vertices[4] = vec2[] (
        vec2(-1.0, 1.0),
        vec2(-1.0, -1.0),
        vec2(1.0, 1.0),
        vec2(1.0, -1.0)
    );

    void main() {
        gl_Position = vec4(vertices[gl_VertexID], 0.0, 1.0);
    }
`;

SamplerProgram.fragShaderSource = `
    #version 300 es

    precision highp float;
    precision highp sampler2DArray;

    const float PI = radians(180.0);
    const float TWOPI = 2.0*PI;
    const float FOURPI = 4.0*PI;
    const float MAX_VALUE = intBitsToFloat(0x7f7fffff);
    const float INF = intBitsToFloat(0x7f800000);

    struct Material {
        uint type;
        vec4 data;
    };

    const uint M_EMISSIVE = 0u;
    const uint M_DIFFUSE = 1u;
    const uint M_DIELECTRIC = 3u;

    #define EMISSIVE(r, g, b, i) Material(M_EMISSIVE, i*vec4(r, g, b, 0.0))
    #define DIFFUSE(r, g, b) Material(M_DIFFUSE, vec4(r, g, b, 0.0))
    #define GLASS(ior) Material(M_DIELECTRIC, vec4(1.0, 1.0, 1.0, ior))

    struct Ray {
        vec3 origin;
        vec3 direction;
        float t;
    };

    struct Interaction {
        vec3 outgoing;
        vec3 location;
        vec3 normal;
        Material material;
    };

    struct Primitive {
        uint type;
        vec3[4] data;
        Material material;
    };


    const uint P_HALFSPACE = 1u;
    const uint P_SPHERE = 2u;

    #define HALFSPACE(location, normal, material) Primitive(P_HALFSPACE, vec3[](location, normal, vec3(0.0), vec3(0.0)), material)
    #define SPHERE(center, radius, material) Primitive(P_SPHERE, vec3[](center, vec3(radius, 0.0, 0.0), vec3(0.0), vec3(0.0)), material)

    uniform uint bounces;

    uniform CameraGeometry {
        vec3 cameraLocation;
        float focalLength;
        vec3 cameraX;
        vec3 cameraY;
        vec3 cameraZ;
    };

    uniform uint nSample;
    uniform sampler2DArray pixelData;

    out vec4[2] outputData;

    bool IMPORTANCE_SAMPLING;



    highp uint intHash(highp uint x) {
        x ^= x >> 17;
        x *= 0xED5AD4BBu;
        x ^= x >> 11;
        x *= 0xAC4C1B51u;
        x ^= x >> 15;
        x *= 0x31848BABu;
        x ^= x >> 14;
        return x;
    } // Source: https://nullprogram.com/blog/2018/07/31/

    highp float toFloat(highp uint m) {
        m &= 0x007FFFFFu; // Extract mantissa bits
        m |= 0x3F800000u; // Add mantissa bits to 1.0
        return uintBitsToFloat(m) - 1.0;
    } // Source: https://stackoverflow.com/a/17479300

    highp uint rngIndex;
    highp uint rngOffset;
    highp float nextRandom() {
        return toFloat(intHash(rngOffset + rngIndex++));
    }

    vec3 orientedNormal(Interaction i) {
        return dot(i.outgoing, i.normal) > 0.0 ? -i.normal : i.normal;
    }

    vec3 transformLocal(vec3 v, vec3 normal) {
        float s = normal.z < 0.0 ? -1.0 : 1.0; // can't use sign() because of value 0
        float a = -1.0/(s + normal.z);
        float b = a*normal.x*normal.y;

        vec3 b1 = vec3(1.0 + s*a*normal.x*normal.x, s*b, -s*normal.x);
        vec3 b2 = vec3(b, s + a*normal.y*normal.y, -normal.y);

        return v.x*b1 + v.y*b2 + v.z*normal;
    } // Source: Tom Duff et al., "Building an Orthonormal Basis, Revisited"

    vec3 sampleHemisphereCosine(vec3 normal) {
        float theta = TWOPI*nextRandom();
        float u = 2.0*nextRandom() - 1.0;
        return normalize(normal + vec3(sqrt(1.0 - u*u) * vec2(cos(theta), sin(theta)), u));
    } // Source: http://www.amietia.com/lambertnotangent.html

    void offsetLocation(inout vec3 location, vec3 normal) {
        const float threshold = 1.0/32.0;
        const float floatScale = 1.0/65536.0;
        const float intScale = 256.0;

        vec3 o = intScale*normal;
        highp ivec3 intOffset = ivec3(mix(o, -o, lessThan(location, vec3(0.0))));
        vec3 intOffsetLocation = intBitsToFloat(floatBitsToInt(location) + intOffset);
        vec3 floatOffsetLocation = location + floatScale*normal;

        location = mix(intOffsetLocation, floatOffsetLocation, lessThan(abs(location), vec3(threshold)));
    } // Source: Carsten Wächter and Nikolaus Binder, "A Fast and Robust Method for Avoiding Self-Intersection"

    vec2 intersectHalfspace(vec3 location, vec3 normal, Ray r, out vec3[2] locations, out vec3[2] normals) {
        float d = dot(r.direction, normal);
        if (d == 0.0) return vec2(0.0, 0.0);

        float t = dot(location - r.origin, normal)/d;
        vec3 l = r.origin + t*r.direction;

        locations = vec3[](l, l);
        normals = vec3[](normal, normal);
        return d < 0.0 ? vec2(-INF, t) : vec2(t, INF);
    }

    vec2 intersectSphere(vec3 center, float radius, Ray r, out vec3[2] locations, out vec3[2] normals) {
        vec3 f = r.origin - center;

        float r2 = radius*radius;
        float b = -dot(f, r.direction); // b < 0 => center is behind the ray
        vec3 g = f + b*r.direction;

        float D = r2 - dot(g, g);
        if (D <= 0.0) return vec2(0.0);

        float c = dot(f, f) - r2; // c > 0 => ray starts outside
        float sqrtD = sqrt(D);
        float q = b > 0.0 ? b + sqrtD : b - sqrtD;
        vec3 v = sqrtD*r.direction;

        /*normals = vec3[](normalize(g - v), normalize(g + v));
        locations = vec3[](center + radius*normals[0], center + radius*normals[1]);
        vec2 interval = vec2(c/q, q);

        bool swap = interval[0] > interval[1];
        //normals = vec3[](normals[int(swap)], normals[int(!swap)]);
        //locations = vec3[](locations[int(swap)], locations[int(!swap)]);
        return swap ? interval.ts : interval;*/

        //FIXME: For testing purposes.
        vec2 interval = vec2(min(c/q, q), max(c/q, q));
        locations = vec3[](r.origin + interval[0]*r.direction, r.origin + interval[1]*r.direction);
        normals = vec3[](normalize(locations[0] - center), normalize(locations[1] - center));
        return interval;
    } // Source: Eric Haines et al., "Precision Improvements for Ray/Sphere Intersection"

    vec3 sampleSphere(vec3 center, float radius, out vec3 normal, out float area) {
        float u1 = nextRandom();
        float u2 = nextRandom();
        float a = 2.0*sqrt(u1*(1.0 - u1));
        float b = TWOPI*u2;

        normal = vec3(a*cos(b), a*sin(b), 1.0 - 2.0*u1);
        area = FOURPI*radius*radius;

        return center + radius*normal;
    }

    void intersectCylinder(vec3 pa, vec3 pb, float radius, Material m, inout Ray r, inout Interaction i) {
        vec3 ca = pb - pa;
        vec3 oc = r.origin - pa;

        float caca = dot(ca, ca);
        float card = dot(ca, r.direction);
        float caoc = dot(ca, oc);
        float a = caca - card*card;
        float b = caca*dot(oc, r.direction) - caoc*card;
        float c = caca*dot(oc, oc) - caoc*caoc - radius*radius*caca;
        float h = b*b - a*c;

        if (h < 0.0) return;
        h = sqrt(h);

        vec2 T = vec2((-b - h)/a, (-b + h)/a);
        vec2 Y = vec2(caoc) + T*card;

        bvec2 valid = bvec2(
            0.0 < T[0] && 0.0 < Y[0] && Y[0] < caca,
            0.0 < T[1] && 0.0 < Y[1] && Y[1] < caca
        );

        int index = int(!valid[0]);
        float t = valid[index] ? T[index] : 0.0;

        if (0.0 < t && t < r.t) {
            r.t = t;
            i = Interaction(r.direction, r.origin + t*r.direction, (oc + t*r.direction - ca*Y[index]/caca)/radius, m);
        };
    } // Adapted from https://www.iquilezles.org/www/articles/intersectors/intersectors.htm

    void intersectBox(vec3 location, vec3 dimensions, Material m, inout Ray r, inout Interaction i) {
        vec3 a = 1.0/r.direction;
        vec3 b = a*(r.origin - location);
        vec3 c = 0.5*abs(a)*dimensions;

        vec3 t1 = -b - c;
        vec3 t2 = -b + c;
        float tN = max(max(t1.x, t1.y), t1.z);
        float tF = min(min(t2.x, t2.y), t2.z);

        if (tN > tF || tF < 0.0) return;

        float t = tN < 0.0 ? tF : tN;
        vec3 normal = -sign(r.direction)*step(t1.yzx, t1.xyz)*step(t1.zxy, t1.xyz); //FIXME: This is probably wrong when ray starts inside.

        if (t < r.t) {
            r.t = t;
            i = Interaction(r.direction, r.origin + t*r.direction, normal, m);
        }
    } // Adapted from https://www.iquilezles.org/www/articles/intersectors/intersectors.htm

    vec2 intersectPrimitive(Primitive p, Ray r, out vec3[2] locations, out vec3[2] normals) {
        switch (p.type) {
            case P_HALFSPACE:
                return intersectHalfspace(p.data[0], p.data[1], r, locations, normals);
            case P_SPHERE:
                return intersectSphere(p.data[0], p.data[1][0], r, locations, normals);
            default:
                // Invalid primitive type.
                return vec2(INF);
        }
    }

    vec3 samplePrimitive(Primitive p, out vec3 normal, out float area) {
        switch (p.type) {
            case P_HALFSPACE:
                return vec3(0.0); // Not implemented.
            case P_SPHERE:
                return sampleSphere(p.data[0], p.data[1][0], normal, area);
            default:
                // Invalid primitive type.
                return vec3(0.0);
        }
    }

    void intersectObject(Primitive p, inout Ray r, inout Interaction i) {
        vec3[2] locations;
        vec3[2] normals;
        vec2 interval = intersectPrimitive(p, r, locations, normals);

        int index = int(interval[0] <= 0.0);
        float t = interval[index];

        if (0.0 < t && t < r.t) {
            r.t = t;
            i = Interaction(r.direction, locations[index], normals[index], p.material);
        }
    }

    int skipNonPositive(vec4 v) {
        bvec3 b = bvec3(v[0] <= 0.0, v[1] <= 0.0, v[2] <= 0.0);
        return int(b.s) + int(b.s && b.t) + int(b.s && b.t && b.p);
    }

    void intersectUnion(Primitive[2] p, inout Ray r, inout Interaction i) {
        vec3[2] l1;
        vec3[2] l2;
        vec3[2] n1;
        vec3[2] n2;
        vec2[2] T = vec2[](
            intersectPrimitive(p[0], r, l1, n1),
            intersectPrimitive(p[1], r, l2, n2)
        );

        // Sort intersections: min(T[0].s, T[1].s), min(T[0].t, T[1].t), max(T[0].s, T[1].s), max(T[0].t, T[1].t)
        bool lower = T[1].s < T[0].s;
        bool upper = T[1].t < T[0].t;
        ivec4 pMap = ivec4(int(lower), int(upper), int(!lower), int(!upper));
        ivec4 iMap = ivec4(0, 1, 0, 1);
        vec4 tValues = vec4(T[pMap[0]][iMap[0]], T[pMap[1]][iMap[1]], T[pMap[2]][iMap[2]], T[pMap[3]][iMap[3]]);

        bool skipInner = tValues[1] >= tValues[2];
        if (skipInner) {
            tValues[1] = 0.0;
            tValues[2] = 0.0;
        }

        int index = skipNonPositive(tValues);
        int primitive = pMap[index];
        int intersection = iMap[index];

        float t = tValues[index];
        vec3 location = vec3[2](l1[intersection], l2[intersection])[primitive];
        vec3 normal = vec3[2](n1[intersection], n2[intersection])[primitive];

        if (0.0 < t && t < r.t) {
            r.t = t;
            i = Interaction(r.direction, location, normal, p[primitive].material);
        }
    }

    void intersectIntersection(Primitive[2] p, inout Ray r, inout Interaction i) {
        vec3[2] l1;
        vec3[2] l2;
        vec3[2] n1;
        vec3[2] n2;
        vec2[2] T = vec2[](
            intersectPrimitive(p[0], r, l1, n1),
            intersectPrimitive(p[1], r, l2, n2)
        );

        ivec2 pMap = ivec2(int(T[1].s > T[0].s), int(T[1].t < T[0].t));
        ivec2 iMap = ivec2(0, 1);
        vec2 tValues = vec2(T[pMap[0]][iMap[0]], T[pMap[1]][iMap[1]]);

        if (tValues[0] >= tValues[1]) return;

        int index = int(tValues[0] <= 0.0);
        int primitive = pMap[index];
        int intersection = iMap[index];

        float t = tValues[index];
        vec3 location = vec3[2](l1[intersection], l2[intersection])[primitive];
        vec3 normal = vec3[2](n1[intersection], n2[intersection])[primitive];

        if (0.0 < t && t < r.t) {
            r.t = t;
            i = Interaction(r.direction, location, normal, p[primitive].material);
        }
    }

    Primitive[2] createLensCSG(vec3 location, vec3 direction, float r1, float r2, float thickness, float refractiveIndex) {
        float d = r1 + r2 - thickness;
        float a = 0.5*(r1*r1 - r2*r2 + d*d)/d;
        float b = d - a;

        vec3 l1 = location + a*direction;
        vec3 l2 = location - b*direction;

        Primitive s1 = SPHERE(l1, r1, GLASS(refractiveIndex));
        Primitive s2 = SPHERE(l2, r2, GLASS(refractiveIndex));
        return Primitive[](s1, s2);
    }

    void intersectLensInterface(vec3 location, vec3 direction, float diameter, float curvature, float eta, inout Ray r, inout Interaction i) {
        vec3 center = location + curvature*direction;
        vec3[2] locations;
        vec3[2] normals;
        vec2 I = intersectSphere(center, abs(curvature), r, locations, normals);

        float s = sign(curvature);
        vec2 CosTheta = -s*vec2(dot(normals[0], direction), dot(normals[1], direction));
        vec2 CosTheta2 = CosTheta*CosTheta;
        float a = 0.5*diameter/abs(curvature);
        float b = 1.0 - a*a;

        bvec2 valid = bvec2(
            I[0] > 0.0 && CosTheta[0] >= 0.0 && CosTheta2[0] >= b,
            I[1] > 0.0 && CosTheta[1] >= 0.0 && CosTheta2[1] >= b
        );

        int index = int(!valid[0]);
        float t = valid[index] ? I[index] : 0.0;

        if (0.0 < t && t < r.t) {
            r.t = t;
            i = Interaction(r.direction, locations[index], s*normals[index], Material(M_DIELECTRIC, vec4(1.0, 1.0, 1.0, eta)));
        }
    }

    void intersectLens(vec3 location, vec3 direction, float diameter, float thickness, vec2 radii, float refractiveIndex, inout Ray r, inout Interaction i) {
        intersectLensInterface(location, direction, diameter, radii[0], refractiveIndex, r, i);
        intersectLensInterface(location + thickness*direction, direction, diameter, radii[1], 1.0/refractiveIndex, r, i);
    }

    const Primitive light = SPHERE(vec3(0.0, 0.0, 5.0), 0.5, EMISSIVE(1.0, 1.0, 1.0, 15.0));

    Interaction traceRay(inout Ray r) {
        Interaction i;

        intersectObject(light, r, i);
        intersectObject(HALFSPACE(vec3(0.0), vec3(0.0, 0.0, 1.0), DIFFUSE(0.95, 0.95, 0.95)), r, i); // Floor
        intersectObject(HALFSPACE(vec3(-3.5, 0.0, 0.0), vec3(1.0, 0.0, 0.0), DIFFUSE(0.9, 0.05, 0.05)), r, i); // Left wall
        intersectObject(HALFSPACE(vec3(3.5, 0.0, 0.0), vec3(-1.0, 0.0, 0.0), DIFFUSE(0.05, 0.9, 0.05)), r, i); // Right wall
        intersectObject(HALFSPACE(vec3(0.0, 3.0, 0.0), vec3(0.0, -1.0, 0.0), DIFFUSE(0.95, 0.95, 0.95)), r, i); // Back wall
        intersectObject(HALFSPACE(vec3(0.0, 0.0, 6.0), vec3(0.0, 0.0, -1.0), DIFFUSE(0.95, 0.95, 0.95)), r, i); // Ceiling

        intersectObject(SPHERE(vec3(0.0, -2.5, 1.5), 1.5, GLASS(1.5)), r, i);
        intersectBox(vec3(-1.5, 1.0, 2.0), vec3(1.0, 3.0, 4.0), DIFFUSE(0.2, 0.2, 0.8), r, i);
        intersectBox(vec3(1.5, 1.0, 1.5), vec3(1.0, 1.0, 3.0), DIFFUSE(0.2, 0.2, 0.2), r, i);
        //intersectObject(SPHERE(vec3(0.0, 1.0, 0.75), 0.75, DIFFUSE(0.95, 0.95, 0.95)), r, i);

        return i;
    } // Scene is hard coded here.

    vec3 sampleLights(out float area, out vec3 emission, out vec3 normal) {
        emission = light.material.data.rgb; // For multiple lights, multiply by the number of lights.
        return samplePrimitive(light, normal, area);
    }

    float visibility(vec3 location, vec3 direction, float d) {
        Ray ray = Ray(location, direction, d);
        traceRay(ray);
        return ray.t < d ? 0.0 : 1.0;
    }

    float dielectricReflectance(float cosThetaI, float eta) {
        float sinThetaI = sqrt(max(0.0, 1.0 - cosThetaI*cosThetaI));
        float sinThetaT = eta*sinThetaI;
        if (sinThetaT >= 1.0) return 1.0; // Total internal reflection.

        float cosThetaT = sqrt(max(0.0, 1.0 - sinThetaT*sinThetaT));
        float etaCosThetaT = eta*cosThetaT;
        float etaCosThetaI = eta*cosThetaI;

        float rP = (cosThetaI - etaCosThetaT)/(cosThetaI + etaCosThetaT);
        float rS = (etaCosThetaI - cosThetaT)/(etaCosThetaI + cosThetaT);

        return 0.5*(rP*rP + rS*rS);
    }

    /*float conductorReflectance() {
        //TODO
    }*/

    bool isDelta(Interaction i) {
        return i.material.type == M_DIELECTRIC;
    }

    Ray sampleBSDF(Interaction i, inout vec3 attenuation) {
        vec3 oNormal = orientedNormal(i);
        Ray incoming;

        switch (i.material.type) {
            case M_EMISSIVE:
                incoming = Ray(i.location, i.outgoing, MAX_VALUE);
                // Emissive materials are treated as non reflective.
                attenuation = vec3(0.0);
                break;

            case M_DIFFUSE:
                offsetLocation(i.location, oNormal);
                incoming = Ray(i.location, sampleHemisphereCosine(oNormal), MAX_VALUE);
                // bsdf = i.material.data.rgb/PI
                // pdf = abs(dot(incoming.direction, i.normal))/PI
                // attenuation *= abs(dot(incoming.direction, i.normal))*bsdf/pdf
                attenuation *= i.material.data.rgb;
                break;

            case M_DIELECTRIC:
                float cosTheta = -dot(i.outgoing, i.normal);
                float eta = cosTheta > 0.0 ? 1.0/i.material.data[3] : i.material.data[3];
                cosTheta = abs(cosTheta);

                float f = dielectricReflectance(cosTheta, eta);
                if (IMPORTANCE_SAMPLING) {
                    if (nextRandom() < f) {
                        offsetLocation(i.location, oNormal);
                        incoming = Ray(i.location, reflect(i.outgoing, oNormal), MAX_VALUE);
                        attenuation *= i.material.data.rgb;
                    } else {
                        offsetLocation(i.location, -oNormal);
                        incoming = Ray(i.location, refract(i.outgoing, oNormal, eta), MAX_VALUE);
                        attenuation *= abs(dot(incoming.direction, i.normal))*eta*eta*i.material.data.rgb/cosTheta;
                    }
                } else {
                    if (nextRandom() < 0.5) {
                        offsetLocation(i.location, oNormal);
                        incoming = Ray(i.location, reflect(i.outgoing, oNormal), MAX_VALUE);
                        attenuation *= 2.0*f*i.material.data.rgb;
                    } else {
                        offsetLocation(i.location, -oNormal);
                        incoming = Ray(i.location, refract(i.outgoing, oNormal, eta), MAX_VALUE);
                        attenuation *= 2.0*(1.0 - f)*abs(dot(incoming.direction, i.normal))*eta*eta*i.material.data.rgb/cosTheta;
                    }
                }
                break;

            default:
                // Invalid material type.
                break;
        }

        return incoming;
    }

    vec3 evaluateBSDF(Interaction i, vec3 incoming) {
        vec3 bsdf;

        switch (i.material.type) {
            case M_DIFFUSE:
                bsdf = i.material.data.rgb/PI;
                break;
            case M_EMISSIVE:
            case M_DIELECTRIC:
                bsdf = vec3(0.0);
                break;
            default:
                // Invalid material type.
                break;
        }

        return bsdf;
    }

    vec3 sampleDirectIllumination(Interaction i) {
        offsetLocation(i.location, orientedNormal(i));

        float area;
        vec3 emission;
        vec3 normal;
        vec3 location = sampleLights(area, emission, normal);
        offsetLocation(location, dot(location - i.location, normal) > 0.0 ? -normal : normal);

        vec3 incoming = location - i.location;
        float d = length(incoming);
        incoming /= d;

        vec3 bsdf = evaluateBSDF(i, incoming);
        float cosThetaL = abs(dot(incoming, normal));
        if (bsdf == vec3(0.0) || cosThetaL == 0.0) return vec3(0.0);

        float pdf = (d*d)/(cosThetaL*area);
        return abs(dot(incoming, i.normal))*bsdf*visibility(i.location, incoming, d)*emission/pdf;
    }

    vec3 sampleGlobalIllumination(vec3 x, vec3 v) {
        vec3 radiance = vec3(0.0);
        vec3 beta = vec3(1.0);

        Ray ray = Ray(x, v, MAX_VALUE);
        bool specularBounce = false;

        for (uint b = 0u; b <= bounces && beta != vec3(0.0); ++b) {
            Interaction i = traceRay(ray);
            if (ray.t >= MAX_VALUE) break;

            if (i.material.type == M_EMISSIVE && (b == 0u || specularBounce)) {
                radiance += beta*i.material.data.rgb;
            }
            radiance += beta*sampleDirectIllumination(i);

            specularBounce = isDelta(i);

            ray = sampleBSDF(i, beta);
        }

        return radiance;
    }



    void main() {
        vec2 resolution = vec2(textureSize(pixelData, 0).xy);

        vec4[] storedData = vec4[](
            texelFetch(pixelData, ivec3(gl_FragCoord.xy, 0), 0),
            texelFetch(pixelData, ivec3(gl_FragCoord.xy, 1), 0)
        );

        vec2 screenLocation = (2.0*gl_FragCoord.xy - resolution)/min(resolution.x, resolution.y);
        IMPORTANCE_SAMPLING = screenLocation.x >= 0.0;

        highp uint nPx = uint(floor(gl_FragCoord.y)*resolution.x + floor(gl_FragCoord.x));
	    //rngOffset = nPx; // This is a piece of artwork!
        rngOffset = intHash(nPx);
        rngIndex = nSample == 0u ? 0u : floatBitsToUint(storedData[0][3]);

        vec2 pxLocation = vec2(nextRandom(), nextRandom()) - vec2(0.5);
        vec2 px = (2.0*(gl_FragCoord.xy + pxLocation) - resolution)/min(resolution.x, resolution.y);

        //vec2 px = (2.0*gl_FragCoord.xy - resolution)/min(resolution.x, resolution.y);

        vec3 rayDirection = normalize(px.x*cameraX + px.y*cameraY - focalLength*cameraZ);

        outputData[0].rgb = storedData[0].rgb + sampleGlobalIllumination(cameraLocation, rayDirection);
        outputData[0][3] = uintBitsToFloat(rngIndex);
        outputData[1][0] = storedData[1][0] + 1.0; // Update sample count.
    }
`;
