// Original C++ code from the minibook
//   Ray Tracing in One Weekend
//   Peter Shirley: http://in1weekend.blogspot.com 
//

// Chapter 12: Where Next?

class Camera {
	constructor(lookFrom, lookAt, vup, vfov, aspect, aperture, focus_dist) {
		const theta = vfov * Math.PI / 180.0;
		const half_height = Math.tan(theta / 2);
		const half_width  = aspect * half_height;
		const w = unitVec(subVec(lookFrom, lookAt));	// Backward
		const u = unitVec(crossVec(vup, w));			// Right
		const v = crossVec(w, u);						// Up
		const su = scaleVec(half_width*focus_dist, u);
		const sv = scaleVec(half_height*focus_dist, v);
		const sw = scaleVec(focus_dist, w);
		this.origin = lookFrom;
		this.lower_left_corner = subVec(lookFrom, addVec(sw,addVec(su,sv)));
		this.horizontal = scaleVec(2, su);
		this.vertical = scaleVec(2, sv);
		this.u = u;
		this.v = v;
		this.w = w;
		this.lens_radius = aperture / 2;
	}
	getRay = function(s, t) {
		const rd = scaleVec(this.lens_radius,randomInUnitDisk());
		const offset = addVec(scaleVec(rd[0],this.u),scaleVec(rd[1],this.v));
		const p = addVec(addVec(scaleVec(s, this.horizontal), scaleVec(t, this.vertical)), this.lower_left_corner);
		const o = addVec(this.origin,offset);

		return new Ray(o, subVec(p, o));
	}
}

class Ray {
	constructor(org, dir) {
		this.origin = org;
		this.direction = dir;
	}
	pointAtParameter(t) {
		return [this.origin[0] + this.direction[0] * t,
		this.origin[1] + this.direction[1] * t,
		this.origin[2] + this.direction[2] * t];
	}
}

class Sphere {
	constructor(center, radius, material) {
		this.center = center;
		this.radius = radius;
	    this.material = material;
	}
	hit(r, tmin, tmax, rec) {
		const ocn = subVec(r.origin, this.center);
		const dir = r.direction;
		const a = dotVec(dir, dir);
		const b = 2 * dotVec(ocn, dir);
		const c = dotVec(ocn, ocn) - this.radius * this.radius;
		const discriminant = b * b - 4 * a * c;
		if (discriminant > 0) {
			const sqrtdisc = Math.sqrt(discriminant);
			let temp = (-b - sqrtdisc) / (2 * a);
			if (temp < tmax && temp > tmin) {
				rec.t = temp;
				rec.p = r.pointAtParameter(temp);
				// Unit vector from center to p (i.e., normal)
				rec.normal = scaleVec(1.0 / this.radius, subVec(rec.p, this.center));
				rec.material = this.material;
				return true;
			}
			temp = (-b + sqrtdisc) / (2 * a);
			if (temp < tmax && temp > tmin) {
				rec.t = temp;
				rec.p = r.pointAtParameter(temp);
				// Unit vector from center to p (i.e., normal)
				rec.normal = scaleVec(1.0 / this.radius, subVec(rec.p, this.center));
				rec.material = this.material;
				return true;
			}
		}
		return false;
	}
}

class Lambertian {
	constructor(albedo) {
		this.albedo = albedo;
	}
	scatter(r, rec, sctr) {
		const target = addVec(addVec(rec.p, rec.normal), randomInUnitSphere());
		sctr.scattered = new Ray(rec.p, subVec(target, rec.p));
		sctr.attenuation = this.albedo;
		return true;
	}
}

class Metal {
	constructor(albedo) {
		this.albedo = albedo;
	}
	scatter(r, rec, sctr) {
		const reflected = reflect(unitVec(r.direction), rec.normal);
		sctr.scattered = new Ray(rec.p, reflected);
		sctr.attenuation = this.albedo;
		return (dotVec(reflected, rec.normal) > 0);
	}
}

class Dielectric {
	constructor(ri) {
		this.refindex = ri;
	}
	scatter(r, rec, sctr) {
		let outward_normal;
		let cosine;
		let ni_over_nt;
		let refracted = [0, 0, 0];
		let reflected = reflect(r.direction, rec.normal);
		let reflect_prob;
		sctr.attenuation = [1, 1, 1];
		if (dotVec(r.direction, rec.normal) > 0) {
			outward_normal = [-rec.normal[0], -rec.normal[1], -rec.normal[2]];
			ni_over_nt = this.refindex;
			cosine = this.refindex * dotVec(r.direction, rec.normal) / lenVec(r.direction);
		} else {
			outward_normal = rec.normal;
			ni_over_nt = 1.0 / this.refindex;
			cosine = -dotVec(r.direction, rec.normal) / lenVec(r.direction);
		}

		if (refract(r.direction, outward_normal, ni_over_nt, refracted))
			reflect_prob = schlick(cosine, this.refindex);

		else
			reflect_prob = 1;

		if (Math.random() < reflect_prob)
			sctr.scattered = new Ray(rec.p, reflected);

		else
			sctr.scattered = new Ray(rec.p, refracted);

		return true;
	}
}

function schlick(cosine, refindex) {
    var r0 = (1 - refindex) / (1 + refindex);
    r0 = r0 * r0;
    return r0 + (1 - r0) * Math.pow(1 - cosine, 5);
}

function reflect(v,n) {
    // Should v be a unit vector?
    return subVec(v, scaleVec(2 * dotVec(v,n), n));
}

function refract(v,n,ni_over_nt,refracted) {
    // dt is cos(theta_i) < 0, guaranteed by callee
	const uv = unitVec(v);
    const dt = dotVec(uv,n);
    const discriminant = 1.0 - ni_over_nt * ni_over_nt * (1 - dt*dt);
    if (discriminant > 0) {
        const sqrtdisc = Math.sqrt(discriminant);
        refracted[0] = ni_over_nt * (uv[0] - n[0] * dt) - n[0] * sqrtdisc;
        refracted[1] = ni_over_nt * (uv[1] - n[1] * dt) - n[1] * sqrtdisc;
        refracted[2] = ni_over_nt * (uv[2] - n[2] * dt) - n[2] * sqrtdisc;
        return true;
    }
    return false;
}

// Vector functions
function lenVec(v) { return Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]); }
function dotVec(u,v) { return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]; }
function scaleVec(s,v) { return [s*v[0], s*v[1], s*v[2]]; }
function addVec(u,v) { return [u[0]+v[0], u[1]+v[1], u[2]+v[2]]; }
function subVec(u,v) { return [u[0]-v[0], u[1]-v[1], u[2]-v[2]]; }
function mulVec(u,v) { return [u[0]*v[0], u[1]*v[1], u[2]*v[2]]; }

function crossVec(u,v) {
    return [u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0]];
}

function affineVec(u,v,s) {
    return [(1.0 - s) * u[0] + s * v[0],
            (1.0 - s) * u[1] + s * v[1],
            (1.0 - s) * u[2] + s * v[2]];
}

function unitVec(v) {
    const len = Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (len == 0) return [0, 0, 0];
    return [v[0]/len, v[1]/len, v[2]/len];
}

function worldHit(world, r, tmin, tmax, rec) {
    let hitAnything = false;
    let temprec = {};
    let closest = tmax;
    for (let i = 0; i < world.length; ++i) {
        if (world[i].hit(r,tmin,closest,temprec)) {
            hitAnything = true;
			rec.t = temprec.t;
			rec.p = temprec.p;
			rec.normal = temprec.normal;
			rec.material = temprec.material;
			closest = temprec.t;
        }
    }
    return hitAnything;
}

function randomInUnitDisk() {
    let x, y;

    do {
        x = Math.random() * 2.0 - 1.0;
        y = Math.random() * 2.0 - 1.0;
    } while (x*x + y*y >= 1.0);

    return [x, y, 0];
}

function randomInUnitSphere() {
    let x, y, z;

    do {
        x = Math.random() * 2.0 - 1.0;
        y = Math.random() * 2.0 - 1.0;
        z = Math.random() * 2.0 - 1.0;
    } while (x*x + y*y + z*z >= 1.0);

    return [x, y, z];
}

function color(r,world,depth) {
    // Hit record:
    //    p - point where ray intercepts scene object
    //    t - parameter corresponding to p
    //    normal - normal vector at p
    //    material - material of the scene object
    let rec = {};

    // Scatter record:
    //     attenuation - color coefficients
    //     scattered - ray in scattered direction
    let sctr = {};

    if (worldHit(world, r, 0.001, 99999.9, rec)) {
        if (depth < 50 && rec.material.scatter(r,rec,sctr)) {
            const attenuation = sctr.attenuation;
            const scattered   = sctr.scattered;
            const col = color(scattered,world,depth+1);
            return mulVec(attenuation,col);
        } else {
            return [0,0,0];
        }
    }

    // Ray didn't hit any object, return ambient gradient
    // Whitish at the bottom, bluish at the top
    const uni = unitVec(r.direction);
    const t = 0.5 * (uni[1] + 1);
    return affineVec([1,1,1], [0.5, 0.7, 1.0], t);
}

function createWorld() {
    var grid = 7;	// Try smaller values for faster rendering
    var a,b;
    var center,choose;
    var rad = 0.2;
    var base = [4, 0.2, 0];
    var list = [new Sphere([0,-1000,0],1000, new Lambertian([0.5,0.5,0.5]))];

    for (a = -grid; a < grid; ++a) {
        for (b = -grid; b < grid; ++b) {
            choose = Math.random();
            center = [a+0.9*Math.random(), 0.2, b+0.9*Math.random()];
            if (lenVec(subVec(center, base)) > 0.9) {
                if (choose < 0.8) {  // Diffuse
//                   list.push(new Sphere(center, rad, new Lambertian([Math.random()*Math.random,Math.random()*Math.random(),Math.random()*Math.random()])));
                   list.push(new Sphere(center, rad, new Lambertian([Math.random(),0.6*Math.random(),0.6*Math.random()])));
				} else if (choose < 0.95) { // Metal
                    list.push(new Sphere(center, rad, new Metal([0.5*(1+Math.random()),0.5*(1+Math.random()),0.5*(1+Math.random())])));
                } else { // Glass
                    list.push(new Sphere(center, rad, new Dielectric(1.5)));
                }
            }
        }
    }

    list.push(new Sphere([0,1,0], 1.0, new Dielectric(1.5)));
    list.push(new Sphere([-4,1,0], 1.0, new Lambertian([0.4, 0.2, 0.1])));
    list.push(new Sphere([4,1,0], 1.0, new Metal([0.7, 0.6, 0.5])));

    return list;
}

function main() {
	// JavaScript is too slow for us to take 100 samples per pixel,
	// so we use a more modest number here. Higher values produce
	// better pictures but take longer to render.
	const numSamples = 50;	
	const adjust = 255 / numSamples;

	const canvas = document.getElementById('myCanvas');
	const ctx = canvas.getContext('2d');

	// Canvas size
	const cv_lx = canvas.width;
	const cv_ly = canvas.height;

	// Maximum x and y coordinates
	const cv_mx = cv_lx - 1;
	const cv_my = cv_ly - 1;

	// Camera
	const lookFrom = [13,2,3];
	const lookAt = [0,0,0];
	const dist_to_focus = lenVec(subVec(lookFrom,lookAt));
	const aperture = 0.1;
	const camera = new Camera(lookFrom,lookAt,[0,1,0], 20, cv_lx / cv_ly, aperture, dist_to_focus);

	// Scene objects
	const world = createWorld();

	// imageData.data is a Uint8ClampedArray[] so any set
	// with a float value < 0 or > 255 is automatically clamped
	// at these limits.
	const imageData = ctx.getImageData(0, 0, cv_lx, cv_ly);
	const pixels = imageData.data;

	for (let cv_y = cv_my; cv_y >= 0; --cv_y) {
		const b = (cv_my - cv_y) * cv_lx;	// Base pixel offset for this line
	    for (let cv_x = 0; cv_x <= cv_mx; ++cv_x) {
			let col = [0,0,0];
			for (let s = 1; s <= numSamples; ++s) {
				const u = (cv_x + Math.random()) / cv_lx;
				const v = (cv_y + Math.random()) / cv_ly;
				col = addVec(col, color(camera.getRay(u,v), world, 0));
			}

	        const offset = (b + cv_x) * 4;
			// With gamma 2 correction (col^(1/2) = sqrt(col))
	        pixels[offset  ] = Math.sqrt(col[0]/numSamples) * 255;
	        pixels[offset+1] = Math.sqrt(col[1]/numSamples) * 255;
	        pixels[offset+2] = Math.sqrt(col[2]/numSamples) * 255;
	        pixels[offset+3] = 255;
			// Without gamma correction
	        // pixels[offset  ] = col[0] * adjust;
	        // pixels[offset+1] = col[1] * adjust;
	        // pixels[offset+2] = col[2] * adjust;
	        // pixels[offset+3] = 255;
	    }
	}
	
	ctx.putImageData(imageData, 0, 0);
}

const isNode = typeof process !== "undefined" && process?.versions?.node;

main();



