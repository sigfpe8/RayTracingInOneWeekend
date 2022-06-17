// Original C++ code from the minibook
//   Ray Tracing in One Weekend
//   Peter Shirley: http://in1weekend.blogspot.com 
//

// Chapter 6: Antialiasing

class Camera {
	constructor() {
		this.lower_left_corner = [-2, -1, -1];
		this.origin = [0, 0, 0];
		this.horizontal = [4, 0, 0];
		this.vertical = [0, 2, 0];
	}
	getRay(u, v) {
		const p = addVec(scaleVec(u, this.horizontal), scaleVec(v, this.vertical));

		return new Ray(this.origin, addVec(this.lower_left_corner, p));
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
	constructor(center, radius) {
		this.center = center;
		this.radius = radius;
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
				return true;
			}
			temp = (-b + sqrtdisc) / (2 * a);
			if (temp < tmax && temp > tmin) {
				rec.t = temp;
				rec.p = r.pointAtParameter(temp);
				// Unit vector from center to p (i.e., normal)
				rec.normal = scaleVec(1.0 / this.radius, subVec(rec.p, this.center));
				return true;
			}
		}
		return false;
	}
}

// Vector functions
function lenVec(v) { return Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]); }
function dotVec(u,v) { return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]; }
function scaleVec(s,v) { return [s*v[0], s*v[1], s*v[2]]; }
function addVec(u,v) { return [u[0]+v[0], u[1]+v[1], u[2]+v[2]]; }
function subVec(u,v) { return [u[0]-v[0], u[1]-v[1], u[2]-v[2]]; }

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
			closest = temprec.t;
        }
    }
    return hitAnything;
}

function Color(r,world) {
    let rec = {};

    if (worldHit(world, r, 0.0, 99999.9, rec)) {
        const nor = rec.normal;
        return [(nor[0]+1)*0.5, (nor[1]+1)*0.5, (nor[2]+1)*0.5];
    }

    var uni = unitVec(r.direction);
    var t = 0.5 * (uni[1] + 1);
    return affineVec([0.5, 0.7, 1.0], [1,1,1], t);
}

function main() {
	// JavaScript is too slow for us to take 100 samples per pixel,
	// so we use a more modest number here. Higher values produce
	// better pictures but take longer to render.
	const numSamples = 40;
	const adjust = 255 / numSamples;
	const camera = new Camera();
	const world = [new Sphere([0,0,-1],0.5), new Sphere([0,-100.5,-1],100)];

	const canvas = document.getElementById('myCanvas');
	const ctx = canvas.getContext('2d');

	// Canvas size
	const cv_lx = canvas.width;
	const cv_ly = canvas.height;

	// Maximum x and y coordinates
	const cv_mx = cv_lx - 1;
	const cv_my = cv_ly - 1;

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
				col = addVec(col, Color(camera.getRay(u,v), world));
			}

	        const offset = (b + cv_x) * 4;
	        pixels[offset  ] = col[0] * adjust;
	        pixels[offset+1] = col[1] * adjust;
	        pixels[offset+2] = col[2] * adjust;
	        pixels[offset+3] = 255;
	    }
	}
	
	ctx.putImageData(imageData, 0, 0);
}

main();
