// Original C++ code from the minibook
//   Ray Tracing in One Weekend
//   Peter Shirley: http://in1weekend.blogspot.com 
//

// Chapter 5: Surface Normals and Multiple Objects

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

const sphere_center = [0,0,-1];
const sphere_radius = 0.5;

function main() {
	const lower_left_corner = [-2, -1, -1];
	const origin = [0,0,0];
	const horizontal = [4,0,0];
	const vertical = [0,2,0];
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
		const v = cv_y / cv_my;				// Vertical scale (0 <= v <= 1)
		const b = (cv_my - cv_y) * cv_lx;	// Base pixel offset for this line
	    for (let cv_x = 0; cv_x <= cv_mx; ++cv_x) {
			const u = cv_x / cv_mx;			// Horizontal scale (0 <= u <= 1)
			const p = addVec(scaleVec(u, horizontal), scaleVec(v, vertical));
			const r = new Ray(origin, addVec(lower_left_corner, p));
	        const col = Color(r, world);

	        const offset = (b + cv_x) * 4;
	        pixels[offset  ] = col[0] * 255;
	        pixels[offset+1] = col[1] * 255;
	        pixels[offset+2] = col[2] * 255;;
	        pixels[offset+3] = 255;
	    }
	}
	
	ctx.putImageData(imageData, 0, 0);
}

main();
