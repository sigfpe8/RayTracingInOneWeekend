// Original C++ code from the minibook
//   Ray Tracing in One Weekend
//   Peter Shirley: http://in1weekend.blogspot.com 
//

// Chapter 3: Rays, a Simple Camera and Background

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

// Vector functions
function lenVec(v) { return Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]); }
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

function Color(r) {
    const uni = unitVec(r.direction);
    const t = 0.5 * (uni[1] + 1);
	// Gradient from white [1,1,1] to light blue [0.5, 0.7,1]
    return affineVec([1,1,1],[0.5, 0.7, 1.0], t);
}

function main() {
	const lower_left_corner = [-2, -1, -1];
	const origin = [0,0,0];
	const horizontal = [4,0,0];
	const vertical = [0,2,0];
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
	        const col = Color(r);

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
