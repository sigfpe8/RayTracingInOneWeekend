// Original C++ code from the minibook
//   Ray Tracing in One Weekend
//   Peter Shirley: http://in1weekend.blogspot.com 
//

// Chapter 1: Output an Image

function main() {
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

	for (let cv_y = 0; cv_y <= cv_my; ++cv_y) {
	    const clr_g = (cv_y / cv_my) * 255;
	    for (let cv_x = 0; cv_x <= cv_mx; ++cv_x) {

	        const offset = ((cv_my - cv_y) * cv_lx + cv_x) * 4;

	        const clr_r = (cv_x / cv_mx) * 255;

	        pixels[offset  ] = clr_r;
	        pixels[offset+1] = clr_g;
	        pixels[offset+2] = 50;
	        pixels[offset+3] = 255;
	    }
	}
	
	ctx.putImageData(imageData, 0, 0);
}

main();
