class PathTracer {

    constructor(gl) {
        this.MODE_SAMPLE = 0;
        this.MODE_RENDER = 1;

        this.gl = gl;
        this.glUtils = WebGL2Utils.from(gl);

        this.framebuffer = gl.createFramebuffer();
        gl.bindFramebuffer(gl.DRAW_FRAMEBUFFER, this.framebuffer);
        gl.drawBuffers([gl.COLOR_ATTACHMENT0, gl.COLOR_ATTACHMENT1]);

        this.framebufferComplete = false;
        this.pixelDataInvalid = false;
        this.inputTextureIndex = 0;

        this.samplerProgram = new SamplerProgram(gl);
        this.rendererProgram = new RendererProgram(gl);

        this.nSample = 0;

        this.camera = {
            target: [0, 0, 2.5],
            distance: 9.5,
            theta: 0.5 * Math.PI,
            phi: 0,
            focalLength: 2,
        };

        this.updateCamera();
        //this.samplerProgram.uniforms.bounces = 10; //TODO: This should be configurable.
        this.samplerProgram.uniforms.bounces = 10;
    }

    updateResolution(resolution) {
        const {css: cssResolution, buffer: bufferResolution, changed} = this.glUtils.updateResolution(resolution);
        this.cssResolution = cssResolution;
        this.bufferResolution = bufferResolution;

        if (changed || this.pixelDataTextures === undefined) {
            this.pixelDataTextures = [this.gl.createTexture(), this.gl.createTexture()];

            this.pixelDataTextures.forEach((texture, index) => {
                this.gl.activeTexture(this.gl.TEXTURE0 + index);
                this.gl.bindTexture(this.gl.TEXTURE_2D_ARRAY, texture);
                this.gl.texParameteri(this.gl.TEXTURE_2D_ARRAY, this.gl.TEXTURE_MAG_FILTER, this.gl.NEAREST);
                this.gl.texParameteri(this.gl.TEXTURE_2D_ARRAY, this.gl.TEXTURE_MIN_FILTER, this.gl.NEAREST);
                this.gl.texStorage3D(this.gl.TEXTURE_2D_ARRAY, 1, this.gl.RGBA32F, ...this.bufferResolution, 2);
            });

            this.nSample = 0;

            return true;
        }

        return false;
    }

    updateCamera({
        target = this.camera.target,
        distance = this.camera.distance,
        theta = this.camera.theta,
        phi = this.camera.phi,
        focalLength = this.camera.focalLength
    } = {}) {
        this.camera.target = target;
        this.camera.distance = Math.max(1, distance);
        this.camera.theta = PathTracer.clamp(theta, 0, Math.PI);
        this.camera.phi = PathTracer.mod(phi, 2 * Math.PI);
        this.camera.focalLength = Math.max(0, focalLength);

        const sinTheta = Math.sin(this.camera.theta);
        const cosTheta = Math.cos(this.camera.theta);
        const sinPhi = Math.sin(this.camera.phi);
        const cosPhi = Math.cos(this.camera.phi);

        const cameraX = [cosPhi, sinPhi, 0];
        const cameraY = [-cosTheta * sinPhi, cosTheta * cosPhi, sinTheta];
        const cameraZ = [sinTheta * sinPhi, -sinTheta * cosPhi, cosTheta];

        const location = [
            this.camera.target[0] + this.camera.distance * cameraZ[0],
            this.camera.target[1] + this.camera.distance * cameraZ[1],
            this.camera.target[2] + this.camera.distance * cameraZ[2]
        ];

        this.samplerProgram.setCamera(location, this.camera.focalLength, cameraX, cameraY, cameraZ);
        this.clear();
    }

    setMode(mode) {
        if (this.mode === mode) return;

        switch (mode) {
            case this.MODE_SAMPLE:
                this.gl.bindFramebuffer(this.gl.DRAW_FRAMEBUFFER, this.framebuffer);
                this.samplerProgram.useProgram();
                break;

            case this.MODE_RENDER:
                this.gl.bindFramebuffer(this.gl.DRAW_FRAMEBUFFER, null);
                this.rendererProgram.useProgram();
                this.rendererProgram.uniforms.pixelData = this.inputTextureIndex; // Last output as input
                this.rendererProgram.uniforms.samples = this.nSample;
                this.rendererProgram.updateUniforms();
                break;

            default:
            // Will never happen.
        }

        this.mode = mode;
    }

    sample() {
        if (this.pixelDataTextures === undefined) this.updateResolution();
        this.setMode(this.MODE_SAMPLE);

        if (this.pixelDataInvalid) {
            this.gl.clear(this.gl.COLOR_BUFFER_BIT);
            this.pixelDataInvalid = false;
        }

        const outputTexture = this.pixelDataTextures[1 - this.inputTextureIndex];
        this.gl.framebufferTextureLayer(this.gl.DRAW_FRAMEBUFFER, this.gl.COLOR_ATTACHMENT0, outputTexture, 0, 0);
        this.gl.framebufferTextureLayer(this.gl.DRAW_FRAMEBUFFER, this.gl.COLOR_ATTACHMENT1, outputTexture, 0, 1);
        this.framebufferComplete = true;

        this.samplerProgram.uniforms.pixelData = this.inputTextureIndex;
        this.samplerProgram.uniforms.nSample = this.nSample;
        this.samplerProgram.updateUniforms();
        this.gl.drawArrays(this.gl.TRIANGLE_STRIP, 0, 4);

        this.nSample += 1;
        this.inputTextureIndex = 1 - this.inputTextureIndex;
    }

    render() {
        if (this.pixelDataTextures === undefined) this.updateResolution();
        this.setMode(this.MODE_RENDER);
        this.gl.drawArrays(this.gl.TRIANGLE_STRIP, 0, 4);
    }

    clear() {
        if (this.framebufferComplete) this.pixelDataInvalid = true;
        this.nSample = 0;
    }



    static clamp(x, lower, upper) {
        return Math.min(upper, Math.max(lower, x));
    }

    static mod(a, b) {
        return a - b * Math.floor(a / b);
    } // Mathematical modulo operation yields the correct sign.

}
