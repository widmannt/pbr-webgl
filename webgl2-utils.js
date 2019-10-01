class WebGL2Utils {

    constructor() {
        throw TypeError("Illegal constructor");
    }

    /**
     * Create and compile a shader.
     *
     * @param {string} source The GLSL source code for the shader.
     * @param {number} type The type of shader, VERTEX_SHADER or FRAGMENT_SHADER.
     *
     * @return {!WebGLShader} The shader.
     */
    compileShader(source, type) {
        const shader = this.gl.createShader(type);

        this.gl.shaderSource(shader, source);
        this.gl.compileShader(shader);

        if (!this.gl.getShaderParameter(shader, this.gl.COMPILE_STATUS)) {
            const error = Error("Could not compile shader:\n" + this.gl.getShaderInfoLog(shader));
            this.gl.deleteShader(shader);
            throw error;
        }

        return shader;
    }

    /**
     * Create and compile a shader from the content of a script element.
     *
     * @param {string} shaderId The id of the script element.
     *
     * @return {!WebGLShader} The shader.
     */
    compileShaderFromScript(shaderId) {
        const shaderElement = document.getElementById(shaderId);
        if (!shaderElement) {
            throw Error("Could not find script element: " + shaderId);
        }

        const shaderSource = shaderElement.firstChild.data.trim();

        let shaderType;
        switch (shaderElement.type) {
            case "x-shader/x-vertex":
                shaderType = this.gl.VERTEX_SHADER;
                break;

            case "x-shader/x-fragment":
                shaderType = this.gl.FRAGMENT_SHADER;
                break;

            default:
                throw Error("Unknown shader type: " + shaderElement.type);
        }

        return this.compileShader(shaderSource, shaderType);
    }

    /**
     * Create a program from vertex shader and fragment shader.
     *
     * @param {!WebGLShader} vertShader The vertex shader.
     * @param {!WebGLShader} fragShader The fragment shader.
     *
     * @return {!WebGLProgram} The program.
     */
    createProgram(vertShader, fragShader) {
        const program = this.gl.createProgram();

        this.gl.attachShader(program, vertShader);
        this.gl.attachShader(program, fragShader);

        this.gl.linkProgram(program);

        if (!this.gl.getProgramParameter(program, this.gl.LINK_STATUS)) {
            const error = Error("Could not create program:\n" + this.gl.getProgramInfoLog(program));
            this.gl.deleteProgram(program);
            throw error;
        }

        return program;
    }

    /**
     * Create a program from GLSL source code.
     *
     * @param {string} vertShaderSource The GLSL source code for the vertex shader.
     * @param {string} fragShaderSource The GLSL source code for the fragment shader.
     *
     * @return {!WebGLProgram} The program.
     */
    createProgramFromSources(vertShaderSource, fragShaderSource) {
        const vertexShader = this.compileShader(vertShaderSource, this.gl.VERTEX_SHADER);
        const fragmentShader = this.compileShader(fragShaderSource, this.gl.FRAGMENT_SHADER);

        return this.createProgram(vertexShader, fragmentShader);
    }

    /**
     * Create a program from the contents of two script elements.
     *
     * @param {string} vertShaderId The vertex shader element id.
     * @param {string} fragShaderId The fragment shader element id.
     *
     * @return {!WebGLProgram} The program.
     */
    createProgramFromScripts(vertShaderId, fragShaderId) {
        const vertexShader = this.compileShaderFromScript(vertShaderId);
        const fragmentShader = this.compileShaderFromScript(fragShaderId);

        return this.createProgram(vertexShader, fragmentShader);
    }

    /**
     * Resize the drawing buffer of the canvas to match device resolution.
     *
     * @param {[number, number]} [resolution] The resolution in device pixels.
     *
     * @return {{css: [number, number], buffer: [number, number], changed: boolean}} The updated resolution values.
     */
    updateResolution(resolution) {
        const pixelRatio = window.devicePixelRatio || 1;

        if (typeof resolution === "undefined") {
            const canvasRect = this.gl.canvas.getBoundingClientRect();
            resolution = [Math.floor(canvasRect.width * pixelRatio), Math.floor(canvasRect.height * pixelRatio)];
        }

        let changed = this.gl.canvas.width !== resolution[0] || this.gl.canvas.height !== resolution[1];
        if (changed) {
            this.gl.canvas.width = resolution[0];
            this.gl.canvas.height = resolution[1];
        }

        // Update css size to improve alignment with device pixels
        const cssResolution = [resolution[0] / pixelRatio, resolution[1] / pixelRatio];
        this.gl.canvas.style.width = cssResolution[0] + "px";
        this.gl.canvas.style.height = cssResolution[1] + "px";

        const bufferResolution = [this.gl.drawingBufferWidth, this.gl.drawingBufferHeight];
        if (changed) {
            this.gl.viewport(0, 0, ...bufferResolution);
        }

        return {css: cssResolution, buffer: bufferResolution, changed};
    }

    /**
     * Find the name of a GLenum value.
     *
     * @param {GLenum} value The enum value.
     *
     * @return {string} The corresponding name.
     */
    glEnumToString(value) {
        for (let key in this.gl) {
            if (this.gl[key] === value) return key;
        }

        return `0x${value.toString(16)}`;
    }

    /**
     * Print uniform block info to console.
     *
     * @param {!WebGLProgram} program The program.
     * @param {string} name The name of the uniform block.
     */
    printUniformBlockInfo(program, name) {
        const index = this.gl.getUniformBlockIndex(program, name);
        const dataSize = this.gl.getActiveUniformBlockParameter(program, index, this.gl.UNIFORM_BLOCK_DATA_SIZE);
        const activeUniforms = this.gl.getActiveUniformBlockParameter(program, index, this.gl.UNIFORM_BLOCK_ACTIVE_UNIFORMS);

        const indices = this.gl.getActiveUniformBlockParameter(program, index, this.gl.UNIFORM_BLOCK_ACTIVE_UNIFORM_INDICES);
        const offsets = this.gl.getActiveUniforms(program, indices, this.gl.UNIFORM_OFFSET);

        const uniforms = {};
        for (let i = 0; i < activeUniforms; ++i) {
            const uniform = this.gl.getActiveUniform(program, indices[i]);
            uniforms[uniform.name] = {
                type: this.glEnumToString(uniform.type),
                size: uniform.size,
                offset: offsets[i]
            };
        }

        // BUG: Groups are not filtered by log level in Firefox and Chrome.
        console.group(this.printUniformBlockInfo.name + ": " + name);
        console.log("data size: " + dataSize);
        console.log("active uniforms: " + activeUniforms);
        console.table(uniforms);
        console.groupEnd();
    }



    /**
     * Provide an instance of WebGL2Utils for the specified WebGL2 context.
     *
     * @param {!WebGL2RenderingContext} gl The WebGL2 context.
     *
     * @return {!WebGL2Utils} An instance associated with the specified WebGL2 context.
     */
    static from(gl) {
        if (!(gl instanceof WebGL2RenderingContext)) {
            throw Error("Can't create instance: No WebGL2RenderingContext provided");
        }

        if (typeof WebGL2Utils.instances === "undefined") {
            WebGL2Utils.instances = new WeakMap();
            return createInstance();
        }

        return WebGL2Utils.instances.get(gl) || createInstance();

        function createInstance() {
            const instance = Object.create(WebGL2Utils.prototype);
            WebGL2Utils.instances.set(gl, instance);

            instance.gl = gl;

            instance.UniformBuffer = class {
                constructor(blockBinding, usage) {
                    instance.buffer = gl.createBuffer();
                    instance.blockBinding = blockBinding;
                    instance.usage = usage;
                }

                update(data) {
                    if (typeof instance.data === "undefined") {
                        instance.data = new Float32Array(data);
                    } else {
                        instance.data.set(data);
                    }

                    gl.bindBuffer(gl.UNIFORM_BUFFER, instance.buffer);
                    gl.bufferData(gl.UNIFORM_BUFFER, instance.data, instance.usage);
                    gl.bindBufferBase(gl.UNIFORM_BUFFER, instance.blockBinding, instance.buffer);
                }
            };

            return instance;
        }
    }

}

class WebGL2ProgramContainer {

    constructor(gl, program) {
        this.gl = gl;
        this.program = program;
        this.pendingUniforms = new Map();

        const uniforms = Object.create(null);

        const activeUniforms = gl.getProgramParameter(program, gl.ACTIVE_UNIFORMS);
        for (let i = 0; i < activeUniforms; ++i) {
            const uniformInfo = gl.getActiveUniform(program, i);
            const uniformSetter = uniformSetterFor(uniformInfo);

            Object.defineProperty(uniforms, uniformInfo.name, {
                set: value => this.pendingUniforms.set(uniformSetter, value),
                get: () => this.pendingUniforms.get(uniformSetter)
            });
        }

        Object.defineProperty(this, "uniforms", {
            set: value => Object.assign(uniforms, value),
            get: () => uniforms
        });

        function uniformSetterFor(uniform) {
            const isArray = uniform.size > 1;
            const location = gl.getUniformLocation(program, uniform.name);
            return {
                [gl.FLOAT]: isArray ? gl.uniform1fv : gl.uniform1f,
                [gl.FLOAT_VEC2]: gl.uniform2fv,
                [gl.FLOAT_VEC3]: gl.uniform3fv,
                [gl.FLOAT_VEC4]: gl.uniform4fv,
                [gl.INT]: isArray ? gl.uniform1iv : gl.uniform1i,
                [gl.INT_VEC2]: gl.uniform2iv,
                [gl.INT_VEC3]: gl.uniform3iv,
                [gl.INT_VEC4]: gl.uniform4iv,
                [gl.UNSIGNED_INT]: isArray ? gl.uniform1uiv : gl.uniform1ui,
                [gl.UNSIGNED_INT_VEC2]: gl.uniform2uiv,
                [gl.UNSIGNED_INT_VEC3]: gl.uniform3uiv,
                [gl.UNSIGNED_INT_VEC4]: gl.uniform4uiv,
                [gl.BOOL]: isArray ? gl.uniform1iv : gl.uniform1i,
                [gl.BOOL_VEC2]: gl.uniform2iv,
                [gl.BOOL_VEC3]: gl.uniform3iv,
                [gl.BOOL_VEC4]: gl.uniform4iv,
                [gl.FLOAT_MAT2]: gl.uniformMatrix2fv,
                [gl.FLOAT_MAT3]: gl.uniformMatrix3fv,
                [gl.FLOAT_MAT4]: gl.uniformMatrix4fv,
                [gl.FLOAT_MAT2x3]: gl.uniformMatrix2x3fv,
                [gl.FLOAT_MAT2x4]: gl.uniformMatrix2x4fv,
                [gl.FLOAT_MAT3x2]: gl.uniformMatrix3x2fv,
                [gl.FLOAT_MAT3x4]: gl.uniformMatrix3x4fv,
                [gl.FLOAT_MAT4x2]: gl.uniformMatrix4x2fv,
                [gl.FLOAT_MAT4x3]: gl.uniformMatrix4x3fv,
                [gl.SAMPLER_2D]: isArray ? gl.uniform1iv : gl.uniform1i,
                [gl.SAMPLER_2D_ARRAY]: isArray ? gl.uniform1iv : gl.uniform1i,
                [gl.SAMPLER_2D_ARRAY_SHADOW]: isArray ? gl.uniform1iv : gl.uniform1i,
                [gl.SAMPLER_2D_SHADOW]: isArray ? gl.uniform1iv : gl.uniform1i,
                [gl.SAMPLER_3D]: isArray ? gl.uniform1iv : gl.uniform1i,
                [gl.SAMPLER_CUBE]: isArray ? gl.uniform1iv : gl.uniform1i,
                [gl.SAMPLER_CUBE_SHADOW]: isArray ? gl.uniform1iv : gl.uniform1i,
                [gl.INT_SAMPLER_2D]: isArray ? gl.uniform1iv : gl.uniform1i,
                [gl.INT_SAMPLER_2D_ARRAY]: isArray ? gl.uniform1iv : gl.uniform1i,
                [gl.INT_SAMPLER_3D]: isArray ? gl.uniform1iv : gl.uniform1i,
                [gl.INT_SAMPLER_CUBE]: isArray ? gl.uniform1iv : gl.uniform1i,
                [gl.UNSIGNED_INT_SAMPLER_2D]: isArray ? gl.uniform1iv : gl.uniform1i,
                [gl.UNSIGNED_INT_SAMPLER_2D_ARRAY]: isArray ? gl.uniform1iv : gl.uniform1i,
                [gl.UNSIGNED_INT_SAMPLER_3D]: isArray ? gl.uniform1iv : gl.uniform1i,
                [gl.UNSIGNED_INT_SAMPLER_CUBE]: isArray ? gl.uniform1iv : gl.uniform1i,
            }[uniform.type].bind(gl, location);
        }
    }

    useProgram() {
        this.gl.useProgram(this.program);
    }

    updateUniforms() {
        for (let [setter, value] of this.pendingUniforms) {
            setter(value);
        }
        this.pendingUniforms.clear();
    }

}
