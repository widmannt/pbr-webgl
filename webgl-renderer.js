class RendererProgram extends WebGL2ProgramContainer {

    constructor(gl) {
        const program = WebGL2Utils.from(gl).createProgramFromSources(
            RendererProgram.vertShaderSource.trim(), RendererProgram.fragShaderSource.trim());
        super(gl, program);
    }

}

RendererProgram.vertShaderSource = `
    #version 300 es

    const vec2 vertices[4] = vec2[] (
        vec2(-1.0,  1.0),
        vec2(-1.0, -1.0),
        vec2( 1.0,  1.0),
        vec2( 1.0, -1.0)
    );

    void main() {
        gl_Position = vec4(vertices[gl_VertexID], 0.0, 1.0);
    }
`;

RendererProgram.fragShaderSource = `
    #version 300 es

    precision highp float;
    precision highp sampler2DArray;

    uniform highp sampler2DArray pixelData;
    uniform uint samples;

    //uniform float blackLevel;
    //uniform float whiteLevel;
    //uniform bvec4 indicators;
    const float blackLevel = 0.0;
    const float whiteLevel = 1.0;
    const bvec4 indicators = bvec4(false, false, false, false);

    out vec4 pixelColor;

    vec3 applyIndicators(vec3 color) {
        bvec3 greaterThanOne = greaterThan(color, vec3(1.0));
        bvec3 lessThanZero = lessThan(color, vec3(0.0));
        bool isClipped = false;

        if (indicators[0] && all(greaterThanOne)) {
            color = vec3(1.0, 0.0, 0.0);
            isClipped = true;
        }

        if (indicators[1] && all(lessThanZero)) {
            color = vec3(0.0, 0.0, 1.0);
            isClipped = true;
        }

        if (indicators[2] && !isClipped && (any(greaterThanOne) || any(lessThanZero))) {
            color = vec3(1.0, 1.0, 0.0);
            isClipped = true;
        }

        if (indicators[3] && !isClipped) {
            color = vec3((color.r + color.g + color.b)/3.0);
        }

        return color;
    }

    vec3 sRGB(vec3 color) {
        bvec3 cutoff = lessThan(color, vec3(0.0031308));
        vec3 higher = vec3(1.055)*pow(color, vec3(1.0/2.4)) - vec3(0.055);
        vec3 lower = color*vec3(12.92);
        return mix(higher, lower, cutoff);
    }

    void main() {
        vec4[] storedData = vec4[](
            texelFetch(pixelData, ivec3(gl_FragCoord.xy, 0), 0),
            texelFetch(pixelData, ivec3(gl_FragCoord.xy, 1), 0)
        );

        //vec3 color = storedData[0].rgb/float(samples);
        //pixelColor.rgb = applyIndicators((color - vec3(blackLevel))/vec3(whiteLevel - blackLevel));

        vec3 color = storedData[0].rgb/storedData[1][0];
        float luma = dot(color, vec3(0.2126, 0.7152, 0.0722)); // Source: https://en.wikipedia.org/wiki/Relative_luminance
	    pixelColor.rgb = (1.0/(1.0 + luma))*color; // Basic tone mapping.

        //pixelColor.rgb = pow(pixelColor.rgb, vec3(1.0/2.2)); // Slightly faster approximation.
        pixelColor.rgb = sRGB(pixelColor.rgb);
    }
`;
