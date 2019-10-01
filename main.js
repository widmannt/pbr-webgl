"use strict";

main();

function main() {

    let targetSamples;
    let resolution;

    const paramHandlers = {
        s: (values) => {
            if (values.length > 1) {
                console.warn("Multiple instances of 'samples' parameter, using the first value assigned.");
            }
            targetSamples = Number(values[0] || NaN);
            if (!Number.isInteger(targetSamples) || targetSamples < 0) {
                if (values[0]) console.warn("Invalid samples value: " + values[0]);
                targetSamples = Infinity;
            }
        },

        r: (values) => {
            if (values.length > 1) {
                console.warn("Multiple instances of 'resolution' parameter, using the first value assigned.");
            }
            if (values[0]) {
                const match = values[0].match(/^([1-9]\d*)x([1-9]\d*)$/);
                if (match) {
                    resolution = [+match[1], +match[2]];
                } else {
                    console.warn("Invalid resolution value: " + values[0]);
                }
            }
        }
    };

    const urlParams = new URLSearchParams(window.location.search);
    const urlParamKeys = new Set(urlParams.keys());

    Object.keys(paramHandlers).forEach(key => {
        paramHandlers[key](urlParams.getAll(key));
        urlParamKeys.delete(key);
    });
    urlParamKeys.forEach(key => console.warn("Unknown parameter: " + key));

    document.addEventListener("DOMContentLoaded", renderer);



    function renderer() {
        const {canvas, gl} = initGl();
        const pathTracer = new PathTracer(gl);

        (function updateResolution() {
            matchMedia(`screen and (resolution: ${window.devicePixelRatio}dppx)`)
                .addEventListener("change", updateResolution, {once: true});
            if (pathTracer.updateResolution(resolution)) {
                console.debug("buffer resolution: %dx%d", ...pathTracer.bufferResolution);
            }
        })();

        let samplesPerFrame = 1;
        let sync;

        let requestId = requestAnimationFrame(startRender);
        function startRender() {
            requestId = requestAnimationFrame(nextFrame);
            sync = renderFrame(samplesPerFrame);
        }

        function nextFrame() {
            requestId = requestAnimationFrame(nextFrame);

            if (gl.getSyncParameter(sync, gl.SYNC_STATUS) === gl.SIGNALED) {
                gl.deleteSync(sync);
                sync = renderFrame(samplesPerFrame);
                updateSampleCounter(pathTracer.nSample);
            } else {
                console.debug("can't keep up");
            }
        }

        function renderFrame(samples) {
            for (let i = 0; i < samples && pathTracer.nSample < targetSamples; ++i) {
                pathTracer.sample();
            }
            pathTracer.render();

            const sync = gl.fenceSync(gl.SYNC_GPU_COMMANDS_COMPLETE, 0);
            gl.flush();

            return sync;
        }


        document.getElementById("exportButton").addEventListener("click", () => {
            pathTracer.render();
            exportImage(canvas, "image.png");
        });

        let rendererRunning = true;
        const rendererButton = document.getElementById("rendererButton");
        rendererButton.addEventListener("click", () => toggleRendererStatus());

        function toggleRendererStatus(override) {
            const shouldRun = typeof override !== "undefined" ? override : !rendererRunning;

            if (shouldRun && !rendererRunning) {
                requestId = requestAnimationFrame(startRender);
                rendererButton.firstChild.data = "Pause";
            }
            if (!shouldRun && rendererRunning) {
                cancelAnimationFrame(requestId);
                if (sync !== null) {
                    gl.deleteSync(sync);
                }
                rendererButton.firstChild.data = "Resume";
            }

            rendererRunning = shouldRun;
        }

        initMouseListener(canvas, pathTracer);
    }



    function initGl() {
        const canvas = document.getElementById("glCanvas");
        const gl = canvas.getContext("webgl2", {alpha: false, depth: false, antialias: false});

        if (gl === null) {
            throw Error("Unable to initialize WebGL2. Your browser or machine may not support it.");
        }

        if (gl.getExtension("EXT_color_buffer_float") === null) {
            throw Error("Extension not supported: EXT_color_buffer_float");
        }

        const glUtils = WebGL2Utils.from(gl);

        return {canvas, gl, glUtils};
    }

    function initMouseListener(canvas, pathTracer) {
        const helper = new MouseDragHelper(document, canvas);

        const leftDragListener = {
            onStartDragging(location) {
                this.lastMouseLocation = normalizeCoordinates(pathTracer.cssResolution, location);
            },
            onDrag(location) {
                const newMouseLocation = normalizeCoordinates(pathTracer.cssResolution, location);
                pathTracer.camera.theta += newMouseLocation[1] - this.lastMouseLocation[1];
                pathTracer.camera.phi -= newMouseLocation[0] - this.lastMouseLocation[0];
                pathTracer.updateCamera();

                this.lastMouseLocation = newMouseLocation;
            },
            onStopDragging() {
                printCameraInfo(pathTracer)
            }
        };

        const middleDragListener = {
            buttons: 4,
            onStartDragging(location) {
                this.lastMouseLocation = normalizeCoordinates(pathTracer.cssResolution, location);
            },
            onDrag(location) {
                const newMouseLocation = normalizeCoordinates(pathTracer.cssResolution, location);
                pathTracer.camera.distance -= 10 * (newMouseLocation[1] - this.lastMouseLocation[1]);
                pathTracer.updateCamera();

                this.lastMouseLocation = newMouseLocation;
            },
            onStopDragging() {
                printCameraInfo(pathTracer)
            }
        };

        helper.addMouseDragListener(leftDragListener);
        helper.addMouseDragListener(middleDragListener);
    }

    const toolbar = document.getElementById("toolbar");
    function exportImage(canvas, fileName) {
        canvas.toBlob((blob) => {
            const a = document.createElement("a");
            a.appendChild(document.createTextNode(fileName));
            a.target = "_blank";
            a.rel = "noopener noreferrer";
            a.href = URL.createObjectURL(blob);
            a.download = fileName;
            toolbar.appendChild(a);
        });
    }

    function printCameraInfo(pathTracer) {
        console.debug("distance: %s, theta: %s, phi: %s",
            formatNumber(pathTracer.camera.distance, 2),
            formatNumber(pathTracer.camera.theta, 4),
            formatNumber(pathTracer.camera.phi, 4)
        );
    }

    const sampleCounter = document.getElementById("sampleCounter");
    function updateSampleCounter(sample) {
        sampleCounter.firstChild.data = sample;
    }

    function round(number, precision) {
        function shift(number, decimals) {
            const [value, exponent] = String(number).split("e");
            return Number(value + "e" + (Number(exponent || 0) + decimals));
        }
        return shift(Math.round(shift(number, +precision)), -precision);
    } // Source: https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/round$revision/1383484#A_better_solution

    function formatNumber(number, decimals) {
        return round(number, decimals).toFixed(decimals);
    }

    function normalizeCoordinates(resolution, location) {
        const rcpmin = 1 / Math.min(...resolution);
        return [(2 * location[0] - resolution[0]) * rcpmin, (resolution[1] - 2 * location[1]) * rcpmin];
    }

}
