window.onload = function() {
    const canvas = document.getElementById('glcanvas');
    const gl = canvas.getContext('webgl2');

    if (!gl) {
        alert('WebGL 2.0 not supported');
        return;
    }

    const fileInput = document.getElementById('file-input');
    const fileList = document.getElementById('file-list');

    const renderer = new Renderer(gl, canvas);

    fileInput.addEventListener('change', handleFiles, false);
    canvas.addEventListener('wheel', event => renderer.handleScroll(event), false);
    canvas.addEventListener('mousedown', event => renderer.handleMouseDown(event), false);
    canvas.addEventListener('mousemove', event => renderer.handleMouseMove(event), false);
    canvas.addEventListener('mouseup', event => renderer.handleMouseUp(event), false);
    canvas.addEventListener('contextmenu', e => e.preventDefault());

    function handleFiles(event) {
        const files = event.target.files;
        fileList.innerHTML = '';

        Array.from(files).forEach(file => {
            const li = document.createElement('li');
            li.textContent = file.name;
            li.addEventListener('click', () => loadOBJFile(file));
            fileList.appendChild(li);
        });
    }

    function loadOBJFile(file) {
        const reader = new FileReader();
        reader.onload = function(event) {
            const objText = event.target.result;
            const mesh = parseOBJ(objText);
            renderer.setMesh(mesh);
        };
        reader.readAsText(file);
    }

    function parseOBJ(text) {
        const objLoader = new OBJLoader();
        return objLoader.parse(text);
    }
};


class Renderer {
    constructor(gl, canvas) {
        this.gl = gl;
        this.canvas = canvas;
        this.scale = 1.0;
        this.translation = vec3.create();
        this.dragging = false;
        this.lastX = 0;
        this.lastY = 0;
        this.mesh = null;
        this.lineWidth = 0.01; // 控制线条宽度

        this.initGL();
    }

    initGL() {
        const gl = this.gl;

        const vertexShaderSource = `
            attribute vec4 a_position;
            attribute vec3 a_color;
            uniform mat4 u_matrix;
            varying vec3 v_color;

            void main() {
                gl_Position = u_matrix * a_position;
                v_color = a_color;
            }
        `;

        const fragmentShaderSource = `
            precision mediump float;
            varying vec3 v_color;
            void main() {
                gl_FragColor = vec4(v_color, 1.0); // 白色
            }
        `;
        
        const edgeVertexShaderSource = `
            attribute vec4 a_position;
            attribute vec3 a_color;
            uniform mat4 u_matrix;
            varying vec3 v_color;

            void main() {
                gl_Position = u_matrix * a_position;
                v_color = a_color;
            }
        `;

        const edgeFragmentShaderSource = `
            precision mediump float;
            varying vec3 v_color;

            void main() {
                gl_FragColor = vec4(v_color, 1.0); // 使用插值后的颜色
            }
        `;

        // 创建用于渲染三角形的着色器程序
        const vertexShader = this.createShader(gl, gl.VERTEX_SHADER, vertexShaderSource);
        const fragmentShader = this.createShader(gl, gl.FRAGMENT_SHADER, fragmentShaderSource);
        this.triangleProgram = this.createProgram(gl, vertexShader, fragmentShader);

        // 创建用于渲染边的着色器程序
        const edgeVertexShader = this.createShader(gl, gl.VERTEX_SHADER, edgeVertexShaderSource);
        const edgeFragmentShader = this.createShader(gl, gl.FRAGMENT_SHADER, edgeFragmentShaderSource);
        this.edgeProgram = this.createProgram(gl, edgeVertexShader, edgeFragmentShader);

        this.positionLocation = gl.getAttribLocation(this.triangleProgram, 'a_position');
        this.colorLocation = gl.getAttribLocation(this.triangleProgram, 'a_color');
        this.matrixLocation = gl.getUniformLocation(this.triangleProgram, 'u_matrix');

        this.edgePositionLocation = gl.getAttribLocation(this.edgeProgram, 'a_position');
        this.edgeColorLocation = gl.getAttribLocation(this.edgeProgram, 'a_color');
        this.edgeMatrixLocation = gl.getUniformLocation(this.edgeProgram, 'u_matrix');

        gl.enable(gl.DEPTH_TEST);
        gl.clearColor(1.0, 1.0, 1.0, 1.0);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    }

    createShader(gl, type, source) {
        const shader = gl.createShader(type);
        gl.shaderSource(shader, source);
        gl.compileShader(shader);
        const success = gl.getShaderParameter(shader, gl.COMPILE_STATUS);
        if (success) {
            return shader;
        }
        console.error(gl.getShaderInfoLog(shader));
        gl.deleteShader(shader);
    }

    createProgram(gl, vertexShader, fragmentShader) {
        const program = gl.createProgram();
        gl.attachShader(program, vertexShader);
        gl.attachShader(program, fragmentShader);
        gl.linkProgram(program);
        const success = gl.getProgramParameter(program, gl.LINK_STATUS);
        if (success) {
            return program;
        }
        console.error(gl.getProgramInfoLog(program));
        gl.deleteProgram(program);
    }

    setMesh(mesh) {
        const gl = this.gl;
        
        this.mesh = mesh;

        this.positionBuffer = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, this.positionBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(mesh.vertices), gl.STATIC_DRAW);

        this.colorBuffer = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, this.colorBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(mesh.colors), gl.STATIC_DRAW); 

        this.indexBuffer = gl.createBuffer();
        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.indexBuffer);
        gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(mesh.indices), gl.STATIC_DRAW);
        
        // Create edge vertices and colors for rendering edges with width
        const edgeData = this.createEdgeData(mesh.vertices, mesh.colors, mesh.indices, this.lineWidth);
        this.edgeVertices = edgeData.vertices;
        this.edgeColors = edgeData.colors;

        this.edgeVertexBuffer = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, this.edgeVertexBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(this.edgeVertices), gl.STATIC_DRAW);

        this.edgeColorBuffer = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, this.edgeColorBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(this.edgeColors), gl.STATIC_DRAW);

        // Calculate initial scale and translation based on the bounding box
        if (AUTO_ADJUST_WINDOW) {
            const width = mesh.maxX - mesh.minX;
            const height = mesh.maxY - mesh.minY;
            const maxDim = Math.max(width, height);

            // Scale to fit the max dimension within the range [-1, 1]
            this.scale = 2 / maxDim;
            // Translate to center the model in the view
            this.translation = vec3.fromValues(-(mesh.minX + mesh.maxX) / 2, -(mesh.minY + mesh.maxY) / 2, 0);
        } else {
            this.scale = 1.0;
            this.translation = vec3.create();
        }

        this.drawScene();
    }

    createEdgeData(vertices, colors, indices, lineWidth) {
        const edgeVertices = [];
        const edgeColors = [];
        const vertexMap = {};

        // Map vertex indices to their positions and colors
        for (let i = 0; i < vertices.length; i += 3) {
            vertexMap[i / 3] = {
                position: [vertices[i], vertices[i + 1], vertices[i + 2]],
                color: [colors[i], colors[i + 1], colors[i + 2]]
            };
        }

        // Create quads for each edge
        for (let i = 0; i < indices.length; i += 3) {
            const [a, b, c] = [indices[i], indices[i + 1], indices[i + 2]];

            this.addEdgeQuad(edgeVertices, edgeColors, vertexMap[a], vertexMap[b], lineWidth);
            this.addEdgeQuad(edgeVertices, edgeColors, vertexMap[b], vertexMap[c], lineWidth);
            this.addEdgeQuad(edgeVertices, edgeColors, vertexMap[c], vertexMap[a], lineWidth);
        }

        return { vertices: edgeVertices, colors: edgeColors };
    }

    addEdgeQuad(edgeVertices, edgeColors, v0, v1, lineWidth) {
        const dx = v1.position[0] - v0.position[0];
        const dy = v1.position[1] - v0.position[1];
        const len = Math.sqrt(dx * dx + dy * dy);
        const nx = (dy / len) * lineWidth;
        const ny = (-dx / len) * lineWidth;

        // Calculate color as half of the original vertex colors
        const color0 = v0.color.map(c => c * 0.85);
        const color1 = v1.color.map(c => c * 0.85);

        edgeVertices.push(
            v0.position[0] + nx, v0.position[1] + ny, v0.position[2],
            v1.position[0] + nx, v1.position[1] + ny, v1.position[2],
            v0.position[0] - nx, v0.position[1] - ny, v0.position[2],

            v1.position[0] + nx, v1.position[1] + ny, v1.position[2],
            v1.position[0] - nx, v1.position[1] - ny, v1.position[2],
            v0.position[0] - nx, v0.position[1] - ny, v0.position[2]
        );

        edgeColors.push(
            ...color0, ...color1, ...color0,
            ...color1, ...color1, ...color0
        );
    }

    handleScroll(event) {
        this.scale *= event.deltaY > 0 ? 1.1 : 0.9;
        this.drawScene();
    }

    handleMouseDown(event) {
        if (event.button === 2) {
            this.dragging = true;
            this.lastX = event.clientX;
            this.lastY = event.clientY;
        }
    }

    handleMouseMove(event) {
        if (this.dragging) {
            const base = 0.002 / this.scale;
            const dx = event.clientX - this.lastX;
            const dy = event.clientY - this.lastY;
            this.lastX = event.clientX;
            this.lastY = event.clientY;

            vec3.add(this.translation, this.translation, [dx * base, -dy * base, 0]);
            this.drawScene();
        }
    }

    handleMouseUp() {
        this.dragging = false;
    }

    drawScene() {
        if (!this.mesh) {
            return;
        }

        const gl = this.gl;

        const projectionMatrix = mat4.create();
        const modelMatrix = mat4.create();
        const matrix = mat4.create();

        // Create orthographic projection matrix
        mat4.ortho(projectionMatrix, -1, 1, -1, 1, -1, 1);

        // 先缩放再位移
        // Apply translation first to center the model
        mat4.scale(modelMatrix, modelMatrix, [this.scale, this.scale, this.scale]);
        mat4.translate(modelMatrix, modelMatrix, this.translation);
        // Then apply scaling to fit the model to the view

        // Combine the matrices
        mat4.multiply(matrix, projectionMatrix, modelMatrix);

        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

        // Use edge program and draw edge quads
        gl.useProgram(this.edgeProgram);
        gl.uniformMatrix4fv(this.edgeMatrixLocation, false, matrix);

        // Bind position buffer for edge quads
        gl.enableVertexAttribArray(this.edgePositionLocation);
        gl.bindBuffer(gl.ARRAY_BUFFER, this.edgeVertexBuffer);
        gl.vertexAttribPointer(this.edgePositionLocation, 3, gl.FLOAT, false, 0, 0);

        // Bind color buffer for edge quads
        gl.enableVertexAttribArray(this.edgeColorLocation);
        gl.bindBuffer(gl.ARRAY_BUFFER, this.edgeColorBuffer);
        gl.vertexAttribPointer(this.edgeColorLocation, 3, gl.FLOAT, false, 0, 0);

        // Draw edge quads
        gl.drawArrays(gl.TRIANGLES, 0, this.edgeVertices.length / 3);

        // Use triangle program and draw triangles
        gl.useProgram(this.triangleProgram);
        gl.uniformMatrix4fv(this.matrixLocation, false, matrix);

        gl.enableVertexAttribArray(this.positionLocation);
        gl.bindBuffer(gl.ARRAY_BUFFER, this.positionBuffer);
        gl.vertexAttribPointer(this.positionLocation, 3, gl.FLOAT, false, 0, 0);

        // Bind color buffer
        gl.enableVertexAttribArray(this.colorLocation);
        gl.bindBuffer(gl.ARRAY_BUFFER, this.colorBuffer);
        gl.vertexAttribPointer(this.colorLocation, 3, gl.FLOAT, false, 0, 0);

        // Draw triangles
        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.indexBuffer);
        gl.drawElements(gl.TRIANGLES, this.mesh.indices.length, gl.UNSIGNED_SHORT, 0);

    }
}
