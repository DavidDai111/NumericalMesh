const AUTO_ADJUST_WINDOW = true; // 设为 false 时，不会自动把模型移到中心并调整模型大小

class OBJLoader {
    constructor() {
        this.vertices = [];
        this.indices = [];
        this.minX = AUTO_ADJUST_WINDOW ? Infinity : -1;
        this.maxX = AUTO_ADJUST_WINDOW ? -Infinity : 1;
        this.minY = AUTO_ADJUST_WINDOW ? Infinity : -1;
        this.maxY = AUTO_ADJUST_WINDOW ? -Infinity : 1;
        this.minZ = Infinity;
        this.maxZ = -Infinity;
        this.colors = [];
    }

    parse(text) {
        const lines = text.split('\n');
        const color_ws = [];

        lines.forEach(line => {
            const parts = line.trim().split(/\s+/);
            if (parts[0] === 'v') {
                let vertex = parts.slice(1, 4).map(parseFloat);
                // Update min and max for x and y if auto adjustment is enabled
                if (AUTO_ADJUST_WINDOW) {
                    if (vertex[0] < this.minX) this.minX = vertex[0];
                    if (vertex[0] > this.maxX) this.maxX = vertex[0];
                    if (vertex[1] < this.minY) this.minY = vertex[1];
                    if (vertex[1] > this.maxY) this.maxY = vertex[1];
                }

                if (vertex[2] < this.minZ) this.minZ = vertex[2];
                if (vertex[2] > this.maxZ) this.maxZ = vertex[2];
                color_ws.push(vertex[2]);
                vertex[2] = 0.0;
                this.vertices.push(...vertex);
            } else if (parts[0] === 'f') {
                const indices = parts.slice(1).map(part => parseInt(part.split('/')[0]) - 1);
                this.indices.push(...indices)
            }
        });

        function mapValueToColor(f, minZ, maxZ) {
            // Normalize f to range [0, 1]
            if (minZ === maxZ) return [1.0, 1.0, 1.0];
            let t = (f - minZ) / (maxZ - minZ);
        
            // Define cold and warm colors as vec3```
            const colors = [
                { t: 0.0, color: { r: 0.23, g: 0.30, b: 0.75 } },  // 蓝色
                { t: 0.5, color: { r: 0.87, g: 0.86, b: 0.86 } },  // 灰色
                { t: 1.0, color: { r: 0.71, g: 0.02, b: 0.15 } }   // 红色
            ];
            function interpolate(start, end, t) {
                return start + (end - start) * t;
            }
            // 查找适当的插值区间
            for (let i = 0; i < colors.length - 1; i++) {
                if (t >= colors[i].t && t <= colors[i + 1].t) {
                    // 在两个颜色节点之间插值
                    const range = colors[i + 1].t - colors[i].t;
                    const scale = (t - colors[i].t) / range;
                    const r = interpolate(colors[i].color.r, colors[i + 1].color.r, scale);
                    const g = interpolate(colors[i].color.g, colors[i + 1].color.g, scale);
                    const b = interpolate(colors[i].color.b, colors[i + 1].color.b, scale);
                    return [r, g, b];
                }
            }
        
            return [r, g, b];
        }

        color_ws.forEach(w => this.colors.push(...mapValueToColor(w, this.minZ, this.maxZ)));

        return {
            vertices: this.vertices,
            indices: this.indices,
            minX: this.minX,
            maxX: this.maxX,
            minY: this.minY,
            maxY: this.maxY,
            colors: this.colors
        };
    }
}
