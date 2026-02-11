/* STRANGE ATTRACTOR MINER - MAIN.JS
   Features: 
   - ODE Flow Mining (Chaos Hunter)
   - RK4 Integration (High Accuracy)
   - "Tail Analysis" Filtering (Fixes Dots & Lines)
   - Mobile-First Performance Tuning
   - 16-Bit Float Integration
   - Auto-Matching JSON/PNG Export
   - Multi-Engine: Poly, Symmetric
   - FIXED: Sym Invariant Subspace Bug (Off-diagonal init)
   - FIXED: Density Color Shift (Normalized v_time by u_pointCount)
   - FIXED: Export Length (Decoupled Smoothness from Simulation Duration)
   - UPGRADE: True Physics Sub-stepping
   - FIXED: Initial Render Quality (Auto-refreshes to full quality)
   - FIXED: Export Opacity (Consistent with viewport)
   - FIXED: Slider Sync Glitch (Opacity now locked to GPU data)
   - UI UPDATE: Removed Invert Checkbox (Auto-switches based on Blend Mode)
   - FIXED: Additive Glow (Removed Tone Mapping to allow "blow out" effects)
   - NEW: POD (Print on Demand) Integration via Peecho
*/

// ==========================================
// 0. CONSTANTS & GLOBAL STATE
// ==========================================
const MAX_ALLOCATION = 50000000; 

// --- POD CONFIG ---
const POD_API_URL = "https://script.google.com/macros/s/AKfycbyy2EWZpZ_LofW4JVHesxmaRq5LgPGrlEfKC49U2mVdujPhA0rr2XKqwlrn-vbFR7rt/exec";

const RECAPTCHA_SITE_KEY = "6Le8WWgsAAAAADEo9EQKpu_ZMaGaN0PHcCw0y4cL";

// --- GLOBAL STATE VARIABLES ---
let pointCount = 0;
let gpuRenderedDensity = 1.0;
let gaussianTex = null;

let camZoom = 2.0; 
let camPanX = 0, camPanY = 0;
let currentQuat = [0, 0, 0, 1]; 
let currentCoeffs = null;
let currentGenType = 'poly';
let colorMode = 0;
let blendMode = 'NORMAL'; 
let colorSeed = [0.0, 0.0, 0.0];
let isInverted = true; // Default to true because default mode is Paint
let incBlack = true;
let incWhite = true;
let bgA = [0.9, 0.9, 0.9]; 
let bgB = [1.0, 1.0, 1.0]; 
let bgParams = [0.0, 0.5, 1.0]; 

// High Quality Defaults
let currentPhysicsSteps = 1000000; 
let currentOpacity = 0.02;          
let currentIntensity = 2.0;        
let currentGamma = 1.8;
let currentNoise = 0.20; 
let currentPointSize = 1.0;
let currentDensity = 1;
let currentJitter = 0.5;
let currentVariation = 0.0;
let isExporting = false; 

// --- VIEWPORT FBO REFS ---
let viewFbo = null;
let viewTex = null;

// --- INTERACTION STATE ---
let isDragging = false;
let isPanning = false;
let isRolling = false; 
let lastX = 0, lastY = 0;
let lastDist = 0;
let lastAngle = 0;
let renderScale = 1.0;     // Current resolution scale (1.0 = Full, 0.5 = Half)
let isInteracting = false; // Are we currently touching/dragging?
let pendingResize = false; // Flag to trigger FBO update


// ==========================================
// 1. THE WORKER (TAIL ANALYSIS OPTIMIZATION)
// ==========================================
// ==========================================
// 1. THE WORKER (INTERPOLATION ENGINE)
// ==========================================
const workerCode = `
    function mulberry32(a) {
        return function() {
          var t = a += 0x6D2B79F5;
          t = Math.imul(t ^ t >>> 15, t | 1);
          t ^= t + Math.imul(t ^ t >>> 7, t | 61);
          return ((t ^ t >>> 14) >>> 0) / 4294967296;
        }
    }

    self.onmessage = function(e) {
        if (e.data.type === 'mine') mine(e.data.genType);
        else if (e.data.type === 'mutate') mutate(e.data.coeffs, e.data.genType);
        else if (e.data.type === 'render') renderExisting(e.data);
    };

    // --- CATMULL-ROM SPLINE MATH ---
    // Interpolates between p1 and p2 using p0 and p3 as control points/tangents
    function catmullRom(p0, p1, p2, p3, t) {
        const v0 = (p2 - p0) * 0.5;
        const v1 = (p3 - p1) * 0.5;
        const t2 = t * t;
        const t3 = t * t2;
        return (2 * p1 - 2 * p2 + v0 + v1) * t3 + 
               (-3 * p1 + 3 * p2 - 2 * v0 - v1) * t2 + 
               v0 * t + p1;
    }

    function renderExisting(data) {
        const c = new Float32Array(data.coeffs);
        // physicsSteps now = "Segments". density = "Points per Segment".
        const result = generateTrace(data.physicsSteps, data.density, data.seedOffset||0, c, data.genType); 
        self.postMessage({
            type: 'found', 
            source: 'render', 
            coeffs: c, 
            genType: data.genType,
            density: data.density,
            buffer: result.buffer, 
            metaBuffer: result.metaBuffer
        }, [result.buffer, result.metaBuffer]);
    }

    function getRandomSprott() {
        return (Math.floor(Math.random() * 25) - 12) / 10.0;
    }

    function mine(genType) {
        let attempts = 0;
        let lastReport = Date.now();
        
        while(true) {
            attempts++;
            let coeffs;
            
            if (genType === 'sym') {
                coeffs = new Float32Array(10); 
                for(let i=0; i<10; i++) coeffs[i] = getRandomSprott(); 
            } else {
                coeffs = new Float32Array(30);
                for(let i=0; i<30; i++) coeffs[i] = (Math.random() * 2.4) - 1.2;
            }
            
            if (checkChaosTail(coeffs, genType)) {
                // Mining uses low density for speed (1 point per physics step)
                const result = generateTrace(50000, 1, 0, coeffs, genType);
                self.postMessage({
                    type: 'found', 
                    source: 'mine', 
                    coeffs: coeffs, 
                    genType: genType,
                    density: 1, 
                    buffer: result.buffer, 
                    metaBuffer: result.metaBuffer, 
                    attempts: attempts
                }, [result.buffer, result.metaBuffer]);
                break;
            }
            
            if (attempts % 1000 === 0) { 
                 if (Date.now() - lastReport > 500) {
                    self.postMessage({type: 'status', msg: 'Scanning (' + genType + ')... ' + attempts});
                    lastReport = Date.now();
                 }
            }
        }
    }

    function mutate(parentCoeffs, genType) {
        let attempts = 0;
        const parent = new Float32Array(parentCoeffs);
        
        while(true) {
            attempts++;
            let child = new Float32Array(parent);
            let idx = Math.floor(Math.random() * child.length);
            
            if (genType === 'sym') {
                child[idx] += (Math.random() < 0.5 ? 0.1 : -0.1); 
                child[idx] = Math.round(child[idx]*10)/10;
            } else {
                child[idx] += (Math.random() - 0.5) * 0.1;
            }

            if (checkChaosTail(child, genType)) {
                const result = generateTrace(50000, 1, 0, child, genType);
                self.postMessage({
                    type: 'found', 
                    source: 'mutate', 
                    coeffs: child, 
                    genType: genType,
                    density: 1,
                    buffer: result.buffer, 
                    metaBuffer: result.metaBuffer, 
                    attempts: attempts
                }, [result.buffer, result.metaBuffer]);
                break;
            }
            if (attempts % 1000 === 0) self.postMessage({type: 'status', msg: 'Mutating... ' + attempts});
        }
    }

    function checkChaosTail(c, genType) {
        // [UNCHANGED - Keeping strict physics check for mining]
        let x, y, z;
        if (genType === 'sym') { x = 0.1; y = 0.0; z = -0.1; } 
        else { x = 0.05; y = 0.05; z = 0.05; }

        let sx = x + 0.000001, sy = y, sz = z;
        let dt = (genType === 'sym') ? 0.01 : 0.02; 
        
        let a0,a1,a2,a3,a4,a5,a6,a7,a8,a9;
        let b0,b1,b2,b3,b4,b5,b6,b7,b8,b9;
        let c0,c1,c2,c3,c4,c5,c6,c7,c8,c9;
        if (genType === 'poly') {
           a0=c[0]; a1=c[1]; a2=c[2]; a3=c[3]; a4=c[4]; a5=c[5]; a6=c[6]; a7=c[7]; a8=c[8]; a9=c[9];
           b0=c[10]; b1=c[11]; b2=c[12]; b3=c[13]; b4=c[14]; b5=c[15]; b6=c[16]; b7=c[17]; b8=c[18]; b9=c[19];
           c0=c[20]; c1=c[21]; c2=c[22]; c3=c[23]; c4=c[24]; c5=c[25]; c6=c[26]; c7=c[27]; c8=c[28]; c9=c[29];
        } else {
           a0=c[0]; a1=c[1]; a2=c[2]; a3=c[3]; a4=c[4]; a5=c[5]; a6=c[6]; a7=c[7]; a8=c[8]; a9=c[9];
        }

        function calcD(px, py, pz, res) {
            if (genType === 'sym') {
                res.dx = a0 + a1*px + a2*py + a3*pz + a4*px*px + a5*py*py + a6*pz*pz + a7*px*py + a8*px*pz + a9*py*pz;
                res.dy = a0 + a1*py + a2*pz + a3*px + a4*py*py + a5*pz*pz + a6*px*px + a7*py*pz + a8*py*px + a9*pz*px;
                res.dz = a0 + a1*pz + a2*px + a3*py + a4*pz*pz + a5*px*px + a6*py*py + a7*pz*px + a8*pz*py + a9*px*py;
            } else {
                res.dx = a0 + a1*px + a2*py + a3*pz + a4*px*px + a5*py*py + a6*pz*pz + a7*px*py + a8*px*pz + a9*py*pz;
                res.dy = b0 + b1*px + b2*py + b3*pz + b4*px*px + b5*py*py + b6*pz*pz + b7*px*py + b8*px*pz + b9*py*pz;
                res.dz = c0 + c1*px + c2*py + c3*pz + c4*px*px + c5*py*py + c6*pz*pz + c7*px*py + c8*px*pz + c9*py*pz;
            }
        }
        
        let k1={dx:0,dy:0,dz:0}, k2={dx:0,dy:0,dz:0}, k3={dx:0,dy:0,dz:0}, k4={dx:0,dy:0,dz:0};
        
        for(let i=0; i<1000; i++) {
            calcD(x, y, z, k1);
            calcD(x + k1.dx*dt*0.5, y + k1.dy*dt*0.5, z + k1.dz*dt*0.5, k2);
            calcD(x + k2.dx*dt*0.5, y + k2.dy*dt*0.5, z + k2.dz*dt*0.5, k3);
            calcD(x + k3.dx*dt, y + k3.dy*dt, z + k3.dz*dt, k4);
            x += (k1.dx + 2*k2.dx + 2*k3.dx + k4.dx)*(dt/6);
            y += (k1.dy + 2*k2.dy + 2*k3.dy + k4.dy)*(dt/6);
            z += (k1.dz + 2*k2.dz + 2*k3.dz + k4.dz)*(dt/6);
            if (Math.abs(x) > 100 || isNaN(x)) return false; 
        }

        sx = x + 0.000001; sy = y; sz = z;
        let lyapunovSum = 0;
        let d0 = 0.000001;
        let minX=1e9, maxX=-1e9, minY=1e9, maxY=-1e9, minZ=1e9, maxZ=-1e9;
        let diagonalSum = 0;
        const visited = new Set();
        let voxRes = (genType === 'sym') ? 0.2 : 0.5;
        
        let steps = 2000;
        for(let i=0; i<steps; i++) {
            calcD(x, y, z, k1);
            calcD(x + k1.dx*dt*0.5, y + k1.dy*dt*0.5, z + k1.dz*dt*0.5, k2);
            calcD(x + k2.dx*dt*0.5, y + k2.dy*dt*0.5, z + k2.dz*dt*0.5, k3);
            calcD(x + k3.dx*dt, y + k3.dy*dt, z + k3.dz*dt, k4);
            x += (k1.dx + 2*k2.dx + 2*k3.dx + k4.dx)*(dt/6);
            y += (k1.dy + 2*k2.dy + 2*k3.dy + k4.dy)*(dt/6);
            z += (k1.dz + 2*k2.dz + 2*k3.dz + k4.dz)*(dt/6);

            calcD(sx, sy, sz, k1);
            calcD(sx + k1.dx*dt*0.5, sy + k1.dy*dt*0.5, sz + k1.dz*dt*0.5, k2);
            calcD(sx + k2.dx*dt*0.5, sy + k2.dy*dt*0.5, sz + k2.dz*dt*0.5, k3);
            calcD(sx + k3.dx*dt, sy + k3.dy*dt, sz + k3.dz*dt, k4);
            sx += (k1.dx + 2*k2.dx + 2*k3.dx + k4.dx)*(dt/6);
            sy += (k1.dy + 2*k2.dy + 2*k3.dy + k4.dy)*(dt/6);
            sz += (k1.dz + 2*k2.dz + 2*k3.dz + k4.dz)*(dt/6);

            if (Math.abs(x) > 100) return false;

            let dx = x - sx, dy = y - sy, dz = z - sz;
            let d = Math.sqrt(dx*dx + dy*dy + dz*dz);
            if (d < 1e-15) return false; 
            lyapunovSum += Math.log(d / d0);
            let s = d0 / d;
            sx = x - (dx * s); sy = y - (dy * s); sz = z - (dz * s);
            
            minX = Math.min(minX, x); maxX = Math.max(maxX, x);
            minY = Math.min(minY, y); maxY = Math.max(maxY, y);
            minZ = Math.min(minZ, z); maxZ = Math.max(maxZ, z);

            if (genType === 'sym') diagonalSum += (Math.abs(x-y) + Math.abs(y-z) + Math.abs(z-x));
            if (i % 10 === 0) visited.add(Math.floor(x/voxRes)+","+Math.floor(y/voxRes)+","+Math.floor(z/voxRes));
        }
        
        let lyapunov = lyapunovSum / steps;
        
        if (lyapunov < 0.001) return false;
        if (lyapunov > 2.0) return false;
        
        let wX = maxX - minX, wY = maxY - minY, wZ = maxZ - minZ;
        let minTotal = (genType === 'sym') ? 1.0 : 3.0; 
        if (wX < 0.1 || wY < 0.1 || wZ < 0.1) return false;
        if ((wX + wY + wZ) < minTotal) return false;
        
        let vThreshold = (genType === 'sym') ? 50 : 25; 
        if (visited.size < vThreshold) return false;
        
        if (genType === 'sym') {
            let avgDiag = diagonalSum / steps;
            if (avgDiag < 0.08) return false;
        }

        return true;
    }

    function generateTrace(physicsSteps, density, seedOffset, c, genType) {
        // --- NEW STRATEGY: PHYSICS SKELETON + SPLINE INTERPOLATION ---
        
        const totalPoints = physicsSteps * density; // Total OUTPUT points
        
        // Safety: Limit memory to avoid browser crash (approx 10M point limit here for safety)
        // If density is 50 and steps 2M, that's 100M points -> too big.
        // We trust the user, but let's be aware.
        
        let posData = new Float32Array(totalPoints * 3);
        let metaData = new Float32Array(totalPoints * 2); 
        const rand = mulberry32(seedOffset + 12345);

        let x, y, z;
        if (genType === 'sym') { x = 0.1; y = 0.0; z = -0.1; } 
        else { x = 0.05; y = 0.05; z = 0.05; }

        // [CHANGE] FIXED DT. We do NOT divide by density anymore.
        // This ensures the attractor shape stays constant regardless of smoothness.
        let dt = (genType === 'sym') ? 0.005 : 0.015;

        // Cache Coeffs
        let a0,a1,a2,a3,a4,a5,a6,a7,a8,a9;
        let b0,b1,b2,b3,b4,b5,b6,b7,b8,b9;
        let c0,c1,c2,c3,c4,c5,c6,c7,c8,c9;
        if (genType === 'poly') {
           a0=c[0]; a1=c[1]; a2=c[2]; a3=c[3]; a4=c[4]; a5=c[5]; a6=c[6]; a7=c[7]; a8=c[8]; a9=c[9];
           b0=c[10]; b1=c[11]; b2=c[12]; b3=c[13]; b4=c[14]; b5=c[15]; b6=c[16]; b7=c[17]; b8=c[18]; b9=c[19];
           c0=c[20]; c1=c[21]; c2=c[22]; c3=c[23]; c4=c[24]; c5=c[25]; c6=c[26]; c7=c[27]; c8=c[28]; c9=c[29];
        } else {
           a0=c[0]; a1=c[1]; a2=c[2]; a3=c[3]; a4=c[4]; a5=c[5]; a6=c[6]; a7=c[7]; a8=c[8]; a9=c[9];
        }

        function calcD(px, py, pz, res) {
            if (genType === 'sym') {
                res.dx = a0 + a1*px + a2*py + a3*pz + a4*px*px + a5*py*py + a6*pz*pz + a7*px*py + a8*px*pz + a9*py*pz;
                res.dy = a0 + a1*py + a2*pz + a3*px + a4*py*py + a5*pz*pz + a6*px*px + a7*py*pz + a8*py*px + a9*pz*px;
                res.dz = a0 + a1*pz + a2*px + a3*py + a4*pz*pz + a5*px*px + a6*py*py + a7*pz*px + a8*pz*py + a9*px*py;
            } else {
                res.dx = a0 + a1*px + a2*py + a3*pz + a4*px*px + a5*py*py + a6*pz*pz + a7*px*py + a8*px*pz + a9*py*pz;
                res.dy = b0 + b1*px + b2*py + b3*pz + b4*px*px + b5*py*py + b6*pz*pz + b7*px*py + b8*px*pz + b9*py*pz;
                res.dz = c0 + c1*px + c2*py + c3*pz + c4*px*px + c5*py*py + c6*pz*pz + c7*px*py + c8*px*pz + c9*py*pz;
            }
        }
        
        let k1={dx:0,dy:0,dz:0}, k2={dx:0,dy:0,dz:0}, k3={dx:0,dy:0,dz:0}, k4={dx:0,dy:0,dz:0};

        // Ring buffer for the last 4 physics points: [p0, p1, p2, p3]
        // We interpolate between p1 and p2.
        const history = []; 
        // Helper to run physics one step
        function stepPhysics() {
            calcD(x, y, z, k1);
            
            // Calc meta data for the CURRENT point (Start of segment)
            let speed = Math.sqrt(k1.dx*k1.dx + k1.dy*k1.dy + k1.dz*k1.dz);
            let logVel = Math.log(speed + 1.0);
            
            // Curvature (approximate from previous step k1 vs current k1??)
            // For speed, we'll just use a simple velocity magnitude check or 
            // re-use previous calc. 
            // Let's stick to velocity for now to keep interpolation fast.
            let curv = 0; // Simplified for this interpolation engine
            
            calcD(x + k1.dx*dt*0.5, y + k1.dy*dt*0.5, z + k1.dz*dt*0.5, k2);
            calcD(x + k2.dx*dt*0.5, y + k2.dy*dt*0.5, z + k2.dz*dt*0.5, k3);
            calcD(x + k3.dx*dt, y + k3.dy*dt, z + k3.dz*dt, k4);
            
            let nextX = x + (k1.dx + 2*k2.dx + 2*k3.dx + k4.dx)*(dt/6);
            let nextY = y + (k1.dy + 2*k2.dy + 2*k3.dy + k4.dy)*(dt/6);
            let nextZ = z + (k1.dz + 2*k2.dz + 2*k3.dz + k4.dz)*(dt/6);
            
            x = nextX; y = nextY; z = nextZ;
            
            return { x:x, y:y, z:z, vel: logVel, curv: curv };
        }

        // Warmup: Run and discard 2000 steps to get on track
        for(let i=0; i<2000; i++) stepPhysics();

        // Prime the history buffer with 4 points
        for(let i=0; i<4; i++) history.push(stepPhysics());

        let outIdx = 0;
        let maxVel = 0;

        // --- MAIN LOOP ---
        for(let i=0; i<physicsSteps; i++) {
            
            // Catmull-Rom Points: p0, p1, p2, p3
            let p0 = history[0];
            let p1 = history[1];
            let p2 = history[2];
            let p3 = history[3];
            
            // If any point is invalid (escaped to infinity), break
            if (Math.abs(p2.x) > 1000 || isNaN(p2.x)) break;

            // Generate 'density' points between p1 and p2
            for(let d=0; d<density; d++) {
                if (outIdx >= totalPoints) break;
                
                let t = d / density; // 0.0 to 0.99...
                
                // Interpolate Position
                let px = catmullRom(p0.x, p1.x, p2.x, p3.x, t);
                let py = catmullRom(p0.y, p1.y, p2.y, p3.y, t);
                let pz = catmullRom(p0.z, p1.z, p2.z, p3.z, t);
                
                // Linear Interpolate Metadata (Good enough for visuals)
                let pVel = p1.vel + (p2.vel - p1.vel) * t;
                let pCurv = p1.curv + (p2.curv - p1.curv) * t;
                
                if (pVel > maxVel) maxVel = pVel;

                posData[outIdx*3] = px;
                posData[outIdx*3+1] = py;
                posData[outIdx*3+2] = pz;
                
                metaData[outIdx*2] = pVel;
                metaData[outIdx*2+1] = pCurv;
                
                outIdx++;
            }

            // Slide window: Remove p0, Add new point at end
            history.shift();
            history.push(stepPhysics());
        }
        
        // Normalization & Centering
        if (maxVel > 0) {
            for(let i=0; i<totalPoints; i++) metaData[i*2] /= maxVel;
        }
        
        let sumX=0, sumY=0, sumZ=0;
        for(let i=0; i<totalPoints; i++) {
            sumX += posData[i*3]; sumY += posData[i*3+1]; sumZ += posData[i*3+2];
        }
        let avgX = sumX/totalPoints, avgY = sumY/totalPoints, avgZ = sumZ/totalPoints;
        let sumDistSq = 0;
        for(let i=0; i<totalPoints; i++) {
            let px = posData[i*3]-avgX; let py = posData[i*3+1]-avgY; let pz = posData[i*3+2]-avgZ;
            posData[i*3] = px; posData[i*3+1] = py; posData[i*3+2] = pz;
            sumDistSq += px*px + py*py + pz*pz;
        }
        let rms = Math.sqrt(sumDistSq / totalPoints);
        if (rms > 0) {
            let s = 0.5 / rms;
            for(let i=0; i<posData.length; i++) posData[i] *= s;
        }

        return { buffer: posData.buffer, metaBuffer: metaData.buffer };
    }
`;

// ==========================================
// 2. WEBGL SETUP
// ==========================================
const canvas = document.getElementById('glcanvas');
const gl = canvas.getContext('webgl2', { preserveDrawingBuffer: true });
if (!gl) {
    document.body.innerHTML = "<h1 style='color:white'>WebGL2 Not Supported</h1>";
    throw new Error("WebGL2 missing");
}

const floatExt = gl.getExtension("EXT_color_buffer_float");

gaussianTex = createGaussianTexture(); // <--- GENERATE SPRITE

const vsSource = `#version 300 es
in vec3 a_position;
in float a_vel;
in float a_curv;

uniform mat4 u_rotation; 
uniform vec2 u_pan;
uniform float u_zoom;
uniform float u_aspect;
uniform float u_pointSize;
uniform vec4 u_tileBounds; 
uniform vec2 u_resolution; 
uniform float u_jitter;     
uniform float u_pointCount; 
// uniform float u_density;  <-- REMOVED (Not needed without LOD)

uniform float u_focusDist; 
uniform float u_focusSpan; 
uniform float u_aperture;  

out float v_vel;
out float v_curv;
out float v_time;
out vec3 v_pos;
out float v_attenuation; 

float hash(float n) { return fract(sin(n) * 1e4); }

void main() {
    float safeCount = max(u_pointCount, 1.0);
    v_time = float(gl_VertexID) / safeCount; 
    
    v_vel = a_vel;
    v_curv = a_curv;
    v_pos = a_position;

    // --- 1. POSITION MATH ---
    vec4 rotated = u_rotation * vec4(a_position, 1.0);
    float dist = rotated.z; 

    // Depth Smearing (Kept - Helps with DoF smoothness)
    float zNoise = (hash(float(gl_VertexID) + 77.7) - 0.5) * 0.04;
    float blurredDist = dist + zNoise;

    rotated.xy *= u_zoom;
    rotated.xy += u_pan; 
    rotated.x /= u_aspect;
    
    vec2 normalizedPos = (rotated.xy * 0.5) + 0.5;
    vec2 tilePos = (normalizedPos - u_tileBounds.xy) / u_tileBounds.zw;
    
    gl_Position = vec4(tilePos * 2.0 - 1.0, 0.0, 1.0);
    
    // --- 2. DOF SIZING ---
    float distFromFocus = abs(blurredDist - u_focusDist);
    float effectiveDist = max(0.0, distFromFocus - u_focusSpan);

    // Quadratic Falloff (Kept - Looks natural)
    float rawBlur = (effectiveDist * effectiveDist) * (u_aperture * 0.1);
    
    float scaleFactor = u_resolution.y / 1080.0;
    float scaledBlur = rawBlur * scaleFactor;

    // CAP SIZE (Crash Protection)
    // We keep this. 400px is huge but safe.
    scaledBlur = clamp(scaledBlur, 0.0, 400.0 * scaleFactor);

    float baseSize = max(u_pointSize, 1.0);
    float targetSize = baseSize + scaledBlur;
    
    // Size Dithering (Kept - Hides rounding artifacts)
    float rnd = hash(float(gl_VertexID) + 555.555);
    float sizeJitter = (rnd - 0.5) * (targetSize * 0.15); 
    
    gl_PointSize = max(1.0, targetSize + sizeJitter); 
    
    // --- 3. ENERGY CONSERVATION ---
    // We removed the LOD multiplier. 
    // We stick to "Inverse Linear" (Bokeh Boost) because it looks best for art.
    // If this is too bright, change to: areaRatio * areaRatio;
    
    float areaRatio = baseSize / targetSize;
    v_attenuation = areaRatio; 
    
    // --- 4. POSITION JITTER ---
    float jx = (hash(float(gl_VertexID)) - 0.5) * u_jitter;
    float jy = (hash(float(gl_VertexID) + 123.45) - 0.5) * u_jitter;
    vec2 pixelSize = 2.0 / u_resolution; 
    gl_Position.xy += vec2(jx, jy) * pixelSize;
}`;

const fsSource = `#version 300 es
precision highp float;

in float v_vel;
in float v_curv;
in float v_time;
in vec3 v_pos;
in float v_attenuation; 

uniform sampler2D u_sprite; // <--- The optimization texture
uniform int u_colorMode;
uniform vec3 u_colorSeed; 
uniform float u_opacity;
uniform float u_intensity; 
uniform float u_noise;
uniform int u_invert;
uniform int u_inc_black; 
uniform int u_inc_white; 
uniform float u_variation; 

out vec4 outColor;

vec3 palette(float t, vec3 a, vec3 b, vec3 c, vec3 d) {
    return a + b * cos(6.28318 * (c * t + d));
}

void main() {
    // --- OPTIMIZED SPRITE LOGIC ---
    // Instead of calculating exp() math, we read the texture.
    // gl_PointCoord is 0.0 to 1.0 across the particle.
    // We read the Alpha (.a) from our Gaussian sprite.
    float shapeAlpha = texture(u_sprite, gl_PointCoord).a;
    
    // Discard invisible pixels to save fill-rate.
    // Since we use a texture, we don't need 'distSq' math anymore.
    if (shapeAlpha < 0.01) discard; 

    // --- COLOR LOGIC ---
    float val = 0.0;
    if (u_colorMode == 0) val = pow(v_vel, 0.5); 
    else if (u_colorMode == 1) val = v_time; 
    else if (u_colorMode == 2) val = pow(v_curv, 0.5); 
    else if (u_colorMode == 3) {
        vec3 p = v_pos + (u_colorSeed * 10.0);
        float pattern = sin(p.x * 5.0) + cos(p.y * 5.0) + sin(p.z * 5.0);
        val = 0.5 + (pattern * 0.25);
    } else {
        val = v_time * 2.0; 
    }
    val = clamp(val, 0.0, 1.0);
    
    vec3 phaseDrift = vec3(sin(v_time*3.0), cos(v_time*5.0), sin(v_time*7.0)) * (u_variation*0.2); 
    vec3 phaseShift = u_colorSeed + phaseDrift;

    vec3 rgb;
    if (u_colorMode == 0) rgb = palette(val, vec3(0.5), vec3(0.5), vec3(1.0), vec3(0.0, 0.33, 0.67) + phaseShift);
    else if (u_colorMode == 1) rgb = palette(val, vec3(0.5), vec3(0.5), vec3(1.0), vec3(0.0, 0.10, 0.20) + phaseShift);
    else if (u_colorMode == 2) rgb = palette(val, vec3(0.5), vec3(0.5), vec3(1.0), vec3(0.8, 0.8, 0.5) + phaseShift);
    else if (u_colorMode == 3) rgb = palette(val, vec3(0.5), vec3(0.5), vec3(1.0), vec3(0.0, 0.10, 0.20) + phaseShift);
    else rgb = palette(val, vec3(0.3), vec3(0.5), vec3(1.0), vec3(0.0, 0.33, 0.67) + phaseShift);

    if (u_inc_black == 0) {
        float luma = dot(rgb, vec3(0.299, 0.587, 0.114));
        rgb = mix(rgb, vec3(0.5), (1.0 - luma) * 0.3); 
    }
    if (u_inc_white == 0) rgb *= 0.85; 

    if (u_invert == 1) {
        rgb = 1.0 - rgb; 
        rgb = pow(rgb, vec3(1.5)); 
    }

    // --- APPLY OPACITY ---
    // Combine the Sprite Alpha with the User Opacity and Physics Attenuation
    float finalOpac = u_opacity * v_attenuation * shapeAlpha;
    
    outColor = vec4(rgb * u_intensity * finalOpac, finalOpac); 
}`;

const bgVsSource = `#version 300 es
void main() {
    float x = float((gl_VertexID & 1) << 2);
    float y = float((gl_VertexID & 2) << 1);
    gl_Position = vec4(x - 1.0, y - 1.0, 0, 1);
}`;

const bgFsSource = `#version 300 es
precision highp float;
uniform vec3 u_bg_a;
uniform vec3 u_bg_b;
uniform vec3 u_bg_params; 
uniform vec2 u_total_resolution; 
uniform vec2 u_tile_offset; 
uniform float u_noise; 
out vec4 outColor;

float hash(vec2 p) {
    vec3 p3  = fract(vec3(p.xyx) * .1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.x + p3.y) * p3.z);
}

void main() {
    vec2 globalCoord = gl_FragCoord.xy + u_tile_offset;
    vec2 uv = globalCoord / u_total_resolution;
    float aspect = u_total_resolution.x / u_total_resolution.y;
    vec2 aspectUV = uv; aspectUV.x *= aspect;
    vec2 center = vec2(0.5 * aspect, 0.5) + (vec2(cos(u_bg_params.x), sin(u_bg_params.x)) * u_bg_params.y);
    float dist = distance(aspectUV, center);
    float t = smoothstep(0.0, 1.5 * u_bg_params.z, dist);
    vec3 bg = mix(u_bg_a, u_bg_b, t);
    float noise = hash(gl_FragCoord.xy);
    bg += vec3(noise - 0.5) * (u_noise * 0.1); 
    outColor = vec4(bg, 1.0);
}`;

const compositeVsSource = `#version 300 es
void main() {
    float x = float((gl_VertexID & 1) << 2);
    float y = float((gl_VertexID & 2) << 1);
    gl_Position = vec4(x - 1.0, y - 1.0, 0, 1);
}`;

const compositeFsSource = `#version 300 es
precision highp float;
uniform sampler2D u_tex; 
uniform float u_passes; 
uniform float u_gamma; 
uniform vec3 u_bg_a;
uniform vec3 u_bg_b;
uniform vec3 u_bg_params;
uniform vec2 u_res;
uniform vec2 u_off;
uniform float u_noise; 
uniform int u_invert;
uniform int u_inc_black;
uniform int u_inc_white;
uniform int u_blend_mode; 
uniform float u_scale; // <--- NEW UNIFORM

uniform int u_show_guide;
uniform float u_print_aspect; 

out vec4 c;

float hash(vec2 p) {
    vec3 p3  = fract(vec3(p.xyx) * .1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.x + p3.y) * p3.z);
}

float triangularNoise(vec2 uv) {
    float r = hash(uv);
    float r2 = hash(uv + vec2(1.2, 3.4)); 
    return r - r2; 
}

void main() { 
    // 1. GLOBAL UV (For Backgrounds)
    // This stays relative to the TOTAL image (0..1)
    vec2 globalUV = (gl_FragCoord.xy + u_off) / u_res;
    float aspect = u_res.x / u_res.y;
    
    vec2 auv = globalUV; auv.x *= aspect;
    vec2 cnt = vec2(0.5*aspect,0.5) + (vec2(cos(u_bg_params.x), sin(u_bg_params.x))*u_bg_params.y);
    float t = smoothstep(0.0, 1.5*u_bg_params.z, distance(auv, cnt));
    vec3 bg = mix(u_bg_a, u_bg_b, t);
    
    bg += (hash(gl_FragCoord.xy)-0.5) * (u_noise*0.05);

    // 2. LOCAL FETCH (For Particles)
    // We go back to texelFetch (Integer lookup), but we apply the Render Scale.
    // If we are in "Fast Mode" (0.5 scale), screen pixel 100 -> texture pixel 50.
    // If we are in "Export Mode" (1.0 scale), tile pixel 100 -> texture pixel 100.
    
    ivec2 texCoord = ivec2(gl_FragCoord.xy * u_scale);
    vec4 pVal = texelFetch(u_tex, texCoord, 0) / u_passes;

    vec3 color = pVal.rgb;
    vec3 energy = color;

    // ... (Tone Mapping & Color Logic Unchanged) ...
    // Tone Mapping
    if (u_blend_mode == 0) { 
        float luma = length(color);
        float tone = 1.0 - exp(-luma * 1.0); 
        vec3 mapped = color * (tone / max(0.0001, luma));
        energy = mapped;
        vec3 gray = vec3(dot(energy, vec3(0.33)));
        energy = mix(gray, energy, 1.3);
    } 

    energy = max(vec3(0.0), energy); 
    energy = pow(energy, vec3(1.0/u_gamma));

    vec3 finalRGB;
    if (u_invert == 1) finalRGB = bg * (1.0 - energy); 
    else finalRGB = bg + energy; 
    
    if (u_show_guide == 1) {
        float viewAspect = u_res.x / u_res.y;
        float targetAspect = u_print_aspect;
        if (targetAspect < viewAspect) {
            float safeRatio = targetAspect / viewAspect;
            float margin = (1.0 - safeRatio) * 0.5;
            if (globalUV.x < margin || globalUV.x > (1.0 - margin)) {
                finalRGB *= 0.3; 
                finalRGB += vec3(0.1, 0.1, 0.1); 
            }
            float lineW = 1.0 / u_res.x; 
            if (abs(globalUV.x - margin) < lineW || abs(globalUV.x - (1.0-margin)) < lineW) {
                finalRGB = vec3(0.5, 0.5, 0.5); 
            }
        }
    }

    c = vec4(finalRGB, 1.0);
    
    float ditherStrength = max(u_noise * 0.2, 0.004); 
    c.rgb += triangularNoise(gl_FragCoord.xy) * ditherStrength;
}`;

function createShader(gl, type, source) {
    const shader = gl.createShader(type);
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    
    // ERROR CHECKING
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
        console.error("SHADER COMPILE ERROR:", gl.getShaderInfoLog(shader));
        gl.deleteShader(shader);
        return null;
    }
    return shader;
}

// --- OPTIMIZATION: PRE-CALCULATED GAUSSIAN SPRITE ---
// Instead of running expensive math for every pixel of every particle,
// we look up the brightness in this tiny texture.
function createGaussianTexture() {
    const size = 64; // 64x64 is plenty for a smooth gradient
    const canvas = document.createElement('canvas');
    canvas.width = size;
    canvas.height = size;
    const ctx = canvas.getContext('2d');

    // Draw a radial gradient (white center, black edge)
    const cx = size / 2;
    const cy = size / 2;
    const radius = size / 2;
    
    const grad = ctx.createRadialGradient(cx, cy, 0, cx, cy, radius);
    grad.addColorStop(0, 'rgba(255, 255, 255, 1.0)'); // Center
    grad.addColorStop(0.2, 'rgba(255, 255, 255, 0.8)'); // Core falloff
    grad.addColorStop(0.5, 'rgba(255, 255, 255, 0.2)'); // Mid falloff
    grad.addColorStop(1, 'rgba(0, 0, 0, 0.0)');       // Edge

    ctx.fillStyle = grad;
    ctx.fillRect(0, 0, size, size);

    // Upload to WebGL
    const tex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, canvas);
    
    // Linear filtering gives us "free" interpolation between the 64 pixels
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    
    return tex;
}

function createProgram(gl, vs, fs) {
    const p = gl.createProgram();
    gl.attachShader(p, createShader(gl, gl.VERTEX_SHADER, vs));
    gl.attachShader(p, createShader(gl, gl.FRAGMENT_SHADER, fs));
    gl.linkProgram(p);
    return p;
}

const particleProgram = createProgram(gl, vsSource, fsSource);
const bgProgram = createProgram(gl, bgVsSource, bgFsSource);
const compositeProgram = createProgram(gl, compositeVsSource, compositeFsSource);


// --- BUFFERS ---
const vao = gl.createVertexArray();
gl.bindVertexArray(vao);
const positionBuffer = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, positionBuffer);
const positionLoc = gl.getAttribLocation(particleProgram, "a_position");
gl.enableVertexAttribArray(positionLoc);
gl.vertexAttribPointer(positionLoc, 3, gl.FLOAT, false, 0, 0);

const metaBuffer = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, metaBuffer);
const velLoc = gl.getAttribLocation(particleProgram, "a_vel");
gl.enableVertexAttribArray(velLoc);
gl.vertexAttribPointer(velLoc, 1, gl.FLOAT, false, 8, 0);
const curvLoc = gl.getAttribLocation(particleProgram, "a_curv");
gl.enableVertexAttribArray(curvLoc);
gl.vertexAttribPointer(curvLoc, 1, gl.FLOAT, false, 8, 4);


// ==========================================
// 3. MATH & HELPERS
// ==========================================
function resizeViewportFBO() {
    if (viewFbo) gl.deleteFramebuffer(viewFbo);
    if (viewTex) gl.deleteTexture(viewTex);
    
    viewFbo = gl.createFramebuffer();
    viewTex = gl.createTexture();
    
    // 1. CALCULATE SCALED DIMENSIONS
    // If renderScale is 0.5, we create a texture half the size of the screen.
    const w = Math.floor(canvas.width * renderScale);
    const h = Math.floor(canvas.height * renderScale);

    gl.bindTexture(gl.TEXTURE_2D, viewTex);
    
    // 2. USE RGBA16F (Safe & Compatible)
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA16F, w, h, 0, gl.RGBA, gl.HALF_FLOAT, null);
    
    // 3. USE LINEAR FILTERING
    // Essential for Dynamic Res so the 0.5 scale image looks "soft" not "blocky"
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

    gl.bindFramebuffer(gl.FRAMEBUFFER, viewFbo);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, viewTex, 0);
    
    const status = gl.checkFramebufferStatus(gl.FRAMEBUFFER);
    if (status !== gl.FRAMEBUFFER_COMPLETE) {
        uiStatus.innerText = "FBO Error: " + status;
        uiStatus.style.color = "red";
    }
    
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
}

function hexToRgb(hex) { return [parseInt(hex.slice(1,3), 16)/255, parseInt(hex.slice(3,5), 16)/255, parseInt(hex.slice(5,7), 16)/255]; }
function rgbToHex(r,g,b) { return "#" + ((1<<24) + (Math.floor(r*255)<<16) + (Math.floor(g*255)<<8) + Math.floor(b*255)).toString(16).slice(1); }
function qIdentity() { return [0, 0, 0, 1]; }
function qMult(q1, q2) {
    let ax = q1[0], ay = q1[1], az = q1[2], aw = q1[3];
    let bx = q2[0], by = q2[1], bz = q2[2], bw = q2[3];
    return [ax * bw + aw * bx + ay * bz - az * by, ay * bw + aw * by + az * bx - ax * bz, az * bw + aw * bz + ax * by - ay * bx, aw * bw - ax * bx - ay * by - az * bz];
}
function qFromAxisAngle(axis, angle) {
    let half = angle * 0.5; let s = Math.sin(half);
    return [axis[0] * s, axis[1] * s, axis[2] * s, Math.cos(half)];
}
function qNormalize(q) {
    let len = Math.sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    if (len === 0) return [0,0,0,1];
    return [q[0]/len, q[1]/len, q[2]/len, q[3]/len];
}
function qToMatrix(q) {
    let x = q[0], y = q[1], z = q[2], w = q[3];
    let x2 = x + x, y2 = y + y, z2 = z + z;
    let xx = x * x2, xy = x * y2, xz = x * z2;
    let yy = y * y2, yz = y * z2, zz = z * z2;
    let wx = w * x2, wy = w * y2, wz = w * z2;
    return [1 - (yy + zz), xy + wz, xz - wy, 0, xy - wz, 1 - (xx + zz), yz + wx, 0, xz + wy, yz - wx, 1 - (xx + yy), 0, 0, 0, 0, 1];
}

// SERIALIZATION HELPERS
function generateID() {
    return Math.random().toString(36).substr(2, 6).toUpperCase();
}

function serializeState() {
    return {
        id: generateID(), // Generate new ID on save
        timestamp: Date.now(),
        genType: currentGenType, // SAVE GENERATOR TYPE
        camera: { zoom: camZoom, panX: camPanX, panY: camPanY, quat: Array.from(currentQuat) },
        colors: { 
            mode: colorMode, 
            seed: Array.from(colorSeed), 
            blend: blendMode, 
            invert: isInverted, 
            black: incBlack, 
            white: incWhite 
        },
        background: { a: Array.from(bgA), b: Array.from(bgB), params: Array.from(bgParams) },
        render: {
            steps: currentPhysicsSteps,
            density: currentDensity,
            opacity: currentOpacity,
            gamma: currentGamma,
            intensity: currentIntensity,
            noise: currentNoise,
            size: currentPointSize,
            jitter: currentJitter,
            variation: currentVariation
        },
        dof: {
            focus: document.getElementById('ui-focus').value,
            aperture: document.getElementById('ui-aperture').value,
            span: document.getElementById('ui-focus-span').value
        }
    };
}

function applyState(data) {
    if (!data) return;
    
    // Apply Settings
    if (data.settings) {
        const s = data.settings;
        
        if(s.genType) currentGenType = s.genType;
        if(s.camera) { camZoom = s.camera.zoom; camPanX = s.camera.panX; camPanY = s.camera.panY; currentQuat = s.camera.quat; }
        if(s.colors) { 
            colorMode = s.colors.mode; colorSeed = s.colors.seed; blendMode = s.colors.blend; 
            isInverted = (blendMode === 'NORMAL'); 
            incBlack = s.colors.black; incWhite = s.colors.white;
        }
        if(s.background) { bgA = s.background.a; bgB = s.background.b; bgParams = s.background.params; }
        if(s.render) {
            currentPhysicsSteps = s.render.steps; currentDensity = s.render.density;
            currentOpacity = s.render.opacity; currentGamma = s.render.gamma;
            currentIntensity = s.render.intensity; currentNoise = s.render.noise;
            currentPointSize = s.render.size; currentJitter = s.render.jitter;
            currentVariation = s.render.variation || 0;
        }

        // --- NEW: LOAD DOF SETTINGS ---
        // If the file has DoF data, use it. Otherwise reset to 0 (default).
        if (s.dof) {
            document.getElementById('ui-focus').value = s.dof.focus || 0;
            document.getElementById('ui-aperture').value = s.dof.aperture || 0;
            document.getElementById('ui-focus-span').value = s.dof.span || 0;
        } else {
            // Legacy file support: Reset DoF so it doesn't get stuck on previous settings
            document.getElementById('ui-focus').value = 0;
            document.getElementById('ui-aperture').value = 0;
            document.getElementById('ui-focus-span').value = 0;
        }

    } else {
        // Legacy file (reset camera)
        camZoom = 2.0; camPanX=0; camPanY=0; currentQuat=qIdentity();
    }
    
    // Update UI Elements to match state
    selectGenType.value = currentGenType;
    selectBlend.value = blendMode;
    selectColor.value = colorMode;
    checkIncBlack.checked = incBlack;
    checkIncWhite.checked = incWhite;
    inputBg1.value = rgbToHex(bgA[0],bgA[1],bgA[2]);
    inputBg2.value = rgbToHex(bgB[0],bgB[1],bgB[2]);
    
    // Force UI sliders to match the loaded values visually
    sliderLength.value = currentPhysicsSteps;
    sliderDensity.value = currentDensity;
    sliderOpacity.value = currentOpacity * 1000;
    sliderIntensity.value = currentIntensity * 10;
    sliderSize.value = currentPointSize * 10;
    sliderGamma.value = currentGamma * 10;
    sliderJitter.value = currentJitter * 10;
    sliderGamma.value = currentGamma * 10;
    sliderSmooth.value = currentNoise * 200;
    sliderVariation.value = currentVariation * 100;
}

// ==========================================
// 4. UI GENERATION
// ==========================================
const div = document.createElement('div');
div.id = 'colorControls';
div.style.position = 'absolute';
div.style.top = '20px';
div.style.right = '20px';
div.style.zIndex = '999'; 
div.style.background = 'rgba(0,0,0,0.85)';
div.style.width = '240px';
div.style.border = '1px solid #0f0';
div.style.fontFamily = 'monospace';
div.style.maxHeight = '90vh';
div.style.overflowY = 'auto';

function createSection(title, contentHTML) {
    const section = document.createElement('div');
    const header = document.createElement('div');
    header.innerHTML = `▶ ${title}`;
    header.style.background = '#222';
    header.style.color = '#fff';
    header.style.padding = '5px 10px';
    header.style.cursor = 'pointer';
    header.style.userSelect = 'none';
    header.style.borderBottom = '1px solid #444';
    const content = document.createElement('div');
    content.style.padding = '10px';
    content.style.display = 'none'; 
    content.innerHTML = contentHTML;
    header.onclick = () => {
        const isClosed = content.style.display === 'none';
        content.style.display = isClosed ? 'block' : 'none';
        header.innerHTML = `${isClosed ? '▼' : '▶'} ${title}`;
    };
    section.appendChild(header);
    section.appendChild(content);
    return section;
}

div.appendChild(createSection("GENERATION", `
    <div style="color:#0f0; margin-bottom:5px;">GENERATOR TYPE</div>
    <select id="ui-gen-type" style="width:100%; margin-bottom:10px;">
        <option value="poly">Polynomial (Cloud/Wire)</option>
        <option value="sym">Symmetric (CodeParade)</option>
    </select>
    <div style="display:flex; gap:5px; margin-bottom:10px;">
        <button id="ui-btn-mine" style="flex:1; cursor:pointer; background:#440000; color:#fff; border:1px solid #f00; padding:10px;">⛏️ MINE</button>
        <button id="ui-btn-mutate" style="flex:1; cursor:pointer; background:#440044; color:#fff; border:1px solid #f0f; padding:10px;">🧬 MUTATE</button>
    </div>
    <div style="display:flex; gap:5px; margin-bottom:10px;">
        <button id="ui-btn-save" style="flex:1; background:#222; color:#fff; border:1px solid #555; padding:5px;">💾 SAVE</button>
        <button id="ui-btn-load" style="flex:1; background:#222; color:#fff; border:1px solid #555; padding:5px;">📂 LOAD</button>
        <input type="file" id="ui-file-input" style="display:none" accept=".json">
    </div>
    <div style="color:#0f0; margin-bottom:5px;">LENGTH (Trail)</div>
    <input type="range" id="ui-length" min="10000" max="2000000" step="10000" value="1000000" style="width:100%;">
    <div style="color:#0f0; margin-bottom:5px;">SMOOTHNESS (Precision)</div>
    <input type="range" id="ui-density" min="1" max="50" value="1" style="width:100%;">
    <div style="color:#0f0; margin-bottom:5px;">VARIATION (Texture)</div>
    <input type="range" id="ui-variation" min="0" max="100" value="0" style="width:100%;">
`));

div.appendChild(createSection("APPEARANCE", `
    <div style="color:#0f0; margin-bottom:5px;">BLEND MODE</div>
    <select id="ui-blend-mode" style="width:100%; margin-bottom:10px;">
        <option value="NORMAL">Normal (Paint)</option>
        <option value="ADD">Additive (Glow)</option>
        <option value="SOFT">Soft (Dye)</option>
    </select>
    <div style="color:#0f0; margin-bottom:5px;">COLOR MODE</div>
    <select id="ui-color-mode" style="width:100%; margin-bottom:10px;">
        <option value="0">Velocity (Magma)</option>
        <option value="1">Time (Rainbow)</option>
        <option value="2">Curvature (Electric)</option>
        <option value="3">Spatial (Marble)</option>
        <option value="4">Ink (Spectral)</option>
    </select>
    <div style="display:flex; gap:10px; margin-bottom:10px; color:#fff;">
        <label><input type="checkbox" id="ui-inc-black" checked> Black</label>
        <label><input type="checkbox" id="ui-inc-white" checked> White</label>
    </div>
    <div style="display:flex; gap:5px; margin-bottom:10px;">
        <button id="ui-btn-color" style="flex:1; background:#222; color:#fff; border:1px solid #555; padding:5px;">🎲 COLOR</button>
        <button id="ui-btn-reset" style="flex:1; background:#222; color:#fff; border:1px solid #555; padding:5px;">👁 RESET</button>
    </div>
    
    <div style="color:#0f0; margin-top:10px; margin-bottom:5px;">BACKGROUND</div>
    <div style="display:flex; gap:5px; margin-bottom:5px;">
        <input type="color" id="ui-bg-color-1" value="#e0e0e0" style="flex:1;">
        <input type="color" id="ui-bg-color-2" value="#ffffff" style="flex:1;">
    </div>
    <button id="ui-btn-reroll-bg" style="width:100%; background:#222; color:#fff; border:1px solid #555; padding:5px;">🎲 REROLL BG</button>
`));

div.appendChild(createSection("ADJUSTMENT", `
    <div style="color:#0f0; margin-bottom:5px;">INTENSITY (Color)</div>
    <input type="range" id="ui-intensity" min="0" max="150" value="100" style="width:100%;">
    
    <div style="color:#0f0; margin-bottom:5px;">OPACITY (Alpha)</div>
    <input type="range" id="ui-opacity" min="1" max="300" value="20" style="width:100%;">
    
    <div style="color:#0f0; margin-bottom:5px;">SIZE (Thickness)</div>
    <input type="range" id="ui-size" min="1" max="50" value="10" style="width:100%;">
    <div style="color:#0f0; margin-bottom:5px;">GAMMA (Contrast)</div>
    <input type="range" id="ui-gamma" min="1" max="35" value="18" style="width:100%;">
    <div style="color:#0f0; margin-bottom:5px;">JITTER (Anti-Alias)</div>
    <input type="range" id="ui-jitter" min="0" max="20" value="5" style="width:100%;">
    <div style="color:#0f0; margin-bottom:5px;">DEBAND (Dither)</div>
    <input type="range" id="ui-deband" min="0" max="100" value="20" style="width:100%;">
    <div style="color:#0f0; margin-top:10px; border-top:1px solid #444; padding-top:5px;">DEPTH OF FIELD</div>
    <div style="display:flex; gap:5px; margin-bottom:5px;">
        <span style="color:#aaa; font-size:10px; width:40px; padding-top:4px;">FOCUS</span>
        <input type="range" id="ui-focus" min="-2000" max="2000" value="0" style="flex:1;">
    </div>
    <div style="display:flex; gap:5px; margin-bottom:5px;">
        <span style="color:#aaa; font-size:10px; width:40px; padding-top:4px;">SPAN</span>
        <input type="range" id="ui-focus-span" min="0" max="1000" value="0" style="flex:1;">
    </div>
    <div style="display:flex; gap:5px; margin-bottom:5px;">
        <span style="color:#aaa; font-size:10px; width:40px; padding-top:4px;">BLUR</span>
        <input type="range" id="ui-aperture" min="0" max="1000" value="0" style="flex:1;">
    </div>
`));

div.appendChild(createSection("EXPORT", `
    <label style="display:block; margin-bottom:10px; color:#aaa; font-size:11px; cursor:pointer;">
        <input type="checkbox" id="ui-show-guide"> Show Print Crop Guide
    </label>

    <div style="display:flex; gap:5px; margin-bottom:5px;">
        <input type="number" id="ui-print-w" value="24" style="width:40%" placeholder="W">
        <span style="color:#0f0; padding-top:4px;">in</span>
        <input type="number" id="ui-print-h" value="36" style="width:40%" placeholder="H">
    </div>
    <div style="display:flex; gap:5px; margin-bottom:10px;">
        <input type="number" id="ui-print-dpi" value="300" style="width:40%" placeholder="DPI">
        <span style="color:#0f0; padding-top:4px;">DPI</span>
        <input type="number" id="ui-print-passes" value="30" style="width:40%" placeholder="Passes">
    </div>
    <button id="ui-btn-snap" style="width:100%; background:#fff; color:#000; border:none; padding:10px; cursor:pointer; font-weight:bold; margin-bottom:5px;">📸 SNAPSHOT</button>
    <button id="ui-btn-order" style="width:100%; background:#0f0; color:#000; border:none; padding:10px; cursor:pointer; font-weight:bold;">🛒 ORDER PRINT</button>
    <div id="ui-export-status" style="color:#fff; font-size:10px; margin-top:5px; text-align:center;"></div>
`));

const footer = document.createElement('div');
footer.innerHTML = `<div id="ui-main-status" style="color:#ffff00; font-size:12px; margin-top:10px; text-align:center;">Initialized</div><div id="ui-error-log" style="color:red; font-size:10px; margin-top:5px;"></div>`;
div.appendChild(footer);

const oldUI = document.getElementById('colorControls');
if(oldUI) oldUI.remove();
document.body.appendChild(div);

const uiStatus = document.getElementById('ui-main-status');
const uiError = document.getElementById('ui-error-log');
const selectGenType = document.getElementById('ui-gen-type');
const btnMine = document.getElementById('ui-btn-mine');
const btnMutate = document.getElementById('ui-btn-mutate');
const btnSave = document.getElementById('ui-btn-save');
const btnLoad = document.getElementById('ui-btn-load');
const fileInput = document.getElementById('ui-file-input');
const inputBg1 = document.getElementById('ui-bg-color-1');
const inputBg2 = document.getElementById('ui-bg-color-2');
const btnRerollBg = document.getElementById('ui-btn-reroll-bg');
const selectBlend = document.getElementById('ui-blend-mode'); 
const selectColor = document.getElementById('ui-color-mode');
// REMOVED CHECK INVERT
const checkIncBlack = document.getElementById('ui-inc-black');
const checkIncWhite = document.getElementById('ui-inc-white');
const btnColor = document.getElementById('ui-btn-color');
const btnReset = document.getElementById('ui-btn-reset');
const sliderLength = document.getElementById('ui-length'); 
const sliderIntensity = document.getElementById('ui-intensity'); 
const sliderOpacity = document.getElementById('ui-opacity');
const sliderVariation = document.getElementById('ui-variation');
const sliderSize = document.getElementById('ui-size');
const sliderDensity = document.getElementById('ui-density');
const sliderJitter = document.getElementById('ui-jitter');
const sliderGamma = document.getElementById('ui-gamma');
const sliderSmooth = document.getElementById('ui-deband');
const btnSnap = document.getElementById('ui-btn-snap');
const btnOrder = document.getElementById('ui-btn-order');
const uiExport = document.getElementById('ui-export-status');
const inpW = document.getElementById('ui-print-w');
const inpH = document.getElementById('ui-print-h');
const inpDPI = document.getElementById('ui-print-dpi');
const inpPasses = document.getElementById('ui-print-passes');

// --- EVENT HANDLERS ---
window.onerror = function(msg, url, line) {
    uiError.innerText = `JS Error: ${msg} (Line ${line})`;
    return false;
};

selectGenType.onchange = (e) => { currentGenType = e.target.value; };
inputBg1.oninput = (e) => { bgA = hexToRgb(e.target.value); };
inputBg2.oninput = (e) => { bgB = hexToRgb(e.target.value); };
btnRerollBg.onclick = () => {
    if (blendMode === 'ADD') {
        // --- DARK SATURATED MODE (For Glow) ---
        // Pick a dominant channel for saturation
        const dominant = Math.floor(Math.random() * 3); 
        
        // Base intensity for the deep color (0.1 to 0.25 is "Dark" but visible)
        const intensity = 0.1 + Math.random() * 0.15; 
        
        const c = [0, 0, 0];
        c[dominant] = intensity; // Set dominant channel
        
        // Add a little bit to the other channels so it's not purely primary (more complex color)
        // but keep them low to maintain saturation.
        c[(dominant + 1) % 3] = Math.random() * (intensity * 0.3);
        c[(dominant + 2) % 3] = Math.random() * (intensity * 0.3);
        
        bgA = c; // Center color (Deep Hue)
        
        // Outer color: Fade to almost black (Vignette feel)
        bgB = [0.02, 0.02, 0.02]; 
        
    } else {
        // --- LIGHT MODE (Original Logic) ---
        const val = 0.8 + Math.random() * 0.2; 
        bgA = [val, val, val]; 
        const r = val - (Math.random() * 0.1);
        const g = val - (Math.random() * 0.1);
        const b = val - (Math.random() * 0.1);
        bgB = [r,g,b];
    }

    // Update Inputs & Params
    inputBg1.value = rgbToHex(bgA[0], bgA[1], bgA[2]);
    inputBg2.value = rgbToHex(bgB[0], bgB[1], bgB[2]);
    bgParams[0] = Math.random() * 6.28; 
    bgParams[1] = 0.2 + Math.random() * 0.6; 
    bgParams[2] = 0.8 + Math.random() * 0.5; 
};
selectBlend.onchange = (e) => { 
    blendMode = e.target.value; 
    // AUTO-SET INVERT: Paint=True, Glow/Dye=False
    isInverted = (blendMode === 'NORMAL');
}; 
selectColor.onchange = (e) => { colorMode = parseInt(e.target.value); };
checkIncBlack.onchange = (e) => { incBlack = e.target.checked; };
checkIncWhite.onchange = (e) => { incWhite = e.target.checked; };
btnColor.onclick = () => { colorSeed = [Math.random(), Math.random(), Math.random()]; };
btnReset.onclick = () => { currentQuat = qIdentity(); camPanX=0; camPanY=0; camZoom=2.0; };
sliderOpacity.oninput = (e) => { currentOpacity = parseInt(e.target.value) / 1000.0; };
sliderIntensity.oninput = (e) => { currentIntensity = parseInt(e.target.value) / 10.0; };
sliderSize.oninput = (e) => { currentPointSize = parseInt(e.target.value) / 10.0; };
sliderVariation.oninput = (e) => { currentVariation = parseInt(e.target.value) / 100.0; };
sliderLength.oninput = (e) => {
    currentPhysicsSteps = parseInt(e.target.value);
    if(currentCoeffs) worker.postMessage({ type: 'render', coeffs: currentCoeffs, physicsSteps: currentPhysicsSteps, density: currentDensity, genType: currentGenType });
};
sliderDensity.oninput = (e) => { 
    currentDensity = parseInt(e.target.value);
    if(currentCoeffs) worker.postMessage({ type: 'render', coeffs: currentCoeffs, physicsSteps: currentPhysicsSteps, density: currentDensity, genType: currentGenType });
};
sliderJitter.oninput = (e) => { currentJitter = parseInt(e.target.value) / 10.0; };
sliderGamma.oninput = (e) => { currentGamma = parseInt(e.target.value) / 10.0; };
sliderSmooth.oninput = (e) => { currentNoise = parseInt(e.target.value) / 200.0; };
btnSnap.onclick = () => { startTiledExport('download'); };
btnOrder.onclick = () => { startTiledExport('pod'); };
btnMine.onclick = () => { uiStatus.innerText = "Scanning..."; uiStatus.style.color = "#ffff00"; worker.postMessage({type: 'mine', genType: currentGenType}); };
btnMutate.onclick = () => {
    if (!currentCoeffs) { uiStatus.innerText = "Mine first!"; return; }
    uiStatus.innerText = "Mutating..."; uiStatus.style.color = "#ff00ff";
    worker.postMessage({ type: 'mutate', coeffs: currentCoeffs, genType: currentGenType });
};

// UPDATED SAVE: Now includes settings & ID
btnSave.onclick = () => { 
    if (!currentCoeffs) return; 
    const meta = serializeState();
    const data = { coeffs: Array.from(currentCoeffs), settings: meta };
    const blob = new Blob([JSON.stringify(data, null, 2)], {type: "application/json"});
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a'); 
    a.href = url; 
    a.download = `attractor_${meta.id}.json`; 
    document.body.appendChild(a); a.click(); document.body.removeChild(a); 
};

// UPDATED LOAD: Now applies settings
btnLoad.onclick = () => fileInput.click();
fileInput.onchange = (e) => { 
    const r = new FileReader(); 
    r.onload = (ev) => { 
        try { 
            const d = JSON.parse(ev.target.result); 
            currentCoeffs = Object.values(d.coeffs || d).map(Number); 
            applyState(d);
            worker.postMessage({ type: 'render', coeffs: currentCoeffs, physicsSteps: currentPhysicsSteps, density: currentDensity, genType: currentGenType }); 
            uiStatus.innerText = "Loaded " + (d.settings ? d.settings.id : "Legacy"); 
        } catch(e){console.error(e);} 
    }; 
    if(e.target.files[0]) r.readAsText(e.target.files[0]); 
};


canvas.oncontextmenu = (e) => e.preventDefault();

function applyTrackballRotation(dx, dy) {
    const dist = Math.hypot(dx, dy);
    if (dist < 0.01) return;
    const axis = [-dy, -dx, 0]; 
    const angle = dist * 0.005;
    const qDelta = qFromAxisAngle(axis, angle);
    currentQuat = qMult(qDelta, currentQuat);
    currentQuat = qNormalize(currentQuat);
}

function applyRoll(angleDelta) {
    const qDelta = qFromAxisAngle([0,0,-1], angleDelta);
    currentQuat = qMult(qDelta, currentQuat); 
    currentQuat = qNormalize(currentQuat);
}

canvas.addEventListener('touchstart', (e) => {
    if(isExporting) return;
    isInteracting = true;
    e.preventDefault();
    if (e.touches.length === 1) {
        lastX = e.touches[0].clientX; lastY = e.touches[0].clientY;
        const cx = canvas.width / 2;
        const cy = canvas.height / 2;
        const dist = Math.hypot(lastX - cx, lastY - cy);
        const maxDist = Math.min(canvas.width, canvas.height) / 2;
        if (dist > maxDist * 0.8) {
            isRolling = true;
            isDragging = false;
            lastAngle = Math.atan2(lastY - cy, lastX - cx);
        } else {
            isRolling = false;
            isDragging = true;
        }
        isPanning = false;
    } else if (e.touches.length === 2) {
        isDragging = false; isPanning = true; isRolling = false;
        lastDist = Math.hypot(e.touches[0].clientX - e.touches[1].clientX, e.touches[0].clientY - e.touches[1].clientY);
        lastX = (e.touches[0].clientX + e.touches[1].clientX) / 2;
        lastY = (e.touches[0].clientY + e.touches[1].clientY) / 2;
    }
}, {passive:false});

canvas.addEventListener('touchmove', (e) => {
    if(isExporting) return;
    e.preventDefault();
    if (e.touches.length === 1) {
        const cx = canvas.width / 2;
        const cy = canvas.height / 2;
        const x = e.touches[0].clientX;
        const y = e.touches[0].clientY;
        if (isRolling) {
            const angle = Math.atan2(y - cy, x - cx);
            const delta = angle - lastAngle;
            applyRoll(delta);
            lastAngle = angle;
            lastX = x; lastY = y;
        } else if (isDragging) {
            const dx = x - lastX;
            const dy = y - lastY;
            applyTrackballRotation(dx, dy);
            lastX = x; lastY = y;
        }
    } else if (e.touches.length === 2 && isPanning) {
        const newDist = Math.hypot(e.touches[0].clientX - e.touches[1].clientX, e.touches[0].clientY - e.touches[1].clientY);
        const cx = (e.touches[0].clientX + e.touches[1].clientX) / 2;
        const cy = (e.touches[0].clientY + e.touches[1].clientY) / 2;
        const zoomDelta = newDist / lastDist;
        camZoom *= zoomDelta;
        const pdx = (cx - lastX) / canvas.height * 2;
        const pdy = (cy - lastY) / canvas.height * 2;
        camPanX += pdx; camPanY -= pdy; 
        lastDist = newDist; lastX = cx; lastY = cy;
    }
}, {passive:false});

canvas.addEventListener('touchend', (e) => {
  isInteracting = false;
  isDragging=false;
  isPanning=false;
  isRolling=false;
});

canvas.onmousedown = (e) => {
    isInteracting = true;
    if(isExporting) return;
    lastX = e.clientX; lastY = e.clientY;
    if (e.button === 2) {
        isPanning = true; isDragging = false; isRolling = false;
    } else {
        isPanning = false;
        const cx = canvas.width / 2;
        const cy = canvas.height / 2;
        const dist = Math.hypot(lastX - cx, lastY - cy);
        const maxDist = Math.min(canvas.width, canvas.height) / 2;
        if (dist > maxDist * 0.8) {
            isRolling = true;
            isDragging = false;
            lastAngle = Math.atan2(lastY - cy, lastX - cx);
        } else {
            isRolling = false;
            isDragging = true;
        }
    }
};

window.onmouseup = () => {
  isInteracting = false;
  isDragging = false; isPanning = false; isRolling = false;
};

window.onmousemove = (e) => {
    if (isExporting) return;
    const x = e.clientX; const y = e.clientY;
    if (isRolling) {
        const cx = canvas.width / 2;
        const cy = canvas.height / 2;
        const angle = Math.atan2(y - cy, x - cx);
        const delta = angle - lastAngle;
        applyRoll(delta);
        lastAngle = angle;
    } else if (isDragging) {
        const dx = x - lastX;
        const dy = y - lastY;
        applyTrackballRotation(dx, dy);
    } else if (isPanning) {
        const dx = x - lastX;
        const dy = y - lastY;
        camPanX += (dx / canvas.height) * 2;
        camPanY -= (dy / canvas.height) * 2;
    }
    lastX = x; lastY = y;
};

canvas.onwheel = (e) => {
    isInteracting = true;

    // Debounce: Reset to false 200ms after scrolling stops
    clearTimeout(window.scrollTimeout);
    window.scrollTimeout = setTimeout(() => { isInteracting = false; }, 200);
  
    if(isExporting) return;
    e.preventDefault();

    // Use multiplicative zoom for fine-grained control at any scale.
    // e.deltaY is usually +/- 100 for a mouse wheel notch.
    // 1.0005 ^ 100 ≈ 1.05 (5% change per notch), which is smooth but responsive.
    const factor = Math.pow(1.00025, -e.deltaY);
    
    camZoom *= factor;
    
    // Cap the minimum zoom to prevent inverting or getting stuck
    camZoom = Math.max(0.01, camZoom); 
};

// ==========================================
// 5. EXPORT & WORKER INITIALIZATION (FINAL)
// ==========================================
const blob = new Blob([workerCode], { type: 'application/javascript' });
worker = new Worker(URL.createObjectURL(blob)); // Assigned to top-level var

worker.onerror = function(e) {
    uiStatus.innerText = "Worker Error: " + e.message;
    uiStatus.style.color = "red";
    uiError.innerText = "Worker Crash details: " + JSON.stringify(e);
};

worker.onmessage = (e) => {
    if (e.data.type === 'status') {
        uiStatus.innerText = e.data.msg;
    } 
    else if (e.data.type === 'found') {
        const data = new Float32Array(e.data.buffer);
        gl.bindBuffer(gl.ARRAY_BUFFER, positionBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, data, gl.STATIC_DRAW);
        const meta = new Float32Array(e.data.metaBuffer);
        gl.bindBuffer(gl.ARRAY_BUFFER, metaBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, meta, gl.STATIC_DRAW);
        pointCount = data.length / 3;
        
        // [FIX] Sync GPU density to prevent slider glitch
        // If density wasn't passed (legacy), default to 1.0
        gpuRenderedDensity = e.data.density || 1.0; 

        if (isExporting && exportResolve) {
            exportResolve();
            return;
        }

        uiStatus.innerText = `FOUND! (${e.data.attempts} attempts)`;
        uiStatus.style.color = "#00ff00";
        currentCoeffs = new Float32Array(e.data.coeffs);
        
        if (e.data.source === 'mine') {
            colorSeed = [Math.random(), Math.random(), Math.random()];
            camPanX = 0; camPanY = 0; camZoom = 2.0; 
            currentQuat = qIdentity();
            
            // [FIX] Auto-refresh to full quality immediately
            uiStatus.innerText = "Refining...";
            worker.postMessage({ 
                type: 'render', 
                coeffs: currentCoeffs, 
                physicsSteps: currentPhysicsSteps, 
                density: currentDensity, 
                genType: currentGenType 
            });
        }
    }
};

let exportResolve = null;

// --- UPDATED POD FUNCTION WITH RECAPTCHA ---
async function startPrintCheckout(blob) {
  const uiExport = document.getElementById('ui-export-status');
  const container = document.getElementById('colorControls');
  
  if (!blob) {
      uiExport.innerText = "Error: No image data generated.";
      return;
  }

  // Check if reCAPTCHA is loaded
  if (typeof grecaptcha === 'undefined') {
      uiExport.innerText = "Error: reCAPTCHA not loaded. Refresh page.";
      return;
  }

  uiExport.innerText = "🤖 Verifying you are human...";
  const allButtons = document.querySelectorAll('button');
  allButtons.forEach(b => b.disabled = true);

  try {
    // 1. GET FRESH RECAPTCHA TOKEN
    // We do this here (not on click) so the token is fresh even if rendering took 5 mins.
    const token = await new Promise((resolve) => {
        grecaptcha.ready(() => {
            grecaptcha.execute(RECAPTCHA_SITE_KEY, {action: 'print_order'}).then(resolve);
        });
    });

    uiExport.innerText = "⏳ Requesting Cloud Storage...";

    // 2. Upload Sequence (Pass token in URL)
    // We append the token as a query parameter
    const authResp = await fetch(POD_API_URL + "?recaptcha_token=" + token);
    const authData = await authResp.json();
    
    if (authData.error) throw new Error(authData.error);

    uiExport.innerText = "☁️ Uploading High-Res Image...";
    const uploadResp = await fetch(authData.uploadUrl, {
      method: "PUT",
      body: blob,
      headers: { "Content-Type": "image/png" }
    });

    if (!uploadResp.ok) throw new Error("Upload failed: " + uploadResp.statusText);

    // We MUST keep the ?Signature=... part
    const signedUrl = authData.publicUrl;
    
    uiExport.innerText = "✅ Ready.";

    // 3. Prepare Container
    const actionContainerId = 'pod-action-container';
    let actionContainer = document.getElementById(actionContainerId);
    if (actionContainer) actionContainer.remove(); 

    actionContainer = document.createElement('div');
    actionContainer.id = actionContainerId;
    actionContainer.style.marginTop = "10px";
    actionContainer.style.borderTop = "1px solid #555";
    actionContainer.style.paddingTop = "10px";

    const inchesW = parseFloat(inpW.value);
    const inchesH = parseFloat(inpH.value);
    const dpi = parseInt(inpDPI.value);
    
    const totalW = Math.floor(inchesW * dpi);
    const totalH = Math.floor(inchesH * dpi);
    
    const masterCanvas = document.createElement('canvas');
    masterCanvas.width = totalW;
    masterCanvas.height = totalH;

    // --- 4. REPLICATE THE WORKING POC ---
    
    // Create the <a> wrapper
    const peechoLink = document.createElement('a');
    peechoLink.href = "https://www.peecho.com"; 
    peechoLink.className = "peecho-print-button";    
    
    // Set Data Attributes
    peechoLink.setAttribute('data-src', signedUrl);
    peechoLink.setAttribute('data-thumbnail', signedUrl);
    peechoLink.setAttribute('data-filetype', 'image');
    peechoLink.setAttribute('data-width', masterCanvas.width);
    peechoLink.setAttribute('data-height', masterCanvas.height);
    peechoLink.setAttribute('data-pages', '1');
    peechoLink.setAttribute('data-reference', 'Strange Attractor #' + generateID());

    // --- BUTTON STYLE ---
    const innerBtn = document.createElement('button');
    innerBtn.innerText = "Review order options ↗️"; 
    
    innerBtn.style.backgroundColor = "#222";
    innerBtn.style.color = "#fff"; 
    innerBtn.style.padding = "10px";
    innerBtn.style.border = "1px solid #0f0"; 
    innerBtn.style.fontFamily = "monospace"; 
    innerBtn.style.fontSize = "12px";
    innerBtn.style.cursor = "pointer";
    innerBtn.style.width = "100%"; 
    innerBtn.style.fontWeight = "bold";

    peechoLink.appendChild(innerBtn);
    actionContainer.appendChild(peechoLink);

    // Insert into DOM
    const footer = document.getElementById('ui-main-status')?.parentNode;
    if(footer) container.insertBefore(actionContainer, footer);
    else container.appendChild(actionContainer);

    // --- 5. INJECT SCRIPT ---
    const scriptId = 'peecho-sdk-script';
    if (!document.getElementById(scriptId)) {
        const script = document.createElement('script');
        script.id = scriptId;
        script.src = "//d3aln0nj58oevo.cloudfront.net/button/script/177021676045966460.js";
        document.body.appendChild(script);
    } else {
        if (window.Peecho && window.Peecho.refresh) {
            console.log("Refreshing Peecho buttons...");
            window.Peecho.refresh();
        }
    }

    // Re-enable buttons
    allButtons.forEach(b => b.disabled = false);

  } catch (err) {
    console.error(err);
    uiExport.innerText = "Error: " + err.message;
    allButtons.forEach(b => b.disabled = false);
  }
}

// ADDED: overridePanX argument
// ADDED: overridePanX argument
function renderTileParticles(totalW, totalH, tileBounds, opac, forcedAspect, jitter, overridePanX) {
    // --- A. BIND OPTIMIZED GAUSSIAN SPRITE ---
    // (Must bind this for exports too!)
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, gaussianTex); 
    gl.uniform1i(gl.getUniformLocation(particleProgram, "u_sprite"), 0);

    const rotMatrix = qToMatrix(currentQuat);
    gl.uniformMatrix4fv(gl.getUniformLocation(particleProgram, "u_rotation"), false, new Float32Array(rotMatrix));
    
    const usePanX = (overridePanX !== undefined) ? overridePanX : camPanX;
    gl.uniform2f(gl.getUniformLocation(particleProgram, "u_pan"), usePanX, camPanY);
    
    gl.uniform1f(gl.getUniformLocation(particleProgram, "u_zoom"), camZoom);
    gl.uniform1f(gl.getUniformLocation(particleProgram, "u_aspect"), forcedAspect);
    gl.uniform4f(gl.getUniformLocation(particleProgram, "u_tileBounds"), tileBounds[0], tileBounds[1], tileBounds[2], tileBounds[3]);

    gl.uniform1i(gl.getUniformLocation(particleProgram, "u_colorMode"), colorMode);
    gl.uniform3f(gl.getUniformLocation(particleProgram, "u_colorSeed"), colorSeed[0], colorSeed[1], colorSeed[2]);
    gl.uniform1i(gl.getUniformLocation(particleProgram, "u_invert"), 0); 
    gl.uniform1i(gl.getUniformLocation(particleProgram, "u_inc_black"), 1); 
    gl.uniform1i(gl.getUniformLocation(particleProgram, "u_inc_white"), 1); 
    gl.uniform1f(gl.getUniformLocation(particleProgram, "u_variation"), currentVariation); 

    // --- DoF Settings ---
    const focusVal = (parseInt(document.getElementById('ui-focus').value) / 1000.0); 
    const focusSpanVal = (parseInt(document.getElementById('ui-focus-span').value) / 1000.0);
    const apertureVal = parseInt(document.getElementById('ui-aperture').value); // No divider (0-1000)

    gl.uniform1f(gl.getUniformLocation(particleProgram, "u_focusDist"), focusVal);
    gl.uniform1f(gl.getUniformLocation(particleProgram, "u_focusSpan"), focusSpanVal);
    gl.uniform1f(gl.getUniformLocation(particleProgram, "u_aperture"), apertureVal);

    // --- Opacity Math (Export) ---
    // For export, we use currentDensity because the worker just re-simulated it.
    let targetOpacity = opac / currentDensity;
    let safeOpacity = Math.max(targetOpacity, 0.000001); 
    
    gl.uniform1f(gl.getUniformLocation(particleProgram, "u_opacity"), safeOpacity);
    gl.uniform1f(gl.getUniformLocation(particleProgram, "u_intensity"), currentIntensity);
    gl.uniform1f(gl.getUniformLocation(particleProgram, "u_noise"), currentNoise);
    
    gl.uniform2f(gl.getUniformLocation(particleProgram, "u_resolution"), totalW, totalH);
    gl.uniform1f(gl.getUniformLocation(particleProgram, "u_jitter"), (jitter !== undefined) ? jitter : currentJitter);
    
    gl.uniform1f(gl.getUniformLocation(particleProgram, "u_pointCount"), pointCount);

    const scalar = totalH / 1080.0; 
    gl.uniform1f(gl.getUniformLocation(particleProgram, "u_pointSize"), currentPointSize * scalar);

    if (pointCount > 0) gl.drawArrays(gl.POINTS, 0, pointCount);
}

// UPDATED EXPORT FUNCTION: Handles both Download and POD
async function startTiledExport(mode = 'download') {
    if (isExporting || !currentCoeffs) return;
    isExporting = true;
    
    const inchesW = parseFloat(inpW.value);
    const inchesH = parseFloat(inpH.value);
    const dpi = parseInt(inpDPI.value);
    const passes = parseInt(inpPasses.value);
    
    const totalW = Math.floor(inchesW * dpi);
    const totalH = Math.floor(inchesH * dpi);
    
    const masterCanvas = document.createElement('canvas');
    masterCanvas.width = totalW;
    masterCanvas.height = totalH;
    const ctx = masterCanvas.getContext('2d');
    
    const TILE_SIZE = 2048; 
    const PADDING = 64; 
    const cols = Math.ceil(totalW / TILE_SIZE);
    const rows = Math.ceil(totalH / TILE_SIZE);
    
    const fbo = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
    const tex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, tex);
    // MOBILE FIX: Use 16F here as well
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA16F, TILE_SIZE+PADDING*2, TILE_SIZE+PADDING*2, 0, gl.RGBA, gl.HALF_FLOAT, null);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tex, 0);

    const resolveFbo = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, resolveFbo);
    const resolveTex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, resolveTex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, TILE_SIZE+PADDING*2, TILE_SIZE+PADDING*2, 0, gl.RGBA, gl.UNSIGNED_BYTE, null);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, resolveTex, 0);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.pixelStorei(gl.PACK_ALIGNMENT, 1); 

    const exportPhysics = currentPhysicsSteps;
    const resolutionRatio = totalH / canvas.height;
    
    // [FIX] Pass RAW opacity. renderTileParticles will divide by density.
    const exportOpacity = currentOpacity * resolutionRatio; 
    const exportJitter = currentJitter * resolutionRatio;

    // --- FIX: ASPECT RATIO PAN COMPENSATION ---
    // If the attractor is off-center and the user panned to fix it,
    // we must scale that pan amount to match the new aspect ratio.
    const screenAspect = canvas.width / canvas.height;
    const printAspect = totalW / totalH;
    
    // Generate ID for this export
    const meta = serializeState(); 
    const exportID = meta.id;

    for (let y = 0; y < rows; y++) {
        for (let x = 0; x < cols; x++) {
            
            const drawX = x * TILE_SIZE;
            const drawY = y * TILE_SIZE;
            const drawW = Math.min(TILE_SIZE, totalW - drawX);
            const drawH = Math.min(TILE_SIZE, totalH - drawY);
            
            const padLeft = (x > 0) ? PADDING : 0;
            const padTop  = (y > 0) ? PADDING : 0;
            const padRight = (x < cols - 1) ? PADDING : 0;
            const padBottom = (y < rows - 1) ? PADDING : 0;
            
            const tileW = drawW + padLeft + padRight;
            const tileH = drawH + padTop + padBottom;
            
            // 1. RENDER TO FLOAT FBO
            gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
            gl.viewport(0, 0, tileW, tileH);
            gl.clearColor(0,0,0,0);
            gl.clear(gl.COLOR_BUFFER_BIT);
            
            const globalX = drawX - padLeft;
            const globalY = drawY - padTop;
            const nX = globalX / totalW;
            const nY = 1.0 - ((globalY + tileH) / totalH);
            const nW = tileW / totalW;
            const nH = tileH / totalH;

            gl.useProgram(particleProgram);
      
            gl.enable(gl.BLEND);
            gl.blendFunc(gl.ONE, gl.ONE); 

            // Bind Sprite Texture
            gl.activeTexture(gl.TEXTURE0);
            gl.bindTexture(gl.TEXTURE_2D, gaussianTex);
            gl.uniform1i(gl.getUniformLocation(particleProgram, "u_sprite"), 0);
          
            for (let p = 0; p < passes; p++) {
                uiExport.innerText = `Tile ${y*cols + x + 1}/${rows*cols} - Pass ${p+1}/${passes}`;
                await new Promise(resolve => {
                    exportResolve = resolve;
                    worker.postMessage({ 
                        type: 'render', 
                        coeffs: currentCoeffs, 
                        physicsSteps: exportPhysics, 
                        density: currentDensity, 
                        seedOffset: p,
                        genType: currentGenType
                    });
                });
                renderTileParticles(totalW, totalH, [nX, nY, nW, nH], exportOpacity, totalW/totalH, exportJitter);
            }
            
            // 2. COMPOSITE TO RESOLVE FBO
            gl.bindFramebuffer(gl.DRAW_FRAMEBUFFER, resolveFbo);
            gl.viewport(0, 0, tileW, tileH);
            
            gl.useProgram(compositeProgram);
            gl.disable(gl.BLEND); 
            gl.activeTexture(gl.TEXTURE0);
            gl.bindTexture(gl.TEXTURE_2D, tex);
            gl.uniform1i(gl.getUniformLocation(compositeProgram, "u_tex"), 0);
            gl.uniform1f(gl.getUniformLocation(compositeProgram, "u_passes"), passes);
            gl.uniform1f(gl.getUniformLocation(compositeProgram, "u_gamma"), currentGamma);
            
            gl.uniform3f(gl.getUniformLocation(compositeProgram, "u_bg_a"), bgA[0], bgA[1], bgA[2]);
            gl.uniform3f(gl.getUniformLocation(compositeProgram, "u_bg_b"), bgB[0], bgB[1], bgB[2]);
            gl.uniform3f(gl.getUniformLocation(compositeProgram, "u_bg_params"), bgParams[0], bgParams[1], bgParams[2]);
            gl.uniform2f(gl.getUniformLocation(compositeProgram, "u_res"), totalW, totalH);
            gl.uniform2f(gl.getUniformLocation(compositeProgram, "u_off"), globalX, totalH - (globalY + tileH));
            gl.uniform1f(gl.getUniformLocation(compositeProgram, "u_noise"), currentNoise); 
            gl.uniform1i(gl.getUniformLocation(compositeProgram, "u_invert"), isInverted?1:0);
            gl.uniform1i(gl.getUniformLocation(compositeProgram, "u_inc_black"), incBlack?1:0);
            gl.uniform1i(gl.getUniformLocation(compositeProgram, "u_inc_white"), incWhite?1:0);
            gl.uniform1f(gl.getUniformLocation(compositeProgram, "u_scale"), 1.0); // <--- FORCE 1.0 FOR TILES
            
            let bMode = 0;
            if (blendMode === 'ADD') bMode = 1;
            else if (blendMode === 'SOFT') bMode = 2;
            gl.uniform1i(gl.getUniformLocation(compositeProgram, "u_blend_mode"), bMode);

            gl.drawArrays(gl.TRIANGLES, 0, 3);
            
            // 3. READ PIXELS
            const pixels = new Uint8Array(tileW * tileH * 4);
            gl.bindFramebuffer(gl.READ_FRAMEBUFFER, resolveFbo);
            gl.readPixels(0, 0, tileW, tileH, gl.RGBA, gl.UNSIGNED_BYTE, pixels);
            
            // 4. FLIP Y
            const flipped = new Uint8ClampedArray(tileW * tileH * 4);
            const rowBytes = tileW * 4;
            for (let r = 0; r < tileH; r++) {
                const srcRow = pixels.subarray(r * rowBytes, (r + 1) * rowBytes);
                const dstOffset = (tileH - r - 1) * rowBytes;
                flipped.set(srcRow, dstOffset);
            }
            
            // 5. PAINT
            const imgData = new ImageData(flipped, tileW, tileH);
            const tempCanvas = document.createElement('canvas');
            tempCanvas.width = tileW; 
            tempCanvas.height = tileH;
            const tempCtx = tempCanvas.getContext('2d');
            tempCtx.putImageData(imgData, 0, 0);
            ctx.drawImage(tempCanvas, padLeft, padTop, drawW, drawH, drawX, drawY, drawW, drawH);
        }
    }
    
    gl.deleteFramebuffer(fbo);
    gl.deleteFramebuffer(resolveFbo);
    gl.deleteTexture(tex);
    gl.deleteTexture(resolveTex);
    
    // --- BRANCH: DOWNLOAD OR POD ---
    if (mode === 'download') {
        uiExport.innerText = "Saving to disk...";
        masterCanvas.toBlob((blob) => {
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url; a.download = `attractor_${exportID}_${inchesW}x${inchesH}in_${dpi}dpi.png`; 
            document.body.appendChild(a); a.click(); document.body.removeChild(a);
            
            // AUTO-SAVE MATCHING JSON
            const data = { coeffs: Array.from(currentCoeffs), settings: meta };
            const jsonBlob = new Blob([JSON.stringify(data, null, 2)], {type: "application/json"});
            const jUrl = URL.createObjectURL(jsonBlob);
            const jA = document.createElement('a');
            jA.href = jUrl; jA.download = `attractor_${exportID}.json`;
            document.body.appendChild(jA); jA.click(); document.body.removeChild(jA);
            
            resetRenderState();
        }, 'image/png'); 

    } else if (mode === 'pod') {
        uiExport.innerText = "Preparing Upload...";
        masterCanvas.toBlob((blob) => {
             // Pass the blob to the POD function, DO NOT reset state yet (the async upload handles that)
             startPrintCheckout(blob);
        }, 'image/png');
    }
}

function resetRenderState() {
    isExporting = false;
    canvas.width = canvas.clientWidth; canvas.height = canvas.clientHeight;
    gl.viewport(0, 0, canvas.width, canvas.height);
    uiExport.innerText = "";
    worker.postMessage({ type: 'render', coeffs: currentCoeffs, physicsSteps: currentPhysicsSteps, density: currentDensity, genType: currentGenType });
}

function renderFrame() {
    // --- A. DYNAMIC RESOLUTION LOGIC ---
    // If dragging/interacting, drop to 50% scale (4x faster).
    // If idle, snap to 100% scale (Crisp).
    const targetScale = isInteracting ? 0.5 : 1.0;
    
    // Only resize if the target scale changed
    if (renderScale !== targetScale) {
        renderScale = targetScale;
        resizeViewportFBO();
    }
    
    // Safety check
    if (!viewFbo) resizeViewportFBO();

    // --- B. RENDER PARTICLES TO FBO ---
    gl.bindFramebuffer(gl.FRAMEBUFFER, viewFbo);
    
    // Important: Viewport matches the FBO size (e.g. 500x500), NOT the screen size (1000x1000)
    gl.viewport(0, 0, Math.floor(canvas.width * renderScale), Math.floor(canvas.height * renderScale));
    
    try {
        gl.clearColor(0,0,0,0);
        gl.clear(gl.COLOR_BUFFER_BIT);
        
        gl.useProgram(particleProgram);
        gl.enable(gl.BLEND);
        gl.blendFunc(gl.ONE, gl.ONE);
        
        // --- C. BIND OPTIMIZED GAUSSIAN SPRITE ---
        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, gaussianTex); // <--- The generated sprite
        gl.uniform1i(gl.getUniformLocation(particleProgram, "u_sprite"), 0);
        
        // --- D. UNIFORMS ---
        const rotMatrix = qToMatrix(currentQuat);
        gl.uniformMatrix4fv(gl.getUniformLocation(particleProgram, "u_rotation"), false, new Float32Array(rotMatrix));
        gl.uniform2f(gl.getUniformLocation(particleProgram, "u_pan"), camPanX, camPanY);
        gl.uniform1f(gl.getUniformLocation(particleProgram, "u_zoom"), camZoom);
        gl.uniform1f(gl.getUniformLocation(particleProgram, "u_aspect"), canvas.width / canvas.height);
        gl.uniform4f(gl.getUniformLocation(particleProgram, "u_tileBounds"), 0.0, 0.0, 1.0, 1.0);
        
        gl.uniform1i(gl.getUniformLocation(particleProgram, "u_colorMode"), colorMode);
        gl.uniform3f(gl.getUniformLocation(particleProgram, "u_colorSeed"), colorSeed[0], colorSeed[1], colorSeed[2]);
        gl.uniform1i(gl.getUniformLocation(particleProgram, "u_invert"), 0); 
        gl.uniform1i(gl.getUniformLocation(particleProgram, "u_inc_black"), 1); 
        gl.uniform1i(gl.getUniformLocation(particleProgram, "u_inc_white"), 1); 
        gl.uniform1f(gl.getUniformLocation(particleProgram, "u_variation"), currentVariation); 

        // --- DoF Settings ---
        const focusVal = (parseInt(document.getElementById('ui-focus').value) / 1000.0); 
        const focusSpanVal = (parseInt(document.getElementById('ui-focus-span').value) / 1000.0);
        const apertureVal = parseInt(document.getElementById('ui-aperture').value); // No divider (0-1000)

        gl.uniform1f(gl.getUniformLocation(particleProgram, "u_focusDist"), focusVal);
        gl.uniform1f(gl.getUniformLocation(particleProgram, "u_focusSpan"), focusSpanVal);
        gl.uniform1f(gl.getUniformLocation(particleProgram, "u_aperture"), apertureVal);

        // Opacity Math (Live View uses gpuRenderedDensity)
        let targetOpacity = currentOpacity / gpuRenderedDensity;
        let safeOpacity = Math.max(targetOpacity, 0.000001); 
        
        gl.uniform1f(gl.getUniformLocation(particleProgram, "u_opacity"), safeOpacity);
        gl.uniform1f(gl.getUniformLocation(particleProgram, "u_intensity"), currentIntensity);
        gl.uniform1f(gl.getUniformLocation(particleProgram, "u_noise"), currentNoise);
        gl.uniform2f(gl.getUniformLocation(particleProgram, "u_resolution"), canvas.width, canvas.height);
        gl.uniform1f(gl.getUniformLocation(particleProgram, "u_jitter"), currentJitter);
        gl.uniform1f(gl.getUniformLocation(particleProgram, "u_pointCount"), pointCount);

        const scalar = Math.max(1.0, canvas.height / 1080.0); 
        gl.uniform1f(gl.getUniformLocation(particleProgram, "u_pointSize"), currentPointSize * scalar);

        if (pointCount > 0) gl.drawArrays(gl.POINTS, 0, pointCount);

        // --- E. COMPOSITE TO SCREEN ---
        gl.bindFramebuffer(gl.FRAMEBUFFER, null);
        gl.viewport(0, 0, canvas.width, canvas.height); // Back to full screen size
        
        gl.useProgram(compositeProgram);
        gl.disable(gl.BLEND); 
        
        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, viewTex);
        gl.uniform1i(gl.getUniformLocation(compositeProgram, "u_tex"), 0);
        
        // PASS THE DYNAMIC SCALE
        // This tells the shader: "The texture is only 50% size, please stretch it"
        gl.uniform1f(gl.getUniformLocation(compositeProgram, "u_scale"), renderScale); 

        gl.uniform1f(gl.getUniformLocation(compositeProgram, "u_passes"), 1.0); 
        gl.uniform1f(gl.getUniformLocation(compositeProgram, "u_gamma"), currentGamma);
        gl.uniform3f(gl.getUniformLocation(compositeProgram, "u_bg_a"), bgA[0], bgA[1], bgA[2]);
        gl.uniform3f(gl.getUniformLocation(compositeProgram, "u_bg_b"), bgB[0], bgB[1], bgB[2]);
        gl.uniform3f(gl.getUniformLocation(compositeProgram, "u_bg_params"), bgParams[0], bgParams[1], bgParams[2]);
        gl.uniform2f(gl.getUniformLocation(compositeProgram, "u_res"), canvas.width, canvas.height);
        gl.uniform2f(gl.getUniformLocation(compositeProgram, "u_off"), 0, 0);
        gl.uniform1f(gl.getUniformLocation(compositeProgram, "u_noise"), currentNoise); 
        gl.uniform1i(gl.getUniformLocation(compositeProgram, "u_invert"), isInverted?1:0);
        gl.uniform1i(gl.getUniformLocation(compositeProgram, "u_inc_black"), incBlack?1:0);
        gl.uniform1i(gl.getUniformLocation(compositeProgram, "u_inc_white"), incWhite?1:0);
        
        let bMode = 0;
        if (blendMode === 'ADD') bMode = 1;
        else if (blendMode === 'SOFT') bMode = 2;
        gl.uniform1i(gl.getUniformLocation(compositeProgram, "u_blend_mode"), bMode);

        const checkGuide = document.getElementById('ui-show-guide');
        const valW = parseFloat(document.getElementById('ui-print-w').value) || 1;
        const valH = parseFloat(document.getElementById('ui-print-h').value) || 1;
        const printAspect = valW / valH;
        
        gl.uniform1i(gl.getUniformLocation(compositeProgram, "u_show_guide"), checkGuide.checked ? 1 : 0);
        gl.uniform1f(gl.getUniformLocation(compositeProgram, "u_print_aspect"), printAspect);
        
        gl.drawArrays(gl.TRIANGLES, 0, 3);
    } catch (e) {
        uiError.innerText = "Render Error: " + e.message;
    }
}

function loop() {
    if (!isExporting) { 
        if (canvas.width !== canvas.clientWidth || canvas.height !== canvas.clientHeight) {
            canvas.width = canvas.clientWidth; canvas.height = canvas.clientHeight;
            resizeViewportFBO(); 
        }
        renderFrame();
    }
    requestAnimationFrame(loop);
}
requestAnimationFrame(loop);
