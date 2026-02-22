/* STRANGE ATTRACTOR MINER - MAIN.JS
   Features: 
   - ODE Flow Mining (Chaos Hunter)
   - RK4 Integration (High Accuracy)
   - Mobile-Responsive UI (Sidebar vs Bottom Drawer)
   - High-DPI / Retina Display Support
   - 16-Bit Float Integration
   - Auto-Matching JSON/PNG Export
   - Multi-Engine: Poly, Symmetric, GRN, Dadras, Thomas, Aizawa, Rikitake, Chua, Hindmarsh-Rose, Moore-Spiegel
   - Power User Mode: Dynamic Search Bounds & Delta-Time control
   - POD (Print on Demand) Integration via Peecho + reCAPTCHA v3
*/

// ==========================================
// 0. CONSTANTS & GLOBAL STATE
// ==========================================
const MAX_ALLOCATION = 50000000; 

// --- POD CONFIG ---
const POD_API_URL = "https://script.google.com/macros/s/AKfycbyy2EWZpZ_LofW4JVHesxmaRq5LgPGrlEfKC49U2mVdujPhA0rr2XKqwlrn-vbFR7rt/exec"; 
const RECAPTCHA_SITE_KEY = "6Le8WWgsAAAAADEo9EQKpu_ZMaGaN0PHcCw0y4cL"; 

// --- DETECT MOBILE ---
const isMobile = window.matchMedia("(max-width: 768px)").matches;

// --- GENERATOR DEFINITIONS (METADATA) ---
const GEN_DEFS = {
    poly: { 
        label: "Polynomial", dt: 0.05, 
        params: [{ name: "Global Range (+/-)", idx: -1, min: 0.1, max: 5.0, valMin: 1.2, valMax: 1.2, defMin: 1.2, defMax: 1.2 }] 
    },
    rikitake: { 
        label: "Rikitake", dt: 0.015,
        params: [
            { name: "μ (Mu)", idx: 0, min: 0, max: 10, valMin: 1.0, valMax: 6.0, defMin: 1.0, defMax: 6.0 },
            { name: "a (Alpha)", idx: 1, min: 0, max: 20, valMin: 2.0, valMax: 8.0, defMin: 2.0, defMax: 8.0 }
        ]
    },
    chua: {
        label: "Chua", dt: 0.015,
        params: [
            { name: "α (Alpha)", idx: 0, min: 0, max: 50, valMin: 9.0, valMax: 18.0, defMin: 9.0, defMax: 18.0 },
            { name: "β (Beta)", idx: 1, min: 0, max: 50, valMin: 20.0, valMax: 35.0, defMin: 20.0, defMax: 35.0 },
            { name: "m0", idx: 2, min: -2, max: 0, valMin: -1.3, valMax: -0.9, defMin: -1.3, defMax: -0.9 },
            { name: "m1", idx: 3, min: -2, max: 0, valMin: -0.9, valMax: -0.5, defMin: -0.9, defMax: -0.5 }
        ]
    },
    hindmarsh: {
        label: "Hindmarsh-Rose", dt: 0.01,
        params: [
            { name: "r (Burst Rate)", idx: 4, min: 0.001, max: 0.05, valMin: 0.001, valMax: 0.01, defMin: 0.001, defMax: 0.01 },
            { name: "I (Current)", idx: 7, min: 0, max: 10, valMin: 2.9, valMax: 3.4, defMin: 2.9, defMax: 3.4 }
        ]
    },
    moore: {
        label: "Moore-Spiegel", dt: 0.0002,
        params: [
            { name: "Γ (Gamma)", idx: 0, min: 0, max: 100, valMin: 20.0, valMax: 35.0, defMin: 20.0, defMax: 35.0 },
            { name: "R (Reynolds)", idx: 1, min: 0, max: 200, valMin: 80.0, valMax: 120.0, defMin: 80.0, defMax: 120.0 }
        ]
    },
    dadras: { 
        label: "Dadras", dt: 0.015, 
        params: [
            { name: "p", idx: 0, min: 1, max: 5, valMin: 2.5, valMax: 3.5, defMin: 2.5, defMax: 3.5 },
            { name: "σ (Sigma)", idx: 1, min: 1, max: 5, valMin: 2.0, valMax: 3.5, defMin: 2.0, defMax: 3.5 },
            { name: "ρ (Rho)", idx: 2, min: 0, max: 5, valMin: 1.5, valMax: 2.5, defMin: 1.5, defMax: 2.5 },
            { name: "β (Beta)", idx: 3, min: 0, max: 5, valMin: 1.5, valMax: 2.5, defMin: 1.5, defMax: 2.5 },
            { name: "γ (Gamma)", idx: 4, min: 5, max: 15, valMin: 8.0, valMax: 10.0, defMin: 8.0, defMax: 10.0 }
        ] 
    },
    thomas: { 
        label: "Thomas", dt: 0.05, 
        params: [
            { name: "b (Dissipation)", idx: 0, min: 0, max: 0.4, valMin: 0.18, valMax: 0.22, defMin: 0.18, defMax: 0.22 }
        ] 
    },
    aizawa: { 
        label: "Aizawa", dt: 0.01, 
        params: [
            { name: "ε (Epsilon)", idx: 0, min: 0, max: 2, valMin: 0.90, valMax: 1.0, defMin: 0.90, defMax: 1.0 },
            { name: "α (Alpha)", idx: 1, min: 0, max: 2, valMin: 0.65, valMax: 0.75, defMin: 0.65, defMax: 0.75 },
            { name: "γ (Gamma)", idx: 2, min: 0, max: 2, valMin: 0.55, valMax: 0.65, defMin: 0.55, defMax: 0.65 },
            { name: "δ (Delta)", idx: 3, min: 0, max: 5, valMin: 3.25, valMax: 3.75, defMin: 3.25, defMax: 3.75 },
            { name: "β (Beta)", idx: 4, min: 0, max: 1, valMin: 0.20, valMax: 0.30, defMin: 0.20, defMax: 0.30 },
            { name: "ζ (Zeta)", idx: 5, min: 0, max: 0.5, valMin: 0.075, valMax: 0.125, defMin: 0.075, defMax: 0.125 }
        ] 
    },
    sym: { label: "Symmetric", dt: 0.015, params: [] },
    grn: { label: "GRN", dt: 0.015, params: [] }
};

// --- GLOBAL STATE VARIABLES ---
let pointCount = 0;
let gpuRenderedDensity = 1.0;
let gaussianTex = null;
let currentConstraints = null;

// Defaults (Will be tuned below if mobile)
let camZoom = 2.0; 
let camPanX = 0, camPanY = 0;
let currentQuat = [0, 0, 0, 1]; 
let currentCoeffs = null;
let currentGenType = 'poly';
let colorMode = 0;
let colorSeed = [0.0, 0.0, 0.0];

// Defaults for Glow Mode
let blendMode = 'ADD';        
let isInverted = false;      
let incBlack = true;
let incWhite = true;

// Default Background to Dark/Black
let bgA = [0.0, 0.0, 0.0];    
let bgB = [0.1, 0.1, 0.1];    
let bgParams = [0.0, 0.5, 1.0]; 

// Performance Defaults
let currentPhysicsSteps = 200000; 
let currentOpacity = 0.06;          
let currentIntensity = 2.0;  
let currentGamma = 1.8;
let currentNoise = 0.05;      
let currentPointSize = 1.0;
let currentDensity = 12;
let currentJitter = 0.5;
let currentVariation = 0.0;
let isExporting = false; 

// Export State
let exportUnit = 'inches'; // 'inches' or 'pixels'
let exportTransparent = false; 

// --- MOBILE SAFE DEFAULTS ---
if (isMobile) {
    currentPhysicsSteps = 200000; // Cap at 250k for safety
    currentPointSize = 2.5;       // Thicker lines for small high-DPI screens
    camZoom = 1.5;                // Zoom out slightly
    currentDensity = 12;           // Keep density low
}

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
let renderScale = 1.0;      
let isInteracting = false; 
let pendingResize = false; 

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
        if (e.data.type === 'mine') mine(e.data.genType, e.data.constraints);
        else if (e.data.type === 'mutate') mutate(e.data.coeffs, e.data.genType);
        else if (e.data.type === 'render') renderExisting(e.data);
    };

    function catmullRom(p0, p1, p2, p3, t) {
        const v0 = (p2 - p0) * 0.5;
        const v1 = (p3 - p1) * 0.5;
        const t2 = t * t;
        const t3 = t * t2;
        return (2 * p1 - 2 * p2 + v0 + v1) * t3 + 
               (-3 * p1 + 3 * p2 - 2 * v0 - v1) * t2 + 
               v0 * t + p1;
    }

    function safePostMessage(msg, transfers) {
        if (transfers) {
            for (let i = 0; i < transfers.length; i++) {
                if (!transfers[i]) {
                    console.error("Worker: Attempted to transfer null/undefined buffer.");
                    return; 
                }
            }
        }
        self.postMessage(msg, transfers);
    }

    function renderExisting(data) {
        const c = new Float32Array(data.coeffs);
        // data.constraints might pass dtOverride if saved
        const dtOverride = (data.constraints && data.constraints.dt) ? data.constraints.dt : null;
        
        const result = generateTrace(data.physicsSteps, data.density, data.seedOffset||0, c, data.genType, dtOverride); 
        safePostMessage({
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

    function mine(genType, constraints) {
        let attempts = 0;
        let lastReport = Date.now();
        
        while(true) {
            attempts++;
            let coeffs;
            
            // --- DYNAMIC PARAMETER GENERATION (POWER MODE) ---
            if (constraints && constraints.params && constraints.params.length > 0) {
                let size = 30; // Poly default
                if (genType === 'rikitake' || genType === 'moore') size = 2;
                if (genType === 'chua') size = 4;
                if (genType === 'hindmarsh') size = 8;
                if (genType === 'dadras') size = 5;
                
                coeffs = new Float32Array(size);

                // Pre-fill defaults for complex constants (e.g. Hindmarsh)
                if (genType === 'hindmarsh') {
                    coeffs[0]=1; coeffs[1]=3; coeffs[2]=1; coeffs[3]=5; coeffs[5]=4; coeffs[6]=-1.6;
                }

                for(let p of constraints.params) {
                    if (p.idx === -1) { 
                        // Poly global symmetric range
                        for(let i=0; i<30; i++) coeffs[i] = (Math.random() * (p.valMax*2)) - p.valMax;
                    } else {
                        // Specific Parameter
                        coeffs[p.idx] = p.valMin + Math.random() * (p.valMax - p.valMin);
                    }
                }
            } 
            else {
                // --- STANDARD DEFAULTS ---
                if (genType === 'sym') {
                    coeffs = new Float32Array(10); for(let i=0; i<10; i++) coeffs[i] = getRandomSprott(); 
                } 
                else if (genType === 'grn') {
                    coeffs = new Float32Array(18);
                    for(let i=0; i<9; i++) coeffs[i] = (Math.random() * 20) - 10; 
                    coeffs[9] = (coeffs[0]+coeffs[1]+coeffs[2])*0.5 + (Math.random()-0.5);
                    coeffs[10] = (coeffs[3]+coeffs[4]+coeffs[5])*0.5 + (Math.random()-0.5);
                    coeffs[11] = (coeffs[6]+coeffs[7]+coeffs[8])*0.5 + (Math.random()-0.5);
                    for(let i=12; i<15; i++) coeffs[i] = 2.5 + Math.random() * 4.5; 
                    let d1 = 0.2+Math.random()*0.2; let d2 = 0.6+Math.random()*0.3; let d3 = 1.0+Math.random()*0.5; 
                    let r = Math.random();
                    if(r<0.33) { coeffs[15]=d1; coeffs[16]=d2; coeffs[17]=d3; }
                    else if(r<0.66) { coeffs[15]=d3; coeffs[16]=d1; coeffs[17]=d2; }
                    else { coeffs[15]=d2; coeffs[16]=d3; coeffs[17]=d1; }
                }
                else if (genType === 'dadras') {
                    coeffs = new Float32Array(5);
                    coeffs[0] = 2.5 + Math.random(); coeffs[1] = 2.0 + Math.random() * 1.5; 
                    coeffs[2] = 1.5 + Math.random(); coeffs[3] = 1.5 + Math.random(); coeffs[4] = 8.0 + Math.random() * 2;   
                }
                else if (genType === 'thomas') {
                    coeffs = new Float32Array(1); coeffs[0] = 0.18 + Math.random() * 0.04; 
                }
                else if (genType === 'aizawa') {
                    coeffs = new Float32Array(6);
                    coeffs[0]=0.95+(Math.random()-0.5)*0.1; coeffs[1]=0.7+(Math.random()-0.5)*0.1;
                    coeffs[2]=0.6+(Math.random()-0.5)*0.1; coeffs[3]=3.5+(Math.random()-0.5)*0.5;
                    coeffs[4]=0.25+(Math.random()-0.5)*0.1; coeffs[5]=0.1+(Math.random()-0.5)*0.05;
                }
                else if (genType === 'rikitake') {
                    coeffs = new Float32Array(2);
                    coeffs[0] = 1.0 + (Math.random() * 6.0); 
                    coeffs[1] = 2.0 + (Math.random() * 8.0); 
                }
                else if (genType === 'chua') {
                    coeffs = new Float32Array(4);
                    coeffs[0] = 8.0 + Math.random() * 12.0;  
                    coeffs[1] = 20.0 + Math.random() * 20.0; 
                    coeffs[2] = -1.143 + (Math.random() - 0.5) * 0.5; 
                    coeffs[3] = -0.714 + (Math.random() - 0.5) * 0.5; 
                }
                else if (genType === 'hindmarsh') {
                    coeffs = new Float32Array(8);
                    coeffs[0]=1.0; coeffs[1]=3.0; coeffs[2]=1.0; coeffs[3]=5.0; 
                    coeffs[4] = 0.001 + (Math.random() * 0.01); 
                    coeffs[5]=4.0; coeffs[6]=-1.6;
                    coeffs[7] = 2.0 + Math.random() * 2.0; 
                }
                else if (genType === 'moore') {
                    coeffs = new Float32Array(2);
                    coeffs[0] = 5.0 + (Math.random() * 45.0); 
                    coeffs[1] = 50.0 + (Math.random() * 100.0); 
                }
                else {
                    coeffs = new Float32Array(30);
                    for(let i=0; i<30; i++) coeffs[i] = (Math.random() * 2.4) - 1.2;
                }
            }
            
            let dtOverride = (constraints && constraints.dt) ? constraints.dt : null;

            if (checkChaosTail(coeffs, genType, dtOverride)) {
                const result = generateTrace(50000, 1, 0, coeffs, genType, dtOverride);
                safePostMessage({
                    type: 'found', 
                    source: 'mine', 
                    coeffs: coeffs, 
                    genType: genType,
                    constraints: constraints, // Pass back so we can save/reload state
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
            } 
            else if (genType === 'grn') {
                const range = (idx < 9) ? 2.0 : 0.5;
                child[idx] += (Math.random() - 0.5) * range;
                if (idx >= 12) child[idx] = Math.abs(child[idx]) + 0.1; 
            }
            else if (genType === 'rikitake') {
                child[idx] += (Math.random() - 0.5) * 0.5;
            }
            else if (genType === 'chua') {
                child[idx] += (Math.random() - 0.5) * 1.0;
            }
            else if (genType === 'hindmarsh') {
                if (idx === 7) child[idx] += (Math.random() - 0.5) * 0.1;
                else if (idx === 4) child[idx] += (Math.random() - 0.5) * 0.001;
            }
            else if (genType === 'moore') {
                child[idx] += (Math.random() - 0.5) * 2.0; 
            }
            else if (genType === 'dadras') {
                child[idx] += (Math.random() - 0.5) * 0.1;
            }
            else if (genType === 'thomas') {
                child[0] += (Math.random() - 0.5) * 0.01;
            }
            else if (genType === 'aizawa') {
                child[idx] += (Math.random() - 0.5) * 0.02;
            }
            else {
                child[idx] += (Math.random() - 0.5) * 0.1;
            }

            // Mutation inherits base parameters, so we pass null dtOverride 
            // (it will assume default dt unless we architect complex inheritance)
            if (checkChaosTail(child, genType, null)) {
                const result = generateTrace(50000, 1, 0, child, genType, null);
                safePostMessage({
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

    function checkChaosTail(c, genType, dtOverride) {
        let x, y, z;
        if (genType === 'sym') { x = 0.1; y = 0.0; z = -0.1; } 
        else if (genType === 'grn') { x = Math.random(); y = Math.random(); z = Math.random(); }
        else if (genType === 'dadras') { x = 1.1; y = 2.1; z = -1.5; }
        else if (genType === 'thomas') { 
            x = (Math.random() - 0.5) * 3.0; 
            y = (Math.random() - 0.5) * 3.0; 
            z = (Math.random() - 0.5) * 3.0; 
        }
        else if (genType === 'aizawa') { x = 0.1; y = 0.0; z = 0.0; }
        else if (genType === 'rikitake') { x = 1.0; y = 0.0; z = 1.0; }
        else if (genType === 'chua') { x = 0.1; y = 0.0; z = 0.0; }
        else if (genType === 'hindmarsh') { x = -1.0; y = 0.0; z = 0.0; }
        else if (genType === 'moore') { x = 0.1; y = 0.0; z = 0.0; }
        else { x = 0.05; y = 0.05; z = 0.05; }

        let sx = x + 0.000001, sy = y, sz = z;
        
        let dt = 0.015;
        if (genType === 'poly') dt = 0.05;
        if (genType === 'aizawa') dt = 0.01;
        if (genType === 'moore') dt = 0.0002; 
        if (genType === 'hindmarsh') dt = 0.02;

        if (dtOverride) dt = dtOverride;

        let p = c; 
        
        let a0,a1,a2,a3,a4,a5,a6,a7,a8,a9;
        let b0,b1,b2,b3,b4,b5,b6,b7,b8,b9;
        let c0,c1,c2,c3,c4,c5,c6,c7,c8,c9;
        if (genType === 'poly') {
           a0=c[0]; a1=c[1]; a2=c[2]; a3=c[3]; a4=c[4]; a5=c[5]; a6=c[6]; a7=c[7]; a8=c[8]; a9=c[9];
           b0=c[10]; b1=c[11]; b2=c[12]; b3=c[13]; b4=c[14]; b5=c[15]; b6=c[16]; b7=c[17]; b8=c[18]; b9=c[19];
           c0=c[20]; c1=c[21]; c2=c[22]; c3=c[23]; c4=c[24]; c5=c[25]; c6=c[26]; c7=c[27]; c8=c[28]; c9=c[29];
        }

        function calcD(px, py, pz, res) {
            if (genType === 'sym') {
                res.dx = p[0] + p[1]*px + p[2]*py + p[3]*pz + p[4]*px*px + p[5]*py*py + p[6]*pz*pz + p[7]*px*py + p[8]*px*pz + p[9]*py*pz;
                res.dy = p[0] + p[1]*py + p[2]*pz + p[3]*px + p[4]*py*py + p[5]*pz*pz + p[6]*px*px + p[7]*py*pz + p[8]*py*px + p[9]*pz*px;
                res.dz = p[0] + p[1]*pz + p[2]*px + p[3]*py + p[4]*pz*pz + p[5]*px*px + p[6]*py*py + p[7]*pz*px + p[8]*pz*py + p[9]*px*py;
            } 
            else if (genType === 'grn') {
                let s1 = p[0]*px + p[1]*py + p[2]*pz - p[9];
                let s2 = p[3]*px + p[4]*py + p[5]*pz - p[10];
                let s3 = p[6]*px + p[7]*py + p[8]*pz - p[11];
                let act1 = 1.0 / (1.0 + Math.exp(-p[12] * s1));
                let act2 = 1.0 / (1.0 + Math.exp(-p[13] * s2));
                let act3 = 1.0 / (1.0 + Math.exp(-p[14] * s3));
                res.dx = act1 - p[15]*px;
                res.dy = act2 - p[16]*py;
                res.dz = act3 - p[17]*pz;
            }
            else if (genType === 'rikitake') {
                res.dx = -p[0]*px + py*pz;
                res.dy = -p[0]*py + px*(p[1] - pz);
                res.dz = 1.0 - px*py;
            }
            else if (genType === 'chua') {
                let h = p[3]*px + 0.5*(p[2]-p[3])*(Math.abs(px+1.0) - Math.abs(px-1.0));
                res.dx = p[0]*(py - px - h);
                res.dy = px - py + pz;
                res.dz = -p[1]*py;
            }
            else if (genType === 'hindmarsh') {
                let x2 = px*px;
                let x3 = x2*px;
                res.dx = py - p[0]*x3 + p[1]*x2 - pz + p[7];
                res.dy = p[2] - p[3]*x2 - py;
                res.dz = p[4]*(p[5]*(px - p[6]) - pz);
            }
            else if (genType === 'moore') {
                res.dx = py;
                res.dy = pz;
                res.dz = -pz - (p[0] - p[1] + p[1]*px*px)*py - p[0]*px;
            }
            else if (genType === 'dadras') {
                res.dx = py - p[0]*px + p[1]*py*pz;
                res.dy = p[2]*py - px*pz + pz;
                res.dz = p[3]*px*py - p[4]*pz;
            }
            else if (genType === 'thomas') {
                res.dx = Math.sin(py) - p[0]*px;
                res.dy = Math.sin(pz) - p[0]*py;
                res.dz = Math.sin(px) - p[0]*pz;
            }
            else if (genType === 'aizawa') {
                let x2 = px*px; let y2 = py*py;
                res.dx = (pz - p[1])*px - p[3]*py;
                res.dy = p[3]*px + (pz - p[1])*py;
                res.dz = p[2] + p[0]*pz - (pz*pz*pz)/3.0 - (x2+y2)*(1.0 + p[4]*pz) + p[5]*pz*px*px*px;
            }
            else {
                res.dx = a0 + a1*px + a2*py + a3*pz + a4*px*px + a5*py*py + a6*pz*pz + a7*px*py + a8*px*pz + a9*py*pz;
                res.dy = b0 + b1*px + b2*py + b3*pz + b4*px*px + b5*py*py + b6*pz*pz + b7*px*py + b8*px*pz + b9*py*pz;
                res.dz = c0 + c1*px + c2*py + c3*pz + c4*px*px + c5*py*py + c6*pz*pz + c7*px*py + c8*px*pz + c9*py*pz;
            }
        }
        
        let k1={dx:0,dy:0,dz:0}, k2={dx:0,dy:0,dz:0}, k3={dx:0,dy:0,dz:0}, k4={dx:0,dy:0,dz:0};
        
        let settleSteps = (genType === 'thomas') ? 5000 : 1500;
        if(genType === 'moore') settleSteps = 5000;

        for(let i=0; i<settleSteps; i++) {
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
        
        let voxRes = 0.5;
        if(genType === 'grn') voxRes = 0.05;
        if(genType === 'thomas') voxRes = 0.2;
        if(genType === 'aizawa') voxRes = 0.1;
        if(genType === 'moore') voxRes = 0.5;

        const visited = new Set();
        let steps = 3000;
        
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

            if (i % 5 === 0) visited.add(Math.floor(x/voxRes)+","+Math.floor(y/voxRes)+","+Math.floor(z/voxRes));
        }
        
        let lyapunov = lyapunovSum / steps;
        
        let minL = 0.001;
        let minVol = 25;
        let minWidth = 1.0;

        if (genType === 'grn') { minL=0.0015; minWidth=0.05; minVol=60; }
        if (genType === 'thomas') { minL=0.001; minWidth=0.5; minVol=25; } 
        if (genType === 'aizawa') { minL=0.0001; minWidth=0.5; minVol=30; }
        if (genType === 'moore') { minL=0.002; minWidth=2.0; minVol=50; }

        if (lyapunov < minL) return false;
        
        let wX = maxX - minX, wY = maxY - minY, wZ = maxZ - minZ;
        if ((wX + wY + wZ) < minWidth) return false;
        if (visited.size < minVol) return false;
        if (lyapunov > 2.0) return false;

        return true;
    }

    function generateTrace(nSteps, density, seedOffset, c, genType, dtOverride) {
        const totalPoints = nSteps * density; 
        
        let posData = new Float32Array(totalPoints * 3);
        let metaData = new Float32Array(totalPoints * 2); 
        const rand = mulberry32(seedOffset + 12345);

        let x, y, z;
        if (genType === 'sym') { x = 0.1; y = 0.0; z = -0.1; } 
        else if (genType === 'grn') { x = 0.1; y = 0.1; z = 0.1; }
        else if (genType === 'dadras') { x = 1.1; y = 2.1; z = -1.5; }
        else if (genType === 'thomas') { 
            x = (rand() - 0.5) * 3.0; 
            y = (rand() - 0.5) * 3.0; 
            z = (rand() - 0.5) * 3.0;
        }
        else if (genType === 'aizawa') { x = 0.1; y = 0.0; z = 0.0; }
        else if (genType === 'rikitake') { x = 1.0; y = 0.0; z = 1.0; }
        else if (genType === 'chua') { x = 0.1; y = 0.0; z = 0.0; }
        else if (genType === 'hindmarsh') { x = -1.0; y = 0.0; z = 0.0; }
        else if (genType === 'moore') { x = 0.1; y = 0.0; z = 0.0; }
        else { x = 0.05; y = 0.05; z = 0.05; }

        let dt = 0.015;
        if (genType === 'poly') dt = 0.05;
        if (genType === 'aizawa') dt = 0.01;
        if (genType === 'moore') dt = 0.0002; 
        if (genType === 'hindmarsh') dt = 0.02;

        if (dtOverride) dt = dtOverride;

        // Cache Coeffs
        let p = c; 
        
        let a0,a1,a2,a3,a4,a5,a6,a7,a8,a9;
        let b0,b1,b2,b3,b4,b5,b6,b7,b8,b9;
        let c0,c1,c2,c3,c4,c5,c6,c7,c8,c9;
        if (genType === 'poly') {
           a0=c[0]; a1=c[1]; a2=c[2]; a3=c[3]; a4=c[4]; a5=c[5]; a6=c[6]; a7=c[7]; a8=c[8]; a9=c[9];
           b0=c[10]; b1=c[11]; b2=c[12]; b3=c[13]; b4=c[14]; b5=c[15]; b6=c[16]; b7=c[17]; b8=c[18]; b9=c[19];
           c0=c[20]; c1=c[21]; c2=c[22]; c3=c[23]; c4=c[24]; c5=c[25]; c6=c[26]; c7=c[27]; c8=c[28]; c9=c[29];
        }

        function calcD(px, py, pz, res) {
            if (genType === 'sym') {
                res.dx = p[0] + p[1]*px + p[2]*py + p[3]*pz + p[4]*px*px + p[5]*py*py + p[6]*pz*pz + p[7]*px*py + p[8]*px*pz + p[9]*py*pz;
                res.dy = p[0] + p[1]*py + p[2]*pz + p[3]*px + p[4]*py*py + p[5]*pz*pz + p[6]*px*px + p[7]*py*pz + p[8]*py*px + p[9]*pz*px;
                res.dz = p[0] + p[1]*pz + p[2]*px + p[3]*py + p[4]*pz*pz + p[5]*px*px + p[6]*py*py + p[7]*pz*px + p[8]*pz*py + p[9]*px*py;
            } 
            else if (genType === 'grn') {
                let s1 = p[0]*px + p[1]*py + p[2]*pz - p[9];
                let s2 = p[3]*px + p[4]*py + p[5]*pz - p[10];
                let s3 = p[6]*px + p[7]*py + p[8]*pz - p[11];
                let act1 = 1.0 / (1.0 + Math.exp(-p[12] * s1));
                let act2 = 1.0 / (1.0 + Math.exp(-p[13] * s2));
                let act3 = 1.0 / (1.0 + Math.exp(-p[14] * s3));
                res.dx = act1 - p[15]*px;
                res.dy = act2 - p[16]*py;
                res.dz = act3 - p[17]*pz;
            }
            else if (genType === 'rikitake') {
                res.dx = -p[0]*px + py*pz;
                res.dy = -p[0]*py + px*(p[1] - pz);
                res.dz = 1.0 - px*py;
            }
            else if (genType === 'chua') {
                let h = p[3]*px + 0.5*(p[2]-p[3])*(Math.abs(px+1.0) - Math.abs(px-1.0));
                res.dx = p[0]*(py - px - h);
                res.dy = px - py + pz;
                res.dz = -p[1]*py;
            }
            else if (genType === 'hindmarsh') {
                let x2 = px*px;
                let x3 = x2*px;
                res.dx = py - p[0]*x3 + p[1]*x2 - pz + p[7];
                res.dy = p[2] - p[3]*x2 - py;
                res.dz = p[4]*(p[5]*(px - p[6]) - pz);
            }
            else if (genType === 'moore') {
                res.dx = py;
                res.dy = pz;
                res.dz = -pz - (p[0] - p[1] + p[1]*px*px)*py - p[0]*px;
            }
            else if (genType === 'dadras') {
                res.dx = py - p[0]*px + p[1]*py*pz;
                res.dy = p[2]*py - px*pz + pz;
                res.dz = p[3]*px*py - p[4]*pz;
            }
            else if (genType === 'thomas') {
                res.dx = Math.sin(py) - p[0]*px;
                res.dy = Math.sin(pz) - p[0]*py;
                res.dz = Math.sin(px) - p[0]*pz;
            }
            else if (genType === 'aizawa') {
                let x2 = px*px; let y2 = py*py;
                res.dx = (pz - p[1])*px - p[3]*py;
                res.dy = p[3]*px + (pz - p[1])*py;
                res.dz = p[2] + p[0]*pz - (pz*pz*pz)/3.0 - (x2+y2)*(1.0 + p[4]*pz) + p[5]*pz*px*px*px;
            }
            else {
                res.dx = a0 + a1*px + a2*py + a3*pz + a4*px*px + a5*py*py + a6*pz*pz + a7*px*py + a8*px*pz + a9*py*pz;
                res.dy = b0 + b1*px + b2*py + b3*pz + b4*px*px + b5*py*py + b6*pz*pz + b7*px*py + b8*px*pz + b9*py*pz;
                res.dz = c0 + c1*px + c2*py + c3*pz + c4*px*px + c5*py*py + c6*pz*pz + c7*px*py + c8*px*pz + c9*py*pz;
            }
        }
        
        let k1={dx:0,dy:0,dz:0}, k2={dx:0,dy:0,dz:0}, k3={dx:0,dy:0,dz:0}, k4={dx:0,dy:0,dz:0};
        
        const history = [];

        function stepPhysics() {
            calcD(x, y, z, k1);
            let speed = Math.sqrt(k1.dx*k1.dx + k1.dy*k1.dy + k1.dz*k1.dz);
            let logVel = Math.log(speed + 1.0);
            let curv = 0; 
            
            calcD(x + k1.dx*dt*0.5, y + k1.dy*dt*0.5, z + k1.dz*dt*0.5, k2);
            calcD(x + k2.dx*dt*0.5, y + k2.dy*dt*0.5, z + k2.dz*dt*0.5, k3);
            calcD(x + k3.dx*dt, y + k3.dy*dt, z + k3.dz*dt, k4);
            
            let nextX = x + (k1.dx + 2*k2.dx + 2*k3.dx + k4.dx)*(dt/6);
            let nextY = y + (k1.dy + 2*k2.dy + 2*k3.dy + k4.dy)*(dt/6);
            let nextZ = z + (k1.dz + 2*k2.dz + 2*k3.dz + k4.dz)*(dt/6);
            
            x = nextX; y = nextY; z = nextZ;
            return { x:x, y:y, z:z, vel: logVel, curv: curv };
        }

        // 1. Settle
        for(let i=0; i<2000; i++) stepPhysics();
        for(let i=0; i<4; i++) history.push(stepPhysics());

        let outIdx = 0;
        let maxVel = 0;

        for(let i=0; i<nSteps; i++) {
            let p0 = history[0]; let p1 = history[1]; let p2 = history[2]; let p3 = history[3];
            if (Math.abs(p2.x) > 1000 || isNaN(p2.x)) break;

            for(let d=0; d<density; d++) {
                if (outIdx >= totalPoints) break;
                let t = d / density; 
                let px = catmullRom(p0.x, p1.x, p2.x, p3.x, t);
                let py = catmullRom(p0.y, p1.y, p2.y, p3.y, t);
                let pz = catmullRom(p0.z, p1.z, p2.z, p3.z, t);
                
                let pVel = p1.vel + (p2.vel - p1.vel) * t;
                let pCurv = p1.curv + (p2.curv - p1.curv) * t;
                if (pVel > maxVel) maxVel = pVel;

                posData[outIdx*3] = px; posData[outIdx*3+1] = py; posData[outIdx*3+2] = pz;
                metaData[outIdx*2] = pVel; metaData[outIdx*2+1] = pCurv;
                outIdx++;
            }
            history.shift(); history.push(stepPhysics());
        }
        
        if (maxVel > 0) {
            for(let i=0; i<totalPoints; i++) metaData[i*2] /= maxVel;
        }
        
        // Auto-Centering & Scaling
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
        
        // GRN needs specific zoom (lives in [0,1])
        let scaleTarget = 0.5;
        if (genType === 'grn') scaleTarget = 0.8;
        if (genType === 'dadras') scaleTarget = 0.35; // Large range
        if (genType === 'thomas') scaleTarget = 0.5;
        if (genType === 'aizawa') scaleTarget = 0.6;
        if (genType === 'chua') scaleTarget = 0.2; // Chua is tiny
        if (genType === 'hindmarsh') scaleTarget = 0.25; 
        if (genType === 'moore') scaleTarget = 0.8; 
        if (genType === 'rikitake') scaleTarget = 2.0;

        if (rms > 0) {
            let s = scaleTarget / rms;
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

gaussianTex = createGaussianTexture();

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

uniform float u_focusDist; 
uniform float u_focusSpan; 
uniform float u_aperture;
uniform float u_opacity;

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

    vec4 rotated = u_rotation * vec4(a_position, 1.0);
    float dist = rotated.z; 

    float zNoise = (hash(float(gl_VertexID) + 77.7) - 0.5) * 0.04;
    float blurredDist = dist + zNoise;

    rotated.xy *= u_zoom;
    rotated.xy += u_pan; 
    rotated.x /= u_aspect;
    
    vec2 normalizedPos = (rotated.xy * 0.5) + 0.5;
    vec2 tilePos = (normalizedPos - u_tileBounds.xy) / u_tileBounds.zw;
    
    gl_Position = vec4(tilePos * 2.0 - 1.0, 0.0, 1.0);
    
    float distFromFocus = abs(blurredDist - u_focusDist);
    float effectiveDist = max(0.0, distFromFocus - u_focusSpan);

    float rawBlur = (effectiveDist * effectiveDist) * (u_aperture * 0.1);
    
    float scaleFactor = u_resolution.y / 1080.0;
    float scaledBlur = rawBlur * scaleFactor;
    scaledBlur = clamp(scaledBlur, 0.0, 400.0 * scaleFactor);

    // FIX 1: Allow baseSize to be sub-pixel (remove the max 1.0 clamp here)
    float baseSize = u_pointSize; 
    float targetSize = baseSize + scaledBlur;
    
    float rnd = hash(float(gl_VertexID) + 555.555);
    float sizeJitter = (rnd - 0.5) * (targetSize * 0.15); 
    
    float finalSize = targetSize + sizeJitter;
    
    // WebGL forces this to 1.0 minimum physically
    gl_PointSize = max(1.0, finalSize); 
    
    // FIX 2: Calculate Area Attenuation (Depth of Field)
    // Protected against divide-by-zero
    float areaRatio = baseSize / max(0.0001, targetSize); 

    // FIX 3: Sub-Pixel Energy Conservation
    // If the logical size is 0.5 but we draw 1.0, we have 4x too much area.
    // We must dim opacity by 0.5*0.5 = 0.25 to maintain correct energy.
    float subPixelComp = 1.0;
    if (finalSize < 1.0) {
        subPixelComp = finalSize * finalSize;
    }
    
    v_attenuation = areaRatio * subPixelComp;
    
    float jx = (hash(float(gl_VertexID)) - 0.5) * u_jitter;
    float jy = (hash(float(gl_VertexID) + 123.45) - 0.5) * u_jitter;
    vec2 pixelSize = 2.0 / u_resolution; 
    gl_Position.xy += vec2(jx, jy) * pixelSize;

    // OPTIMIZATION: Vertex Culling
    // If the particle is effectively invisible due to extreme depth blur and low global opacity,
    // throw it outside the clip space (off-screen) so the rasterizer ignores it completely.
    if (v_attenuation * u_opacity < 0.000001) {
        gl_Position = vec4(2.0, 2.0, 2.0, 1.0); 
    }
}`;

const fsSource = `#version 300 es
precision highp float;

in float v_vel;
in float v_curv;
in float v_time;
in vec3 v_pos;
in float v_attenuation; 

// REMOVED u_sprite uniform
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
    // --- PROCEDURAL SPLAT MATH ---
    // Calculate distance from the center of the point (0.0 to 0.5)
    float dist = length(gl_PointCoord - vec2(0.5));
    
    // Circular crop: if outside the circle, discard immediately
    if (dist > 0.5) discard; 
    
    // Smooth Gaussian-like falloff using exponential decay
    float shapeAlpha = exp(-dist * dist * 16.0); 

    // Calculate final opacity early
    float finalOpac = u_opacity * v_attenuation * shapeAlpha;
    
    // OPTIMIZATION: Discard nearly invisible fragments before doing complex color math
    if (finalOpac < 0.000001) discard;

    float val = 0.0;
    if (u_colorMode == 0) val = pow(v_vel, 0.5); 
    else if (u_colorMode == 2) val = pow(v_curv, 0.5); 
    else if (u_colorMode == 3) {
        vec3 p = v_pos + (u_colorSeed * 10.0);
        float pattern = sin(p.x * 5.0) + cos(p.y * 5.0) + sin(p.z * 5.0);
        val = 0.5 + (pattern * 0.25);
    } 
    val = clamp(val, 0.0, 1.0);
    
    vec3 phaseDrift = vec3(sin(v_time*3.0), cos(v_time*5.0), sin(v_time*7.0)) * (u_variation*0.2); 
    vec3 phaseShift = u_colorSeed + phaseDrift;

    vec3 rgb;
    if (u_colorMode == 0) rgb = palette(val, vec3(0.5), vec3(0.5), vec3(1.0), vec3(0.0, 0.33, 0.67) + phaseShift);
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
uniform float u_scale; 
uniform int u_transparent;

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
    vec2 globalUV = (gl_FragCoord.xy + u_off) / u_res;
    ivec2 texCoord = ivec2(gl_FragCoord.xy * u_scale);
    vec4 pVal = texelFetch(u_tex, texCoord, 0) / u_passes;

    vec3 color = pVal.rgb;
    vec3 energy = color;

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

    // --- TRANSPARENCY EXPORT LOGIC ---
    if (u_transparent == 1) {
        if (u_invert == 1) {
             // INK MODE (Normal): Output Black with Alpha = Darkness
             float alpha = dot(energy, vec3(0.299, 0.587, 0.114)); 
             c = vec4(0.0, 0.0, 0.0, clamp(alpha, 0.0, 1.0));
        } else {
             // GLOW MODE (Add): Fix for "Insubstantial" look
             // 1. Alpha is determined by the BRIGHTEST color channel
             float maxComp = max(energy.r, max(energy.g, energy.b));
             float alpha = clamp(maxComp, 0.0, 1.0);
             
             // 2. Un-multiply RGB so the PNG viewer restores full brightness
             vec3 safeRGB = vec3(0.0);
             if(alpha > 0.0001) {
                 safeRGB = energy / alpha;
             }
             
             c = vec4(clamp(safeRGB, 0.0, 1.0), alpha);
        }
    } else {
        // --- STANDARD BACKGROUND COMPOSITE ---
        float aspect = u_res.x / u_res.y;
        vec2 auv = globalUV; auv.x *= aspect;
        vec2 cnt = vec2(0.5*aspect,0.5) + (vec2(cos(u_bg_params.x), sin(u_bg_params.x))*u_bg_params.y);
        float t = smoothstep(0.0, 1.5*u_bg_params.z, distance(auv, cnt));
        vec3 bg = mix(u_bg_a, u_bg_b, t);
        
        // FIXED: USE GLOBAL COORDS FOR NOISE TO PREVENT TILE SEAMS
        bg += (hash(gl_FragCoord.xy + u_off)-0.5) * (u_noise*0.05);

        vec3 finalRGB;
        if (u_invert == 1) finalRGB = bg * (1.0 - energy); 
        else finalRGB = bg + energy; 
        
        if (u_show_guide == 1) {
            float viewAspect = u_res.x / u_res.y;
            float targetAspect = u_print_aspect;
            
            // 1. Print is NARROWER than Screen (Vertical Bars / Pillarbox)
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
            // 2. Print is WIDER than Screen (Horizontal Bars / Letterbox)
            else if (targetAspect > viewAspect) {
                float safeRatio = viewAspect / targetAspect;
                float margin = (1.0 - safeRatio) * 0.5;
                if (globalUV.y < margin || globalUV.y > (1.0 - margin)) {
                    finalRGB *= 0.3; 
                    finalRGB += vec3(0.1, 0.1, 0.1); 
                }
                float lineH = 1.0 / u_res.y; 
                if (abs(globalUV.y - margin) < lineH || abs(globalUV.y - (1.0-margin)) < lineH) {
                    finalRGB = vec3(0.5, 0.5, 0.5); 
                }
            }
        }
        
        c = vec4(finalRGB, 1.0);
        float ditherStrength = max(u_noise * 0.2, 0.004); 
        
        // FIXED: USE GLOBAL COORDS FOR DITHER TOO
        c.rgb += triangularNoise(gl_FragCoord.xy + u_off) * ditherStrength;
    }
}`;

function createShader(gl, type, source) {
    const shader = gl.createShader(type);
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
        console.error("SHADER COMPILE ERROR:", gl.getShaderInfoLog(shader));
        gl.deleteShader(shader);
        return null;
    }
    return shader;
}

function createGaussianTexture() {
    const size = 64; 
    const canvas = document.createElement('canvas');
    canvas.width = size;
    canvas.height = size;
    const ctx = canvas.getContext('2d');

    const cx = size / 2;
    const cy = size / 2;
    const radius = size / 2;
    
    const grad = ctx.createRadialGradient(cx, cy, 0, cx, cy, radius);
    grad.addColorStop(0, 'rgba(255, 255, 255, 1.0)'); 
    grad.addColorStop(0.2, 'rgba(255, 255, 255, 0.8)');
    grad.addColorStop(0.5, 'rgba(255, 255, 255, 0.2)');
    grad.addColorStop(1, 'rgba(0, 0, 0, 0.0)');

    ctx.fillStyle = grad;
    ctx.fillRect(0, 0, size, size);

    const tex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, canvas);
    
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
    
    // Scale texture to device physical pixels * render scale
    const w = Math.floor(canvas.width * renderScale);
    const h = Math.floor(canvas.height * renderScale);

    gl.bindTexture(gl.TEXTURE_2D, viewTex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA16F, w, h, 0, gl.RGBA, gl.HALF_FLOAT, null);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

    gl.bindFramebuffer(gl.FRAMEBUFFER, viewFbo);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, viewTex, 0);
    
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
        id: generateID(),
        timestamp: Date.now(),
        genType: currentGenType,
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
        },
        export: {
            unit: exportUnit, 
            transparent: exportTransparent, 
            width: document.getElementById('ui-print-w').value,
            height: document.getElementById('ui-print-h').value,
            dpi: document.getElementById('ui-print-dpi').value,
            passes: document.getElementById('ui-print-passes').value
        }
    };
}

function applyState(data) {
    if (!data) return;
    
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
        if (s.dof) {
            document.getElementById('ui-focus').value = s.dof.focus || 0;
            document.getElementById('ui-aperture').value = s.dof.aperture || 0;
            document.getElementById('ui-focus-span').value = s.dof.span || 0;
        }
        if (s.export) {
            const ex = s.export;
            exportUnit = ex.unit || 'inches';
            exportTransparent = ex.transparent || false; 
            document.getElementById('ui-export-unit').value = exportUnit;
            document.getElementById('ui-export-transparent').checked = exportTransparent;
            
            // Fix Labels
            const labelW = document.getElementById('ui-label-w');
            const labelH = document.getElementById('ui-label-h');
            if(exportUnit === 'pixels') {
                labelW.innerText = "Width (px)"; labelH.innerText = "Height (px)";
            } else {
                labelW.innerText = "Width (in)"; labelH.innerText = "Height (in)";
            }

            document.getElementById('ui-print-w').value = ex.width || (exportUnit==='pixels'?1920:24);
            document.getElementById('ui-print-h').value = ex.height || (exportUnit==='pixels'?1080:36);
            document.getElementById('ui-print-dpi').value = ex.dpi || (exportUnit==='pixels'?72:300);
            document.getElementById('ui-print-passes').value = ex.passes || 1;
        }

    } else {
        camZoom = 2.0; camPanX=0; camPanY=0; currentQuat=qIdentity();
    }
    
    selectGenType.value = currentGenType;
    selectBlend.value = blendMode;
    selectColor.value = colorMode;
    checkIncBlack.checked = incBlack;
    checkIncWhite.checked = incWhite;
    inputBg1.value = rgbToHex(bgA[0],bgA[1],bgA[2]);
    inputBg2.value = rgbToHex(bgB[0],bgB[1],bgB[2]);
    
    sliderLength.value = currentPhysicsSteps;
    sliderDensity.value = currentDensity;
    sliderOpacity.value = currentOpacity * 1000;
    sliderIntensity.value = currentIntensity * 10;
    sliderSize.value = currentPointSize * 10;
    sliderGamma.value = currentGamma * 10;
    sliderJitter.value = currentJitter * 10;
    sliderSmooth.value = currentNoise * 200;
    sliderVariation.value = currentVariation * 100;
}

//
function writePngDpi(blob, dpi) {
    // 1. Convert DPI to Pixels Per Meter (PNG standard)
    // 1 inch = 0.0254 meters
    const pixelsPerMeter = Math.round(dpi / 0.0254);
    
    return new Promise(resolve => {
        const reader = new FileReader();
        reader.onload = (e) => {
            const buffer = e.target.result;
            const view = new DataView(buffer);
            
            // Check for PNG signature
            if (view.getUint32(0) !== 0x89504E47) {
                resolve(blob); // Not a PNG, return original
                return;
            }

            // 2. Create the pHYs chunk
            // Length (4) + Type (4) + Data (9) + CRC (4) = 21 bytes
            const physChunk = new Uint8Array(21);
            const physView = new DataView(physChunk.buffer);
            
            physView.setUint32(0, 9); // Length of data
            // Type: "pHYs" (0x70485973)
            physChunk.set([112, 72, 89, 115], 4); 
            
            physView.setUint32(8, pixelsPerMeter); // X axis
            physView.setUint32(12, pixelsPerMeter); // Y axis
            physChunk[16] = 1; // Unit specifier: 1 = meter
            
            // Calculate CRC for Type + Data (bytes 4 to 16)
            const crcInput = physChunk.subarray(4, 17);
            const crc = crc32(crcInput);
            physView.setUint32(17, crc);

            // 3. Construct new Blob: Signature + IHDR + pHYs + rest
            // IHDR is always the first chunk, 13 bytes data + 12 bytes overhead = 25 bytes
            // plus 8 byte signature = 33 bytes total for start.
            const newBlob = new Blob([
                buffer.slice(0, 33), // Sig + IHDR
                physChunk,           // Insert pHYs
                buffer.slice(33)     // Rest of file (IDAT, IEND)
            ], { type: 'image/png' });

            resolve(newBlob);
        };
        reader.readAsArrayBuffer(blob);
    });
}

// CRC32 Table-less implementation (compact)
function crc32(buf) {
    let crc = -1;
    for (let i = 0; i < buf.length; i++) {
        crc ^= buf[i];
        for (let j = 0; j < 8; j++) {
            crc = (crc >>> 1) ^ ((crc & 1) ? 0xEDB88320 : 0);
        }
    }
    return (crc ^ -1) >>> 0;
}

// ==========================================
// 4. UI GENERATION (RESPONSIVE)
// ==========================================

// 1. Inject CSS for Responsive Layout
const style = document.createElement('style');
style.textContent = `
    /* Hide reCAPTCHA Badge but stay compliant */
    .grecaptcha-badge { visibility: hidden; }

    /* --- DESKTOP (DEFAULT) --- */
    #colorControls {
        position: absolute;
        top: 20px;
        right: 20px;
        width: 240px;
        max-height: 90vh;
        background: rgba(0,0,0,0.85);
        border: 1px solid #0f0;
        font-family: monospace;
        overflow-y: auto;
        z-index: 999;
        display: block; /* Visible by default on desktop */
        backdrop-filter: blur(5px);
    }

    /* Section Headers */
    .ui-header {
        background: rgba(34, 34, 34, 0.85); /* Semi-transparent */
        backdrop-filter: blur(5px);
        color: #fff;
        padding: 5px 10px;
        cursor: pointer;
        user-select: none;
        border-bottom: 1px solid #444;
        font-weight: bold;
    }

    .ui-content {
        padding: 10px;
        display: none; /* Collapsed by default */
    }

    /* Toggle Button (Gear) - Fixed Bottom Right */
    #ui-toggle-btn {
        position: absolute;
        bottom: 20px;
        right: 20px;
        background: #000;
        color: #0f0;
        border: 1px solid #0f0;
        width: 40px;
        height: 40px;
        font-size: 20px;
        cursor: pointer;
        z-index: 1000;
        display: flex;
        align-items: center;
        justify-content: center;
        user-select: none;
        box-shadow: 2px 2px 5px rgba(0,0,0,0.5);
    }

/* Quick Search Button - Fixed Bottom Right (Next to Gear) */
    #ui-quick-search-btn {
        position: absolute;
        bottom: 20px;
        right: 70px; /* 20px margin + 40px gear + 10px gap */
        background: #0f0;
        color: #000;
        border: 1px solid #0f0;
        height: 40px;
        padding: 0 20px;
        font-size: 16px;
        font-weight: bold;
        font-family: monospace;
        cursor: pointer;
        z-index: 1000;
        display: flex;
        align-items: center;
        justify-content: center;
        box-shadow: 2px 2px 5px rgba(0,0,0,0.5);
        transition: background 0.2s, color 0.2s;
    }
    #ui-quick-search-btn:hover {
        background: #fff;
        color: #000;
        border-color: #fff;
    }

    /* --- MOBILE OVERRIDES --- */
    @media (max-width: 600px) {
        #colorControls {
            top: auto;
            bottom: 0;
            left: 0;
            right: 0;
            width: 100%;
            max-height: 45vh; /* Bottom sheet */
            border: none;
            border-top: 2px solid #0f0;
        }

        /* Float the button on mobile too */
        #ui-toggle-btn {
            width: 50px;
            height: 50px;
            border-radius: 50%;
        }

        /* Fat Finger Targets */
        input[type=range], button, select {
            min-height: 35px; 
            margin-bottom: 5px;
        }

#ui-quick-search-btn {
            height: 50px;
            right: 80px; /* 20px margin + 50px gear + 10px gap */
            border-radius: 25px; /* Pill shape to match circular gear */
        }
    }
`;
document.head.appendChild(style);

// 2. Create UI Containers
const toggleBtn = document.createElement('div');
toggleBtn.id = 'ui-toggle-btn';
toggleBtn.innerText = '⚙️';
toggleBtn.onclick = (e) => {
    e.stopPropagation(); // Prevent canvas click
    const panel = document.getElementById('colorControls');
    const quickBtn = document.getElementById('ui-quick-search-btn'); // Target the new button
    
    if (panel.style.display === 'none') {
        panel.style.display = 'block';
        if(quickBtn) quickBtn.style.display = 'none'; // Hide quick search when panel opens
    } else {
        panel.style.display = 'none';
        if(quickBtn) quickBtn.style.display = 'flex'; // Show quick search when panel closes
    }
};
document.body.appendChild(toggleBtn);

// --- ADD THE QUICK SEARCH BUTTON ---
const quickSearchBtn = document.createElement('button');
quickSearchBtn.id = 'ui-quick-search-btn';
// On mobile, the panel is hidden, so show the button. On desktop, panel is open, so hide it.
quickSearchBtn.style.display = isMobile ? 'flex' : 'none';
quickSearchBtn.innerText = '⛏️ SEARCH!';
quickSearchBtn.onclick = (e) => {
    e.stopPropagation(); // Prevent canvas interaction
    
    // Update UI Status (mirroring the main MINE button behavior)
    const uiStatus = document.getElementById('ui-main-status');
    if(uiStatus) {
        uiStatus.innerText = "Scanning..."; 
        uiStatus.style.color = "#ffff00"; 
    }
    
    // Trigger the worker to mine a new attractor
    worker.postMessage({
        type: 'mine', 
        genType: currentGenType, 
        constraints: currentConstraints
    }); 
};
document.body.appendChild(quickSearchBtn);
// -----------------------------------

const div = document.createElement('div');
div.id = 'colorControls';
// On mobile, maybe start hidden so they see the art first?
if (isMobile) div.style.display = 'none'; 

// --- POWER MODE MODAL ---
const powerModal = document.createElement('div');
powerModal.id = 'power-modal';
powerModal.style.cssText = `
    display: none; position: fixed; top: 50%; left: 50%; transform: translate(-50%, -50%);
    width: 90%; max-width: 400px; background: rgba(10,10,10,0.95); border: 1px solid #0f0;
    padding: 20px; z-index: 10000; box-shadow: 0 0 20px rgba(0,255,0,0.2);
    font-family: monospace; color: #fff; max-height: 80vh; overflow-y: auto; backdrop-filter: blur(5px);
`;
document.body.appendChild(powerModal);

const overlay = document.createElement('div');
overlay.id = 'power-overlay';
overlay.style.cssText = `
    display: none; position: fixed; top: 0; left: 0; width: 100%; height: 100%;
    background: rgba(0,0,0,0.7); z-index: 9999; backdrop-filter: blur(2px);
`;
overlay.onclick = () => { powerModal.style.display = 'none'; overlay.style.display = 'none'; };
document.body.appendChild(overlay);

function updatePowerUI() {
    const type = selectGenType.value;
    const defs = GEN_DEFS[type];
    
    let html = `<div style="display:flex; justify-content:space-between; align-items:center;">
        <h3 style="margin:0; color:#0f0">⚡ Power Mode: ${defs.label}</h3>
        <button id="pm-reset" style="background:#444; color:#fff; border:1px solid #777; padding:2px 8px; font-size:10px; cursor:pointer;">RESET DEFAULTS</button>
    </div>`;
    html += `<div style="margin-bottom:15px; font-size:12px; color:#aaa;">Define search ranges for the miner.</div>`;
    
    html += `<div style="margin-bottom:10px;">
        <label>Time Step (dt)</label>
        <input type="number" id="pm-dt" value="${defs.dt}" step="0.0001" style="width:100%; background:#222; color:#fff; border:1px solid #444; padding:5px;">
    </div>`;

    if (defs.params.length === 0) {
        html += `<div style="color:#666; padding:10px; border:1px dashed #444;">No configurable parameters for this generator yet.</div>`;
    } else {
        defs.params.forEach((p, i) => {
            html += `<div style="border-top:1px solid #333; padding: 10px 0;">
                <div style="color:#fff; font-weight:bold; margin-bottom:5px;">${p.name}</div>
                <div style="display:flex; gap:10px; align-items:center;">
                    <input type="number" id="pm-min-${i}" value="${p.valMin}" step="0.001" style="width:45%; background:#222; color:#fff; border:1px solid #444; padding:5px;">
                    <span style="color:#888;">to</span>
                    <input type="number" id="pm-max-${i}" value="${p.valMax}" step="0.001" style="width:45%; background:#222; color:#fff; border:1px solid #444; padding:5px;">
                </div>
            </div>`;
        });
    }

    html += `<button id="pm-apply" style="width:100%; padding:12px; background:#0f0; color:#000; font-weight:bold; border:none; margin-top:15px; cursor:pointer;">APPLY & MINE ⛏️</button>`;
    
    powerModal.innerHTML = html;

    // Reset Handler
    document.getElementById('pm-reset').onclick = () => {
        defs.dt = (type === 'moore') ? 0.0002 : (type === 'poly' || type === 'thomas' ? 0.05 : 0.015);
        defs.params.forEach(p => {
            p.valMin = p.defMin;
            p.valMax = p.defMax;
        });
        updatePowerUI(); // Refresh UI
    };

    // Apply Handler
    document.getElementById('pm-apply').onclick = () => {
        // 1. Update Persistent State
        defs.dt = parseFloat(document.getElementById('pm-dt').value);
        defs.params.forEach((p, i) => {
            p.valMin = parseFloat(document.getElementById(`pm-min-${i}`).value);
            p.valMax = parseFloat(document.getElementById(`pm-max-${i}`).value);
        });

        // 2. Build Constraints object for Worker
        const newParams = defs.params.map(p => ({
            idx: p.idx,
            valMin: p.valMin,
            valMax: p.valMax
        }));

        currentConstraints = { dt: defs.dt, params: newParams };
        powerModal.style.display = 'none';
        overlay.style.display = 'none';
        
        uiStatus.innerText = "Mining with Constraints...";
        worker.postMessage({ type: 'mine', genType: type, constraints: currentConstraints });
    };
}

function createSection(title, contentHTML) {
    const section = document.createElement('div');
    const header = document.createElement('div');
    header.className = 'ui-header';
    header.innerHTML = `▶ ${title}`;
    
    const content = document.createElement('div');
    content.className = 'ui-content';
    content.innerHTML = contentHTML;
    
    // FIX: Explicitly set display state so toggle works on first click
    content.style.display = 'none'; 

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
        <option value="grn">Gene Regulatory Network</option>
        <option value="dadras">Dadras (Complex Butterfly)</option>
        <option value="thomas">Thomas (Cyclic Lattice)</option>
        <option value="aizawa">Aizawa (Sphere Tube)</option>
        <option value="rikitake">Rikitake (Double Spiral)</option>
        <option value="chua">Chua (Double Scroll)</option>
        <option value="hindmarsh">Hindmarsh-Rose (Neuron)</option>
        <option value="moore">Moore-Spiegel (Cosmic Knot)</option>
    </select>
    <div style="display:flex; gap:5px; margin-bottom:5px;">
        <button id="ui-btn-mine" style="flex:1; cursor:pointer; background:#440000; color:#fff; border:1px solid #f00; padding:10px;">⛏️ MINE</button>
        <button id="ui-btn-mutate" style="flex:1; cursor:pointer; background:#440044; color:#fff; border:1px solid #f0f; padding:10px;">🧬 MUTATE</button>
    </div>
    <button id="ui-btn-power" style="width:100%; background:#222; color:#0f0; border:1px dashed #0f0; padding:5px; margin-bottom:10px; cursor:pointer;">⚡ POWER MODE</button>

    <div style="display:flex; gap:5px; margin-bottom:10px;">
        <button id="ui-btn-save" style="flex:1; background:#222; color:#fff; border:1px solid #555; padding:5px;">💾 SAVE</button>
        <button id="ui-btn-load" style="flex:1; background:#222; color:#fff; border:1px solid #555; padding:5px;">📂 LOAD</button>
        <input type="file" id="ui-file-input" style="display:none" accept=".json">
    </div>
    <div style="color:#0f0; margin-bottom:5px;">LENGTH (Trail)</div>
    <input type="range" id="ui-length" min="10000" max="2000000" step="10000" value="200000" style="width:100%;">
    <div style="color:#0f0; margin-bottom:5px;">SMOOTHNESS (Precision)</div>
    <input type="range" id="ui-density" min="1" max="50" value="12" style="width:100%;">
    <div style="color:#0f0; margin-bottom:5px;">VARIATION (Texture)</div>
    <input type="range" id="ui-variation" min="0" max="100" value="0" style="width:100%;">
`));

div.appendChild(createSection("APPEARANCE", `
    <div style="color:#0f0; margin-bottom:5px;">BLEND MODE</div>
    <select id="ui-blend-mode" style="width:100%; margin-bottom:10px;">
        <option value="NORMAL">Normal (Paint)</option>
        <option value="ADD">Additive (Glow)</option>
    </select>
    <div style="color:#0f0; margin-bottom:5px;">COLOR MODE</div>
    <select id="ui-color-mode" style="width:100%; margin-bottom:10px;">
        <option value="0">Velocity (Magma)</option>
        <option value="2">Curvature (Electric)</option>
        <option value="3">Spatial (Marble)</option>
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
    <input type="range" id="ui-intensity" min="0" max="150" value="20" style="width:100%;">
    
    <div style="color:#0f0; margin-bottom:5px;">OPACITY (Alpha)</div>
    <input type="range" id="ui-opacity" min="1" max="300" value="60" style="width:100%;">
    
    <div style="color:#0f0; margin-bottom:5px;">SIZE (Thickness)</div>
    <input type="range" id="ui-size" min="1" max="50" value="10" style="width:100%;">
    <div style="color:#0f0; margin-bottom:5px;">GAMMA (Contrast)</div>
    <input type="range" id="ui-gamma" min="1" max="35" value="18" style="width:100%;">
    <div style="color:#0f0; margin-bottom:5px;">JITTER (Anti-Alias)</div>
    <input type="range" id="ui-jitter" min="0" max="20" value="5" style="width:100%;">
    <div style="color:#0f0; margin-bottom:5px;">DEBAND (Dither)</div>
    <input type="range" id="ui-deband" min="0" max="100" value="10" style="width:100%;">
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

    <div style="margin-bottom:5px; font-size:12px; color:#aaa;">Dimensions</div>
    
    <select id="ui-export-unit" style="width:100%; margin-bottom:10px; background:#333; color:#fff; border:1px solid #555;">
        <option value="inches">Inches (Print)</option>
        <option value="pixels">Pixels (Screen)</option>
    </select>

    <div style="display:flex; gap:5px; align-items:flex-end; margin-bottom:10px;">
        <div style="flex:1">
            <span id="ui-label-w" style="color:#0f0; font-size:10px;">Width (in)</span>
            <input type="number" id="ui-print-w" value="24" style="width:100%">
        </div>
        <button id="ui-btn-swap" style="width:24px; height:24px; margin-bottom:2px; background:#222; color:#fff; border:1px solid #555; cursor:pointer; padding:0; line-height:1;">⇄</button>
        <div style="flex:1">
            <span id="ui-label-h" style="color:#0f0; font-size:10px;">Height (in)</span>
            <input type="number" id="ui-print-h" value="36" style="width:100%">
        </div>
    </div>
    
    <label style="display:block; margin-bottom:10px; color:#aaa; font-size:11px; cursor:pointer;">
        <input type="checkbox" id="ui-export-transparent"> Transparent Background
    </label>

    <div style="margin-bottom:5px; font-size:12px; color:#aaa;">Quality</div>
    <div style="display:flex; gap:5px; margin-bottom:10px;">
        <div style="flex:1">
            <span style="color:#0f0; font-size:10px;">DPI</span>
            <input type="number" id="ui-print-dpi" value="300" style="width:100%">
        </div>
        <div style="flex:1">
            <span style="color:#0f0; font-size:10px;">Passes</span>
            <input type="number" id="ui-print-passes" value="1" style="width:100%">
        </div>
    </div>

    <button id="ui-btn-snap" style="width:100%; background:#fff; color:#000; border:none; padding:10px; cursor:pointer; font-weight:bold; margin-bottom:5px;">📸 SNAPSHOT</button>
    <button id="ui-btn-order" style="width:100%; background:#0f0; color:#000; border:none; padding:10px; cursor:pointer; font-weight:bold;">🛒 ORDER PRINT</button>
    <div id="ui-export-status" style="color:#fff; font-size:10px; margin-top:5px; text-align:center;"></div>
`));

const footer = document.createElement('div');
footer.innerHTML = `
    <div id="ui-main-status" style="color:#ffff00; font-size:12px; margin-top:10px; text-align:center;">Initialized</div>
    <div id="ui-error-log" style="color:red; font-size:10px; margin-top:5px;"></div>
    <div style="font-size:10px; color:#777; margin-top:10px; line-height:1.2;">
        This site is protected by reCAPTCHA and the Google <a href="https://policies.google.com/privacy" target="_blank" style="color:#999">Privacy Policy</a> and <a href="https://policies.google.com/terms" target="_blank" style="color:#999">Terms of Service</a> apply.
    </div>
`;
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

selectBlend.value = blendMode; 
selectColor.value = colorMode;

// --- EVENT HANDLERS ---
window.onerror = function(msg, url, line) {
    uiError.innerText = `JS Error: ${msg} (Line ${line})`;
    return false;
};

// Unit & Transparent Toggle Logic
const selectExportUnit = document.getElementById('ui-export-unit');
const checkTransparent = document.getElementById('ui-export-transparent');
const labelW = document.getElementById('ui-label-w');
const labelH = document.getElementById('ui-label-h');
const btnSwap = document.getElementById('ui-btn-swap');

btnSwap.onclick = () => {
    const temp = inpW.value;
    inpW.value = inpH.value;
    inpH.value = temp;
};

selectExportUnit.onchange = (e) => {
    exportUnit = e.target.value;
    const wInput = document.getElementById('ui-print-w');
    const hInput = document.getElementById('ui-print-h');
    const dpiInput = document.getElementById('ui-print-dpi');

    if (exportUnit === 'pixels') {
        labelW.innerText = "Width (px)";
        labelH.innerText = "Height (px)";
        wInput.value = 1920; 
        hInput.value = 1080;
        dpiInput.value = 72; 
    } else {
        labelW.innerText = "Width (in)";
        labelH.innerText = "Height (in)";
        wInput.value = 24;
        hInput.value = 36;
        dpiInput.value = 300; 
    }
};

checkTransparent.onchange = (e) => { exportTransparent = e.target.checked; };

selectGenType.onchange = (e) => { 
    currentGenType = e.target.value;
    currentConstraints = null; // Reset constraints on type switch
};

inputBg1.oninput = (e) => { bgA = hexToRgb(e.target.value); };
inputBg2.oninput = (e) => { bgB = hexToRgb(e.target.value); };
btnRerollBg.onclick = () => {
    if (blendMode === 'ADD') {
        const dominant = Math.floor(Math.random() * 3); 
        const intensity = 0.1 + Math.random() * 0.15; 
        const c = [0, 0, 0];
        c[dominant] = intensity; 
        c[(dominant + 1) % 3] = Math.random() * (intensity * 0.3);
        c[(dominant + 2) % 3] = Math.random() * (intensity * 0.3);
        bgA = c; 
        bgB = [0.02, 0.02, 0.02]; 
    } else {
        const val = 0.8 + Math.random() * 0.2; 
        bgA = [val, val, val]; 
        const r = val - (Math.random() * 0.1);
        const g = val - (Math.random() * 0.1);
        const b = val - (Math.random() * 0.1);
        bgB = [r,g,b];
    }

    inputBg1.value = rgbToHex(bgA[0], bgA[1], bgA[2]);
    inputBg2.value = rgbToHex(bgB[0], bgB[1], bgB[2]);
    bgParams[0] = Math.random() * 6.28; 
    bgParams[1] = 0.2 + Math.random() * 0.6; 
    bgParams[2] = 0.8 + Math.random() * 0.5; 
};
selectBlend.onchange = (e) => { 
    blendMode = e.target.value; 
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

// Bind Power Mode
document.getElementById('ui-btn-power').onclick = () => {
    updatePowerUI();
    powerModal.style.display = 'block';
    overlay.style.display = 'block';
};

btnMine.onclick = () => { 
    uiStatus.innerText = "Scanning..."; 
    uiStatus.style.color = "#ffff00"; 
    worker.postMessage({type: 'mine', genType: currentGenType, constraints: currentConstraints}); 
};

btnMutate.onclick = () => {
    if (!currentCoeffs) { uiStatus.innerText = "Mine first!"; return; }
    uiStatus.innerText = "Mutating..."; uiStatus.style.color = "#ff00ff";
    worker.postMessage({ type: 'mutate', coeffs: currentCoeffs, genType: currentGenType });
};

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

// ==========================================
// CORRECTED EVENT LISTENERS FOR HIGH-DPI
// ==========================================
function getRelativeTouchPos(touch, element) {
    const rect = element.getBoundingClientRect();
    return {
        x: touch.clientX - rect.left,
        y: touch.clientY - rect.top,
        cx: rect.width / 2,
        cy: rect.height / 2,
        width: rect.width,
        height: rect.height
    };
}

canvas.addEventListener('touchstart', (e) => {
    if(isExporting) return;
    isInteracting = true;
    e.preventDefault();
    if (e.touches.length === 1) {
        const p = getRelativeTouchPos(e.touches[0], canvas);
        lastX = p.x; lastY = p.y;
        
        const dist = Math.hypot(p.x - p.cx, p.y - p.cy);
        const maxDist = Math.min(p.width, p.height) / 2;
        
        if (dist > maxDist * 0.8) {
            isRolling = true;
            isDragging = false;
            lastAngle = Math.atan2(p.y - p.cy, p.x - p.cx);
        } else {
            isRolling = false;
            isDragging = true;
        }
        isPanning = false;
    } else if (e.touches.length === 2) {
        isDragging = false; isPanning = true; isRolling = false;
        const p1 = getRelativeTouchPos(e.touches[0], canvas);
        const p2 = getRelativeTouchPos(e.touches[1], canvas);
        lastDist = Math.hypot(p1.x - p2.x, p1.y - p2.y);
        lastX = (p1.x + p2.x) / 2;
        lastY = (p1.y + p2.y) / 2;
    }
}, {passive:false});

canvas.addEventListener('touchmove', (e) => {
    if(isExporting) return;
    e.preventDefault();
    if (e.touches.length === 1) {
        const p = getRelativeTouchPos(e.touches[0], canvas);
        
        if (isRolling) {
            const angle = Math.atan2(p.y - p.cy, p.x - p.cx);
            const delta = angle - lastAngle;
            applyRoll(delta);
            lastAngle = angle;
            lastX = p.x; lastY = p.y;
        } else if (isDragging) {
            const dx = p.x - lastX;
            const dy = p.y - lastY;
            applyTrackballRotation(dx, dy);
            lastX = p.x; lastY = p.y;
        }
    } else if (e.touches.length === 2 && isPanning) {
        const p1 = getRelativeTouchPos(e.touches[0], canvas);
        const p2 = getRelativeTouchPos(e.touches[1], canvas);
        
        const newDist = Math.hypot(p1.x - p2.x, p1.y - p2.y);
        const cx = (p1.x + p2.x) / 2;
        const cy = (p1.y + p2.y) / 2;
        
        const zoomDelta = newDist / lastDist;
        camZoom *= zoomDelta;
        
        const pdx = (cx - lastX) / p1.height * 2;
        const pdy = (cy - lastY) / p1.height * 2;
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
    const rect = canvas.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;
    const cx = rect.width / 2;
    const cy = rect.height / 2;
    
    lastX = mx; lastY = my;
    
    if (e.button === 2) {
        isPanning = true; isDragging = false; isRolling = false;
    } else {
        isPanning = false;
        const dist = Math.hypot(mx - cx, my - cy);
        const maxDist = Math.min(rect.width, rect.height) / 2;
        
        if (dist > maxDist * 0.8) {
            isRolling = true;
            isDragging = false;
            lastAngle = Math.atan2(my - cy, mx - cx);
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
    const rect = canvas.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;
    
    if (isRolling) {
        const cx = rect.width / 2;
        const cy = rect.height / 2;
        const angle = Math.atan2(my - cy, mx - cx);
        const delta = angle - lastAngle;
        applyRoll(delta);
        lastAngle = angle;
    } else if (isDragging) {
        const dx = mx - lastX;
        const dy = my - lastY;
        applyTrackballRotation(dx, dy);
    } else if (isPanning) {
        const dx = mx - lastX;
        const dy = my - lastY;
        camPanX += (dx / rect.height) * 2;
        camPanY -= (dy / rect.height) * 2;
    }
    lastX = mx; lastY = my;
};

canvas.onwheel = (e) => {
    isInteracting = true;
    clearTimeout(window.scrollTimeout);
    window.scrollTimeout = setTimeout(() => { isInteracting = false; }, 200);
    if(isExporting) return;
    e.preventDefault();
    const factor = Math.pow(1.00025, -e.deltaY);
    camZoom *= factor;
    camZoom = Math.max(0.01, camZoom); 
};

// ==========================================
// 5. EXPORT & WORKER INITIALIZATION (FINAL)
// ==========================================
const blob = new Blob([workerCode], { type: 'application/javascript' });
worker = new Worker(URL.createObjectURL(blob)); 

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
        
        gpuRenderedDensity = e.data.density || 1.0; 

        if (isExporting && exportResolve) {
            exportResolve();
            return;
        }

        // FIXED: Only update attempts when finding NEW parameters (mine/mutate)
        // This prevents the 'render' pass from overwriting the success message.
        if (e.data.source === 'mine' || e.data.source === 'mutate') {
            uiStatus.innerText = `FOUND! (${e.data.attempts} attempts)`;
            uiStatus.style.color = "#00ff00";
            currentCoeffs = new Float32Array(e.data.coeffs);
            
            if (e.data.source === 'mine') {
                colorSeed = [Math.random(), Math.random(), Math.random()];
                camPanX = 0; camPanY = 0; camZoom = 2.0; 
                currentQuat = qIdentity();
            }
            
            uiStatus.innerText = `Refining... (${e.data.attempts})`;
            worker.postMessage({ 
                type: 'render', 
                coeffs: currentCoeffs, 
                physicsSteps: currentPhysicsSteps, 
                density: currentDensity, 
                genType: currentGenType,
                // Pass constraints back if they exist, so the UI state remains consistent if we save
                constraints: e.data.constraints 
            });
        }
    }
};

let exportResolve = null;

async function startPrintCheckout(blob) {
  const uiExport = document.getElementById('ui-export-status');
  const container = document.getElementById('colorControls');
  
  if (!blob) {
      uiExport.innerText = "Error: No image data generated.";
      return;
  }

  if (typeof grecaptcha === 'undefined') {
      uiExport.innerText = "Error: reCAPTCHA not loaded. Refresh page.";
      return;
  }

  uiExport.innerText = "🤖 Verifying you are human...";
  const allButtons = document.querySelectorAll('button');
  allButtons.forEach(b => b.disabled = true);

  try {
    const token = await new Promise((resolve) => {
        grecaptcha.ready(() => {
            grecaptcha.execute(RECAPTCHA_SITE_KEY, {action: 'print_order'}).then(resolve);
        });
    });

    uiExport.innerText = "⏳ Requesting Cloud Storage...";

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

    const signedUrl = authData.publicUrl;
    
    uiExport.innerText = "✅ Ready.";

    const actionContainerId = 'pod-action-container';
    let actionContainer = document.getElementById(actionContainerId);
    if (actionContainer) actionContainer.remove(); 

    actionContainer = document.createElement('div');
    actionContainer.id = actionContainerId;
    actionContainer.style.marginTop = "10px";
    actionContainer.style.borderTop = "1px solid #555";
    actionContainer.style.paddingTop = "10px";
    
    // Match the side padding of the other UI elements
    actionContainer.style.margin = "0 10px 10px 10px"; 

    const inchesW = parseFloat(inpW.value);
    const inchesH = parseFloat(inpH.value);
    const dpi = parseInt(inpDPI.value);
    
    let totalW, totalH;
    
    if (exportUnit === 'pixels') {
        totalW = Math.floor(inchesW);
        totalH = Math.floor(inchesH);
    } else {
        totalW = Math.floor(inchesW * dpi);
        totalH = Math.floor(inchesH * dpi);
    }
    
    const masterCanvas = document.createElement('canvas');
    masterCanvas.width = totalW;
    masterCanvas.height = totalH;

    const peechoLink = document.createElement('a');
    peechoLink.href = "https://www.peecho.com"; 
    peechoLink.target = "_blank"; 

    peechoLink.className = "peecho-print-button";     
    peechoLink.setAttribute('data-src', signedUrl);
    peechoLink.setAttribute('data-thumbnail', signedUrl);
    peechoLink.setAttribute('data-filetype', 'image');
    peechoLink.setAttribute('data-width', masterCanvas.width);
    peechoLink.setAttribute('data-height', masterCanvas.height);
    peechoLink.setAttribute('data-pages', '1');
    peechoLink.setAttribute('data-reference', 'Strange Attractor #' + generateID());

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

    const footer = document.getElementById('ui-main-status')?.parentNode;
    if(footer) container.insertBefore(actionContainer, footer);
    else container.appendChild(actionContainer);

    const scriptId = 'peecho-sdk-script';
    if (!document.getElementById(scriptId)) {
        const script = document.createElement('script');
        script.id = scriptId;
        script.src = "https://d3aln0nj58oevo.cloudfront.net/button/script/177021676045966460.js";
        document.body.appendChild(script);
    } else {
        if (window.Peecho && window.Peecho.refresh) {
            console.log("Refreshing Peecho buttons...");
            window.Peecho.refresh();
        }
    }

    allButtons.forEach(b => b.disabled = false);

  } catch (err) {
    console.error(err);
    uiExport.innerText = "Error: " + err.message;
    allButtons.forEach(b => b.disabled = false);
  }
}

function renderTileParticles(totalW, totalH, tileBounds, opac, forcedAspect, jitter, overridePanX, overridePanY, overrideZoom) {
    //gl.activeTexture(gl.TEXTURE0);
    //gl.bindTexture(gl.TEXTURE_2D, gaussianTex); 
    //gl.uniform1i(gl.getUniformLocation(particleProgram, "u_sprite"), 0);

    const rotMatrix = qToMatrix(currentQuat);
    gl.uniformMatrix4fv(gl.getUniformLocation(particleProgram, "u_rotation"), false, new Float32Array(rotMatrix));
    
    // FIX: Allow overrides for PanX, PanY, and Zoom for accurate letterbox framing
    const usePanX = (overridePanX !== undefined) ? overridePanX : camPanX;
    const usePanY = (overridePanY !== undefined) ? overridePanY : camPanY;
    const useZoom = (overrideZoom !== undefined) ? overrideZoom : camZoom;
    
    gl.uniform2f(gl.getUniformLocation(particleProgram, "u_pan"), usePanX, usePanY);
    gl.uniform1f(gl.getUniformLocation(particleProgram, "u_zoom"), useZoom);
    
    gl.uniform1f(gl.getUniformLocation(particleProgram, "u_aspect"), forcedAspect);
    gl.uniform4f(gl.getUniformLocation(particleProgram, "u_tileBounds"), tileBounds[0], tileBounds[1], tileBounds[2], tileBounds[3]);
    
    gl.uniform1i(gl.getUniformLocation(particleProgram, "u_colorMode"), colorMode);
    gl.uniform3f(gl.getUniformLocation(particleProgram, "u_colorSeed"), colorSeed[0], colorSeed[1], colorSeed[2]);
    gl.uniform1i(gl.getUniformLocation(particleProgram, "u_invert"), 0); 
    gl.uniform1i(gl.getUniformLocation(particleProgram, "u_inc_black"), 1); 
    gl.uniform1i(gl.getUniformLocation(particleProgram, "u_inc_white"), 1); 
    gl.uniform1f(gl.getUniformLocation(particleProgram, "u_variation"), currentVariation); 

    const focusVal = (parseInt(document.getElementById('ui-focus').value) / 1000.0); 
    const focusSpanVal = (parseInt(document.getElementById('ui-focus-span').value) / 1000.0);
    const apertureVal = parseInt(document.getElementById('ui-aperture').value); 

    gl.uniform1f(gl.getUniformLocation(particleProgram, "u_focusDist"), focusVal);
    gl.uniform1f(gl.getUniformLocation(particleProgram, "u_focusSpan"), focusSpanVal);
    gl.uniform1f(gl.getUniformLocation(particleProgram, "u_aperture"), apertureVal);

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

async function startTiledExport(mode = 'download') {
    if (isExporting || !currentCoeffs) return;
    isExporting = true;
    
    const inchesW = parseFloat(inpW.value);
    const inchesH = parseFloat(inpH.value);
    const dpi = parseInt(inpDPI.value);
    const passes = parseInt(inpPasses.value);
    
    let totalW, totalH;
    
    if (exportUnit === 'pixels') {
        totalW = Math.floor(inchesW);
        totalH = Math.floor(inchesH);
    } else {
        totalW = Math.floor(inchesW * dpi);
        totalH = Math.floor(inchesH * dpi);
    }
    
    const masterCanvas = document.createElement('canvas');
    masterCanvas.width = totalW;
    masterCanvas.height = totalH;
    const ctx = masterCanvas.getContext('2d');
    
    const TILE_SIZE = 2048; 
    // Increased padding to handle large aperture/blur without clipping
    const PADDING = 512; 
    const cols = Math.ceil(totalW / TILE_SIZE);
    const rows = Math.ceil(totalH / TILE_SIZE);
    
    const fbo = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
    const tex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, tex);
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
    
    const exportOpacity = currentOpacity; 
    const exportJitter = currentJitter * resolutionRatio;

    const screenAspect = canvas.width / canvas.height;
    const printAspect = totalW / totalH;
    
    // FIX: Calculate precise framing adjustments for Letterbox vs Pillarbox
    let exportZoom = camZoom;
    let exportPanX = camPanX;
    let exportPanY = camPanY;

    if (printAspect > screenAspect) {
        // Letterbox mode: Print is relatively wider. Because the shader locks height, 
        // we must multiply the zoom and pan to physically scale the geometry up 
        // so the top and bottom get chopped off exactly as previewed.
        const framingRatio = printAspect / screenAspect;
        exportZoom *= framingRatio;
        exportPanX *= framingRatio;
        exportPanY *= framingRatio;
    }

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

            //gl.activeTexture(gl.TEXTURE0);
            //gl.bindTexture(gl.TEXTURE_2D, gaussianTex);
            //gl.uniform1i(gl.getUniformLocation(particleProgram, "u_sprite"), 0);
          
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
                        genType: currentGenType,
                        constraints: meta.constraints 
                    });
                });
                
                // FIX: Pass the compensated Zoom and Pan variables
                renderTileParticles(totalW, totalH, [nX, nY, nW, nH], exportOpacity, totalW/totalH, exportJitter, exportPanX, exportPanY, exportZoom);
            }
            
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
            gl.uniform1i(gl.getUniformLocation(compositeProgram, "u_transparent"), exportTransparent?1:0);
            gl.uniform1f(gl.getUniformLocation(compositeProgram, "u_scale"), 1.0); 
            
            let bMode = 0;
            if (blendMode === 'ADD') bMode = 1;
            gl.uniform1i(gl.getUniformLocation(compositeProgram, "u_blend_mode"), bMode);

            gl.drawArrays(gl.TRIANGLES, 0, 3);
            
            const pixels = new Uint8Array(tileW * tileH * 4);
            gl.bindFramebuffer(gl.READ_FRAMEBUFFER, resolveFbo);
            gl.readPixels(0, 0, tileW, tileH, gl.RGBA, gl.UNSIGNED_BYTE, pixels);
            
            const flipped = new Uint8ClampedArray(tileW * tileH * 4);
            const rowBytes = tileW * 4;
            for (let r = 0; r < tileH; r++) {
                const srcRow = pixels.subarray(r * rowBytes, (r + 1) * rowBytes);
                const dstOffset = (tileH - r - 1) * rowBytes;
                flipped.set(srcRow, dstOffset);
            }
            
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
    
    // FIX: Inject DPI metadata before Saving or Uploading
    masterCanvas.toBlob(async (rawBlob) => {
        // Inject the correct DPI (pHYs chunk)
        const finalBlob = await writePngDpi(rawBlob, dpi);

        if (mode === 'download') {
            uiExport.innerText = "Saving to disk...";
            const url = URL.createObjectURL(finalBlob); 
            const a = document.createElement('a');
            a.href = url; a.download = `attractor_${exportID}_${inchesW}x${inchesH}in_${dpi}dpi.png`; 
            document.body.appendChild(a); a.click(); document.body.removeChild(a);
            
            // Save JSON 
            const data = { coeffs: Array.from(currentCoeffs), settings: meta };
            const jsonBlob = new Blob([JSON.stringify(data, null, 2)], {type: "application/json"});
            const jUrl = URL.createObjectURL(jsonBlob);
            const jA = document.createElement('a');
            jA.href = jUrl; jA.download = `attractor_${exportID}.json`;
            document.body.appendChild(jA); jA.click(); document.body.removeChild(jA);
            
            resetRenderState();

        } else if (mode === 'pod') {
            uiExport.innerText = "Preparing Upload...";
            // Upload the blob with the correct metadata
            startPrintCheckout(finalBlob); 
        }
    }, 'image/png');

}

function resetRenderState() {
    isExporting = false;
    canvas.width = canvas.clientWidth; canvas.height = canvas.clientHeight;
    gl.viewport(0, 0, canvas.width, canvas.height);
    uiExport.innerText = "";
    worker.postMessage({ type: 'render', coeffs: currentCoeffs, physicsSteps: currentPhysicsSteps, density: currentDensity, genType: currentGenType });
}

function renderFrame() {
    const desktopScale = isInteracting ? 0.5 : 1.0;
    const mobileScale = 0.75; 
    let renderScaleTarget = isMobile ? mobileScale : desktopScale;
    
    if (Math.abs(renderScale - renderScaleTarget) > 0.01) {
        renderScale = renderScaleTarget;
        resizeViewportFBO();
    }
    
    if (!viewFbo) resizeViewportFBO();

    gl.bindFramebuffer(gl.FRAMEBUFFER, viewFbo);
    gl.viewport(0, 0, Math.floor(canvas.width * renderScale), Math.floor(canvas.height * renderScale));
    
    try {
        gl.clearColor(0,0,0,0);
        gl.clear(gl.COLOR_BUFFER_BIT);
        
        gl.useProgram(particleProgram);
        gl.enable(gl.BLEND);
        gl.blendFunc(gl.ONE, gl.ONE);
        
        //gl.activeTexture(gl.TEXTURE0);
        //gl.bindTexture(gl.TEXTURE_2D, gaussianTex); 
        //gl.uniform1i(gl.getUniformLocation(particleProgram, "u_sprite"), 0);
        
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

        const focusVal = (parseInt(document.getElementById('ui-focus').value) / 1000.0); 
        const focusSpanVal = (parseInt(document.getElementById('ui-focus-span').value) / 1000.0);
        const apertureVal = parseInt(document.getElementById('ui-aperture').value); 

        gl.uniform1f(gl.getUniformLocation(particleProgram, "u_focusDist"), focusVal);
        gl.uniform1f(gl.getUniformLocation(particleProgram, "u_focusSpan"), focusSpanVal);
        gl.uniform1f(gl.getUniformLocation(particleProgram, "u_aperture"), apertureVal);

        const interactionDimmer = isInteracting ? 0.4 : 1.0; 

        let targetOpacity = (currentOpacity / gpuRenderedDensity) * interactionDimmer; 
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

        gl.bindFramebuffer(gl.FRAMEBUFFER, null);
        gl.viewport(0, 0, canvas.width, canvas.height); 
        
        gl.useProgram(compositeProgram);
        gl.disable(gl.BLEND); 
        
        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, viewTex);
        gl.uniform1i(gl.getUniformLocation(compositeProgram, "u_tex"), 0);
        
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
        gl.uniform1i(gl.getUniformLocation(compositeProgram, "u_transparent"), 0); // Viewport always opaque

        let bMode = 0;
        if (blendMode === 'ADD') bMode = 1;
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
        // --- HIGH DPI FIX ---
        const dpr = window.devicePixelRatio || 1;
        const displayWidth  = Math.floor(canvas.clientWidth * dpr);
        const displayHeight = Math.floor(canvas.clientHeight * dpr);

        if (canvas.width !== displayWidth || canvas.height !== displayHeight) {
            canvas.width = displayWidth;
            canvas.height = displayHeight;
            resizeViewportFBO(); 
        }
        renderFrame();
    }
    requestAnimationFrame(loop);
}
requestAnimationFrame(loop);

// --- INITIAL AUTO-SEARCH ---
// Trigger the first search automatically so the canvas isn't empty on load.
setTimeout(() => {
    const uiStatus = document.getElementById('ui-main-status');
    if (uiStatus) {
        uiStatus.innerText = "Initializing First Attractor...";
        uiStatus.style.color = "#ffff00";
    }
    // Simulate a click on the hidden main mine button to kick off the worker
    if (btnMine) btnMine.click();
}, 100); // 100ms delay ensures the UI and WebGL context are fully painted and ready