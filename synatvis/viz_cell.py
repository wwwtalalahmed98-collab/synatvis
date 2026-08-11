"""Animated gene->protein->purified-protein journey (the integrative "living cell").

Self-contained HTML. Renders an original, labelled schematic Chlamydomonas cell
plus a downstream purification bench, and animates the construct through every
checkpoint built by journey.build_journey(). At each checkpoint a live readout
panel shows the actual parameters + scales the construct passes through, each with
a pass / warn / fail mark and a citation — so the journey is self-explanatory.

Original schematic (not microscopy); routes and values reflect the model's
predictions and published rules, not an observed movie of the protein.
"""
from __future__ import annotations

import json
from typing import Dict, List, Optional

_CELL_SVG = r"""
<svg id="scene" viewBox="0 0 980 770" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="Schematic Chlamydomonas cell and purification bench">
  <defs>
    <radialGradient id="gCyto" cx="50%" cy="40%" r="75%">
      <stop offset="0%" stop-color="var(--cyto0)"/><stop offset="100%" stop-color="var(--cyto1)"/></radialGradient>
    <radialGradient id="gChloro" cx="50%" cy="60%" r="70%">
      <stop offset="0%" stop-color="var(--chl0)"/><stop offset="100%" stop-color="var(--chl1)"/></radialGradient>
    <radialGradient id="gNuc" cx="45%" cy="40%" r="65%">
      <stop offset="0%" stop-color="var(--nuc0)"/><stop offset="100%" stop-color="var(--nuc1)"/></radialGradient>
    <radialGradient id="gPyr" cx="50%" cy="50%" r="60%">
      <stop offset="0%" stop-color="var(--pyr0)"/><stop offset="100%" stop-color="var(--pyr1)"/></radialGradient>
  </defs>

  <g stroke="var(--org)" stroke-width="3" fill="none" opacity="0.7" stroke-linecap="round">
    <path d="M430,110 C400,64 360,42 330,20"/><path d="M545,110 C575,64 615,42 645,20"/></g>

  <ellipse cx="490" cy="300" rx="452" ry="276" fill="var(--wall)" opacity="0.5"/>
  <ellipse id="membrane" cx="490" cy="300" rx="440" ry="266" fill="url(#gCyto)"
           stroke="var(--memb)" stroke-width="4"/>

  <path id="chloroplast" d="M112,278 C112,530 300,578 490,578 C680,578 868,530 868,278
      C868,445 700,512 490,512 C280,512 112,445 112,278 Z"
      fill="url(#gChloro)" stroke="var(--chlst)" stroke-width="3"/>
  <ellipse id="pyrenoid" cx="490" cy="484" rx="50" ry="42" fill="url(#gPyr)" stroke="var(--chlst)" stroke-width="2.5"/>
  <text x="490" y="489" class="orgtiny" text-anchor="middle">pyrenoid</text>
  <ellipse cx="176" cy="372" rx="16" ry="9" fill="var(--eye)" transform="rotate(-24 176 372)"/>

  <g opacity="0.95">
    <rect x="250" y="158" width="92" height="32" rx="16" fill="var(--mito)"/>
    <path d="M258,174 q7,-9 14,0 t14,0 t14,0 t14,0 t14,0" fill="none" stroke="#fff" stroke-width="1.5" opacity="0.45"/>
  </g>
  <g transform="rotate(18 790 203)" opacity="0.95">
    <rect x="742" y="188" width="96" height="30" rx="15" fill="var(--mito)"/>
    <path d="M750,203 q7,-8 14,0 t14,0 t14,0 t14,0 t14,0" fill="none" stroke="#fff" stroke-width="1.5" opacity="0.45"/>
  </g>
  <!-- chloroplast thylakoid stacks (grana) -->
  <g stroke="var(--chlst)" stroke-width="2.4" opacity="0.5" stroke-linecap="round">
    <line x1="300" y1="505" x2="342" y2="505"/><line x1="300" y1="512" x2="342" y2="512"/><line x1="300" y1="519" x2="342" y2="519"/>
    <line x1="640" y1="505" x2="682" y2="505"/><line x1="640" y1="512" x2="682" y2="512"/><line x1="640" y1="519" x2="682" y2="519"/>
  </g>

  <circle id="nucleus" cx="470" cy="238" r="86" fill="url(#gNuc)" stroke="var(--nucst)" stroke-width="3"/>
  <circle cx="452" cy="230" r="29" fill="var(--nucst)" opacity="0.4"/>
  <!-- nuclear-pore complexes around the envelope -->
  <g fill="var(--memb)" stroke="var(--nucst)" stroke-width="1.4" opacity="0.9">
    <circle cx="470" cy="153" r="4.5"/><circle cx="536" cy="184" r="4.5"/><circle cx="404" cy="184" r="4.5"/>
    <circle cx="391" cy="300" r="4.5"/><circle cx="470" cy="323" r="4.5"/></g>
  <circle id="pore" cx="554" cy="256" r="10" fill="var(--memb)" stroke="var(--nucst)" stroke-width="2"/>

  <path id="er" d="M554,290 C618,290 654,320 670,362 C684,398 674,442 638,460"
        fill="none" stroke="var(--er)" stroke-width="10" stroke-linecap="round" opacity="0.85"/>
  <g id="erRibos" fill="var(--org)"></g>

  <g id="golgi" transform="translate(710,382)">
    <path d="M-46,-20 C-10,-34 22,-34 52,-18" fill="none" stroke="var(--golgi)" stroke-width="8" stroke-linecap="round"/>
    <path d="M-50,-4 C-12,-18 24,-18 56,-2" fill="none" stroke="var(--golgi)" stroke-width="8" stroke-linecap="round"/>
    <path d="M-46,12 C-10,-2 22,-2 52,14" fill="none" stroke="var(--golgi)" stroke-width="8" stroke-linecap="round"/></g>

  <g class="lbl" fill="var(--ink)">
    <text x="470" y="144" text-anchor="middle">nucleus</text>
    <text x="586" y="252">nuclear pore</text>
    <text x="698" y="326">ER + ribosomes</text>
    <text x="710" y="430" text-anchor="middle">Golgi</text>
    <text x="490" y="566" text-anchor="middle">chloroplast</text>
    <text x="250" y="150">mitochondrion</text>
    <text x="150" y="358" text-anchor="end">eyespot</text>
    <text x="488" y="28" text-anchor="middle">flagella</text></g>

  <!-- purification bench -->
  <g id="bench">
    <line x1="60" y1="700" x2="920" y2="700" stroke="var(--benchline)" stroke-width="3"/>
    <text x="70" y="672" class="benchlbl">PURIFICATION BENCH</text>
    <g id="benchStations"></g>
  </g>

  <g id="mrnaTrack" opacity="0">
    <rect x="330" y="342" width="300" height="6" rx="3" fill="var(--mrna)"/>
    <circle id="ribosome" cx="330" cy="345" r="15" fill="var(--org)" opacity="0.95"/>
    <text x="480" y="380" text-anchor="middle" class="orgtiny">translation — ribosome dwell &prop; 1 / codon adaptiveness</text></g>

  <g id="foldMorph" opacity="0">
    <polyline id="foldChain" fill="none" stroke="var(--cargo)" stroke-width="3.4"
              stroke-linecap="round" stroke-linejoin="round"/>
    <text id="foldLabel" x="598" y="392" text-anchor="middle" class="orgtiny">folding — the chain collapses to its native shape</text></g>

  <g id="cargo">
    <circle id="cargoHalo" cx="470" cy="238" r="18" fill="none" stroke="var(--cargo)" stroke-width="2" opacity="0.5"/>
    <circle id="cargoDot" cx="470" cy="238" r="16" fill="var(--panel)" stroke="var(--cargo)" stroke-width="2.5"/>
    <g id="cargoGlyph"></g></g>
</svg>
"""

_TEMPLATE = r"""<div class="cellwrap">
  <div class="chead">
    <div><div class="ctitle">__NAME__</div>
      <div class="csub">Gene &rarr; protein &rarr; purified protein &middot; <span id="locName">-</span></div></div>
    <div class="cgauge"><div class="gnum" id="epiNum">--</div><div class="glab">expression<br>index</div></div>
  </div>

  <div class="lvstrip" id="lvstrip"></div>
  <div class="stagebar" id="stagebar"></div>
  __SVG__
  <div class="orglegend">
    <span><i style="background:var(--nuc1)"></i>Nucleus</span>
    <span><i style="background:var(--er)"></i>ER + ribosomes</span>
    <span><i style="background:var(--golgi)"></i>Golgi</span>
    <span><i style="background:var(--chl1)"></i>Chloroplast</span>
    <span><i style="background:var(--mito)"></i>Mitochondria</span>
    <span><i style="background:var(--benchline)"></i>Purification bench</span>
    <span><i style="background:var(--cargo)"></i>Your molecule</span>
  </div>

  <div class="readout" id="readout">
    <div class="rhead"><span id="rGlyph" class="rglyph"></span>
      <div><div id="rTitle" class="rtitle">Ready</div>
        <div id="rProg" class="rprog"></div></div></div>
    <div id="rParams" class="rparams"></div>
  </div>

  <div class="cfoot">
    <div class="ctrls">
      <button id="btnPlay">&#9654; Play</button>
      <button id="btnPrev">&#8592;</button><button id="btnNext">&#8594;</button>
      <button id="btnRestart">&#8635;</button>
      <label>speed <input id="spd" type="range" min="0.5" max="3" step="0.5" value="1"></label>
    </div>
  </div>
  <p class="disc">Original schematic (not microscopy). Every value is computed from the validated
  scan, the transparent expression ensemble, and published biophysics/purification rules — each
  parameter carries its citation. It is a model and a set of rules, not an observed movie of the protein.</p>
</div>

<style>
:root{
  --bg:#f6f8fb;--panel:#fff;--ink:#1c2733;--sub:#5a6b7b;--line:#e3e9f0;
  --cyto0:#eef6ff;--cyto1:#d7e8fb;--wall:#c7dcf0;--memb:#3f6f9f;
  --chl0:#8fd39a;--chl1:#3f9d63;--chlst:#2e7d4f;--pyr0:#eaf7ee;--pyr1:#bfe6c9;
  --nuc0:#e7dcfb;--nuc1:#c3a9ef;--nucst:#7a5bb0;--eye:#e8663a;
  --mito:#d98cc4;--er:#f0b36b;--golgi:#e0a34e;--org:#c05a2e;--mrna:#8a6bd6;
  --cargo:#1f7ae0;--glyco:#2e9d7c;--good:#1f9d55;--warn:#d98a1f;--bad:#cf4b3a;
  --benchline:#9fb0c2;--track:#eef1f5;}
@media(prefers-color-scheme:dark){:root{
  --bg:#0f151d;--panel:#161f2b;--ink:#e7eef6;--sub:#9fb0c2;--line:#26323f;
  --cyto0:#16283b;--cyto1:#102033;--wall:#1b3450;--memb:#5b90c4;
  --chl0:#2f7a4c;--chl1:#175f38;--chlst:#4fbf7f;--pyr0:#1d3b2a;--pyr1:#255a3c;
  --nuc0:#39295e;--nuc1:#5a3f8f;--nucst:#b79be6;--eye:#ff7a4a;
  --mito:#b95fa3;--er:#d9954a;--golgi:#c98a37;--org:#e07a4a;--mrna:#a487e6;
  --cargo:#4aa0ff;--glyco:#43c79b;--good:#37c07a;--warn:#e6a53a;--bad:#e5674f;
  --benchline:#3a4a5c;--track:#1c2733;}}
:root[data-theme="dark"]{--bg:#0f151d;--panel:#161f2b;--ink:#e7eef6;--sub:#9fb0c2;--line:#26323f;--cyto0:#16283b;--cyto1:#102033;--wall:#1b3450;--memb:#5b90c4;--chl0:#2f7a4c;--chl1:#175f38;--chlst:#4fbf7f;--pyr0:#1d3b2a;--pyr1:#255a3c;--nuc0:#39295e;--nuc1:#5a3f8f;--nucst:#b79be6;--eye:#ff7a4a;--mito:#b95fa3;--er:#d9954a;--golgi:#c98a37;--org:#e07a4a;--mrna:#a487e6;--cargo:#4aa0ff;--glyco:#43c79b;--good:#37c07a;--warn:#e6a53a;--bad:#e5674f;--benchline:#3a4a5c;--track:#1c2733;}
:root[data-theme="light"]{--bg:#f6f8fb;--panel:#fff;--ink:#1c2733;--sub:#5a6b7b;--line:#e3e9f0;--cyto0:#eef6ff;--cyto1:#d7e8fb;--wall:#c7dcf0;--memb:#3f6f9f;--chl0:#8fd39a;--chl1:#3f9d63;--chlst:#2e7d4f;--pyr0:#eaf7ee;--pyr1:#bfe6c9;--nuc0:#e7dcfb;--nuc1:#c3a9ef;--nucst:#7a5bb0;--eye:#e8663a;--mito:#d98cc4;--er:#f0b36b;--golgi:#e0a34e;--org:#c05a2e;--mrna:#8a6bd6;--cargo:#1f7ae0;--glyco:#2e9d7c;--good:#1f9d55;--warn:#d98a1f;--bad:#cf4b3a;--benchline:#9fb0c2;--track:#eef1f5;}
*{box-sizing:border-box}
.cellwrap{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;color:var(--ink);
  background:var(--panel);border:1px solid var(--line);border-radius:16px;padding:18px;
  max-width:1060px;margin:0 auto;box-shadow:0 8px 30px rgba(20,40,70,.08)}
.chead{display:flex;justify-content:space-between;align-items:center;gap:12px;margin-bottom:8px}
.ctitle{font-size:20px;font-weight:700}.csub{color:var(--sub);font-size:13px;margin-top:2px}
.cgauge{display:flex;align-items:center;gap:8px;background:var(--bg);border:1px solid var(--line);border-radius:12px;padding:6px 12px}
.gnum{font-size:26px;font-weight:800;line-height:1}
.glab{font-size:10px;color:var(--sub);text-transform:uppercase;letter-spacing:.4px}
.lvstrip{display:grid;grid-template-columns:repeat(6,1fr);gap:5px;margin:8px 0 6px}
.lvcell{border:1px solid var(--line);border-top:3px solid var(--lc);border-radius:8px;
  padding:5px 7px;background:var(--bg)}
.lvcell .lvh{display:flex;align-items:center;gap:5px}
.lvcell .lvn{font-size:9px;font-weight:700;color:var(--sub)}
.lvcell .lvi{font-size:13px;font-weight:800;margin-left:auto;color:var(--lc)}
.lvcell .lvt{font-size:10px;font-weight:700;margin-top:2px;color:var(--ink);line-height:1.2}
.stagebar{display:flex;gap:5px;margin:6px 0 8px;flex-wrap:wrap}
.stagebar .st{border-left:3px solid var(--chipc,var(--line))}
.stagebar .st{font-size:10.5px;padding:5px 7px;border-radius:7px;background:var(--bg);
  border:1px solid var(--line);color:var(--sub);transition:.2s;cursor:pointer;white-space:nowrap}
.stagebar .st.active{background:var(--cargo);color:#fff;border-color:var(--cargo);font-weight:700}
.stagebar .st.done{color:var(--good);border-color:var(--good)}
#scene{width:100%;height:auto;display:block}
.lbl text{font-size:12px;font-weight:600}.orgtiny{font-size:11px;fill:var(--sub)}
.benchlbl{font-size:12px;font-weight:700;fill:var(--sub);letter-spacing:.5px}
.orglegend{display:flex;flex-wrap:wrap;gap:8px 14px;font-size:11px;color:var(--sub);margin:2px 2px 6px}
.orglegend i{display:inline-block;width:10px;height:10px;border-radius:3px;margin-right:5px;vertical-align:middle}
.readout{background:var(--bg);border:1px solid var(--line);border-radius:12px;padding:12px 14px;margin-top:6px;min-height:120px}
.rhead{display:flex;align-items:center;gap:12px;margin-bottom:8px}
.rglyph{width:44px;height:44px;flex:0 0 44px;display:flex;align-items:center;justify-content:center;
  background:var(--panel);border:1px solid var(--line);border-radius:10px}
.rglyph svg{width:30px;height:30px}
.rtitle{font-size:15px;font-weight:700}.rprog{font-size:11px;color:var(--sub)}
.rparams{display:flex;flex-direction:column;gap:7px}
.prow{display:grid;grid-template-columns:22px minmax(140px,1.4fr) auto 1fr;gap:8px;align-items:center}
.pico{width:18px;height:18px;border-radius:50%;display:inline-flex;align-items:center;justify-content:center;
  font-size:11px;font-weight:800;color:#fff}
.pico.ok{background:var(--good)}.pico.warn{background:var(--warn)}.pico.bad{background:var(--bad)}
.pico.info{background:var(--sub)}
.plabel{font-size:12.5px;color:var(--ink)}.pval{font-size:12.5px;font-weight:650;text-align:right;font-variant-numeric:tabular-nums}
.pscale{height:7px;background:var(--track);border:1px solid var(--line);border-radius:999px;overflow:hidden}
.pfill{height:100%;border-radius:999px}
.pdetail{grid-column:2 / 5;font-size:11px;color:var(--sub);line-height:1.45;margin:-1px 0 3px}
.pref{grid-column:2 / 5;font-size:10px;color:var(--sub);opacity:.8;font-style:italic}
.cfoot{margin-top:10px}
.ctrls{display:flex;align-items:center;gap:8px;font-size:12px;color:var(--sub);flex-wrap:wrap}
.ctrls button{font:inherit;font-size:13px;border:1px solid var(--line);background:var(--bg);color:var(--ink);
  border-radius:9px;padding:6px 11px;cursor:pointer}.ctrls button:hover{border-color:var(--cargo)}
.disc{font-size:11px;color:var(--sub);margin:10px 2px 0;line-height:1.5}
#cargoHalo{animation:pulse 1.7s ease-in-out infinite}
@keyframes pulse{0%,100%{r:16;opacity:.5}50%{r:27;opacity:0}}
@media(prefers-reduced-motion:reduce){#cargoHalo{animation:none}}
@media(max-width:640px){.prow{grid-template-columns:20px 1fr auto;}.pscale{display:none}}
</style>

<script>
/*__DATA__*/
(function(){
  var svg=document.getElementById('scene'); if(!svg) return;
  var NS='http://www.w3.org/2000/svg';
  var s='var(--ink)', a='var(--cargo)';   // glyph colours (used early by the bench)
  var dot=document.getElementById('cargoDot'), halo=document.getElementById('cargoHalo');
  var glyphG=document.getElementById('cargoGlyph');
  // the travelling cargo shows the molecule at each stage; map organelle-ish symbols to a molecule form
  var MOLMAP={pore:'rna',er:'protein',bench:'protein',spliceosome:'rna',thermo:'protein',
    scissors:'protein',bead:'protein',lysis:'protein',vesicle:'protein',glycan:'glycan',polya:'rna'};
  var track=document.getElementById('mrnaTrack'), ribo=document.getElementById('ribosome');
  document.getElementById('epiNum').textContent=Math.round(DATA.epi);
  var eEl=document.getElementById('epiNum');
  eEl.style.color=DATA.epi>=65?'var(--good)':DATA.epi>=40?'var(--warn)':'var(--bad)';
  document.getElementById('locName').textContent=DATA.loc;

  // ER-bound ribosome dots
  var er=document.getElementById('er'),erR=document.getElementById('erRibos'),L=er.getTotalLength();
  for(var t=0.12;t<0.95;t+=0.16){var pt=er.getPointAtLength(L*t);
    var c=document.createElementNS(NS,'circle');c.setAttribute('cx',pt.x);c.setAttribute('cy',pt.y);
    c.setAttribute('r',4.5);c.setAttribute('fill','var(--org)');erR.appendChild(c);}

  // organelle coordinates
  var ORG={nucleus:[470,238],nucleusEdge:[512,290],pore:[554,256],ribo:[480,345],
    fold:[598,342],er:[636,420],golgi:[710,382],secreted:[946,300],membrane:[900,290],
    cyto:[330,300]};

  // lay out bench stations for 'bench' checkpoints
  var benchCps=DATA.checkpoints.filter(function(c){return c.organelle==='bench';});
  var bx0=140,bx1=860,by=700, bStations=document.getElementById('benchStations');
  benchCps.forEach(function(c,i){
    var x=benchCps.length>1?bx0+(bx1-bx0)*i/(benchCps.length-1):(bx0+bx1)/2;
    c._xy=[x,by-4];
    var g=document.createElementNS(NS,'g'); g.setAttribute('transform','translate('+x+','+(by)+')');
    g.setAttribute('id','bs'+i); g.setAttribute('opacity','0.4');
    var circle=document.createElementNS(NS,'circle'); circle.setAttribute('r',18);
    circle.setAttribute('cy',-18); circle.setAttribute('fill','var(--panel)');
    circle.setAttribute('stroke','var(--benchline)'); circle.setAttribute('stroke-width','2'); g.appendChild(circle);
    var gl=document.createElementNS(NS,'g'); gl.setAttribute('transform','translate(-11,-29) scale(0.73)');
    gl.innerHTML=glyphInner(c.symbol); g.appendChild(gl);
    var tx=document.createElementNS(NS,'text'); tx.setAttribute('y',18); tx.setAttribute('text-anchor','middle');
    tx.setAttribute('class','orgtiny'); tx.textContent=c.title.split(/[ (]/)[0]; g.appendChild(tx);
    bStations.appendChild(g);
  });
  function xy(cp){ return cp._xy || ORG[cp.organelle] || [490,300]; }

  // ---- scientific glyph library (inner markup drawn in a 0 0 30 30 box) ----
  function glyphInner(k){
    var G={
     dna:'<path d="M9 3 C21 9 9 21 21 27 M21 3 C9 9 21 21 9 27" stroke="'+a+'" stroke-width="2" fill="none"/><path d="M11 7h8M11 12h8M11 18h8M11 23h8" stroke="'+s+'" stroke-width="1.4"/>',
     rna:'<path d="M6 20 C12 6 18 26 26 10" stroke="'+a+'" stroke-width="2.4" fill="none"/><circle cx="6" cy="20" r="2.4" fill="'+s+'"/>',
     spliceosome:'<circle cx="15" cy="15" r="8" fill="none" stroke="'+a+'" stroke-width="2"/><path d="M15 7 C24 9 24 21 15 23" stroke="'+s+'" stroke-width="1.6" fill="none"/>',
     polya:'<path d="M4 15h12" stroke="'+s+'" stroke-width="2"/><text x="16" y="19" font-size="9" font-weight="700" fill="'+a+'">AAAA</text>',
     pore:'<circle cx="15" cy="15" r="9" fill="none" stroke="'+a+'" stroke-width="2"/><circle cx="15" cy="15" r="3.4" fill="'+s+'"/>',
     ribosome:'<ellipse cx="15" cy="12" rx="9" ry="6" fill="'+a+'" opacity="0.85"/><ellipse cx="15" cy="20" rx="7" ry="5" fill="'+s+'" opacity="0.7"/>',
     fold:'<path d="M5 15 C10 5 12 25 16 15 C20 5 22 25 26 15" stroke="'+a+'" stroke-width="2.4" fill="none"/>',
     er:'<path d="M5 9 C14 4 18 14 26 9 M5 15 C14 10 18 20 26 15 M5 21 C14 16 18 26 26 21" stroke="'+a+'" stroke-width="1.8" fill="none"/>',
     glycan:'<path d="M15 26 V16 M15 16 L9 9 M15 16 L21 9" stroke="'+s+'" stroke-width="1.8" fill="none"/><rect x="12" y="24" width="6" height="4" fill="'+s+'"/><circle cx="9" cy="8" r="3" fill="'+a+'"/><circle cx="21" cy="8" r="3" fill="'+a+'"/><circle cx="15" cy="16" r="3" fill="'+a+'"/>',
     vesicle:'<circle cx="15" cy="15" r="9" fill="none" stroke="'+a+'" stroke-width="2"/><circle cx="15" cy="15" r="3.6" fill="'+s+'"/>',
     protein:'<path d="M9 12 C7 20 16 26 21 19 C26 12 18 6 13 9 C10 10 10 11 9 12Z" fill="'+a+'" opacity="0.85"/>',
     lysis:'<path d="M15 4 L18 12 L26 12 L20 17 L23 26 L15 21 L7 26 L10 17 L4 12 L12 12 Z" fill="none" stroke="'+a+'" stroke-width="1.6"/>',
     thermo:'<rect x="12" y="4" width="6" height="16" rx="3" fill="none" stroke="'+s+'" stroke-width="1.6"/><circle cx="15" cy="23" r="5" fill="'+a+'"/><rect x="13.5" y="10" width="3" height="12" fill="'+a+'"/>',
     bead:'<circle cx="15" cy="15" r="9" fill="'+a+'" opacity="0.3" stroke="'+a+'" stroke-width="1.6"/><path d="M9 15h12M15 9v12" stroke="'+s+'" stroke-width="1.2"/>',
     scissors:'<circle cx="9" cy="9" r="3.4" fill="none" stroke="'+a+'" stroke-width="1.8"/><circle cx="9" cy="21" r="3.4" fill="none" stroke="'+a+'" stroke-width="1.8"/><path d="M11.5 10.5 L25 20 M11.5 19.5 L25 10" stroke="'+s+'" stroke-width="1.8"/>',
     vial:'<path d="M12 4 h6 v7 l4 12 a3 3 0 0 1 -3 4 h-8 a3 3 0 0 1 -3 -4 l4 -12 Z" fill="none" stroke="'+s+'" stroke-width="1.6"/><path d="M11 18 h8 l2 5 a3 3 0 0 1 -3 4 h-6 a3 3 0 0 1 -3 -4 Z" fill="'+a+'" opacity="0.75"/>',
     bench:'<rect x="5" y="18" width="20" height="3" fill="'+s+'"/><rect x="8" y="8" width="4" height="10" fill="'+a+'"/><rect x="18" y="6" width="4" height="12" fill="'+a+'"/>',
     complex:'<circle cx="11" cy="12" r="6" fill="'+a+'" opacity="0.75"/><circle cx="20" cy="11" r="5.4" fill="'+s+'" opacity="0.75"/><circle cx="15" cy="20" r="5.6" fill="'+a+'" opacity="0.55"/>'
    };
    return (G[k]||G.protein);
  }
  function glyph(k){return '<svg viewBox="0 0 30 30" xmlns="http://www.w3.org/2000/svg">'+glyphInner(k)+'</svg>';}
  var ICO={ok:'✓',warn:'!',bad:'✗',info:'i'};

  // ---- six molecular-biology levels ----
  var LV=[['pretranscript','Pre-transcription'],['transcript','Transcription'],
    ['posttranscript','Post-transcription'],['pretranslation','Pre-translation'],
    ['translation','Translation'],['posttranslation','Post-translation']];
  var LVCOL={pretranscript:'#7a5bb0',transcript:'#256abf',posttranscript:'#2e8d6b',
    pretranslation:'#c98a37',translation:'#c05a2e',posttranslation:'#b23b8f'};
  var RANK={info:0,ok:1,warn:2,bad:3}, STMARK={0:['✓','var(--good)'],1:['✓','var(--good)'],
    2:['!','var(--warn)'],3:['✗','var(--bad)']};
  function levelStatus(key){var w=0; DATA.checkpoints.forEach(function(c){
    if(c.level===key)(c.params||[]).forEach(function(p){w=Math.max(w,RANK[p.status||'info']);});}); return w;}
  var strip=document.getElementById('lvstrip'); strip.innerHTML='';
  LV.forEach(function(pair,i){var key=pair[0], w=levelStatus(key), m=STMARK[w];
    var d=document.createElement('div'); d.className='lvcell'; d.style.setProperty('--lc',LVCOL[key]);
    d.innerHTML='<div class="lvh"><span class="lvn">'+(i+1)+'</span>'
      +'<span class="lvi" style="color:'+m[1]+'">'+m[0]+'</span></div>'
      +'<div class="lvt">'+pair[1]+'</div>';
    strip.appendChild(d);});

  // stage rail (chips colour-coded by level)
  var bar=document.getElementById('stagebar'); bar.innerHTML='';
  var chips=DATA.checkpoints.map(function(c,i){var d=document.createElement('div');
    d.className='st'; d.textContent=(i+1)+'. '+c.title.split(/[ (]/)[0];
    d.style.borderLeft='3px solid '+(LVCOL[c.level]||'#9aa6b2');
    d.onclick=function(){goto(i);}; bar.appendChild(d); return d;});

  // ---- render one checkpoint's parameter panel ----
  function scaleColor(st){return st==='bad'?'var(--bad)':st==='warn'?'var(--warn)':'var(--good)';}
  function showParams(c,idx){
    document.getElementById('rGlyph').innerHTML=glyph(c.symbol);
    glyphG.innerHTML=glyphInner(MOLMAP[c.symbol]||c.symbol);   // morph the travelling molecule
    document.getElementById('rTitle').textContent=c.title;
    document.getElementById('rProg').textContent='Checkpoint '+(idx+1)+' of '+DATA.checkpoints.length
      +(c.landscape?'  ·  translation landscape shown in-cell':'');
    var box=document.getElementById('rParams'); box.innerHTML='';
    c.params.forEach(function(p){
      var row=document.createElement('div'); row.className='prow';
      var st=p.status||'info';
      var ico='<span class="pico '+st+'">'+ICO[st]+'</span>';
      var scale=''; if(p.scale!==null && p.scale!==undefined){
        scale='<span class="pscale"><span class="pfill" style="width:'+Math.max(3,Math.min(100,p.scale))
              +'%;background:'+scaleColor(st)+'"></span></span>';
      } else { scale='<span></span>'; }
      row.innerHTML=ico+'<span class="plabel">'+esc(p.label)+'</span>'
        +'<span class="pval">'+esc(String(p.value))+'</span>'+scale;
      box.appendChild(row);
      if(p.detail){var d=document.createElement('div'); d.className='pdetail'; d.textContent=p.detail; box.appendChild(d);}
      if(p.ref){var r=document.createElement('div'); r.className='pref'; r.textContent='ref: '+p.ref; box.appendChild(r);}
    });
  }
  function esc(s){return s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');}

  // ---- animation ----
  var land=(DATA.landscape&&DATA.landscape.length)?DATA.landscape:[1];
  var mx=Math.max.apply(null,land)||1; land=land.map(function(v){return v/mx;});
  var i=0, playing=false, raf=null, speed=1;
  document.getElementById('spd').oninput=function(e){speed=parseFloat(e.target.value);};
  function setDot(x,y){dot.setAttribute('cx',x);dot.setAttribute('cy',y);halo.setAttribute('cx',x);halo.setAttribute('cy',y);
    glyphG.setAttribute('transform','translate('+(x-15)+','+(y-15)+')');}
  function ease(t){return t<.5?4*t*t*t:1-Math.pow(-2*t+2,3)/2;}   // easeInOutCubic (smoother)

  // ---- protein-folding morph: extended chain -> compact native fold ----
  var foldMorph=document.getElementById('foldMorph'), foldChain=document.getElementById('foldChain'), cargoG=document.getElementById('cargo');
  var _fx=ORG.fold[0], _fy=ORG.fold[1], NPT=14;
  var EXT=[], FOLD=[];
  for(var _k=0;_k<NPT;_k++){ EXT.push([_fx-70+140*_k/(NPT-1), _fy+((_k%2)?13:-13)]); }
  var _ang=[0,1.1,2.3,3.4,4.6,5.7,0.6,1.8,2.9,4.1,5.2,1.2,2.4,3.6], _rr=[24,19,25,16,23,18,13,10,15,9,12,6,5,7];
  for(var _j=0;_j<NPT;_j++){ FOLD.push([_fx+Math.cos(_ang[_j])*_rr[_j], _fy+Math.sin(_ang[_j])*_rr[_j]]); }
  function drawFold(e){
    var pts=EXT.map(function(pt,k){return (pt[0]+(FOLD[k][0]-pt[0])*e).toFixed(1)+','+(pt[1]+(FOLD[k][1]-pt[1])*e).toFixed(1);});
    foldChain.setAttribute('points', pts.join(' '));
  }
  function hideFold(){foldMorph.setAttribute('opacity','0'); cargoG.setAttribute('opacity','1');}

  function lightBench(c){
    benchCps.forEach(function(bc,bi){var g=document.getElementById('bs'+bi);
      if(g) g.setAttribute('opacity', bc===c?'1':(bc._done?'0.85':'0.4'));});
  }

  function animateTo(idx,cb){
    var c=DATA.checkpoints[idx], from=i>=0&&DATA.checkpoints[i]?xy(DATA.checkpoints[i]):xy(c), to=xy(c);
    chips.forEach(function(ch,j){ch.classList.toggle('done',j<idx);ch.classList.toggle('active',j===idx);});
    showParams(c,idx);
    if(c.organelle==='bench'){lightBench(c);}
    var xl=c.landscape; if(xl) track.setAttribute('opacity','1');
    var fold=(c.key==='folding');
    if(fold){foldMorph.setAttribute('opacity','1'); cargoG.setAttribute('opacity','0.15'); drawFold(0);} else {hideFold();}
    var start=performance.now(), dur=(fold?1500:1100)/speed;
    function step(now){var t=Math.min(1,(now-start)/dur),e=ease(t);
      setDot(from[0]+(to[0]-from[0])*e, from[1]+(to[1]-from[1])*e);
      if(xl){var n=land.length, idx2=Math.min(n-1,Math.floor(t*n)), ad=land[idx2]||0.5;
        ribo.setAttribute('cx',330+300*t); ribo.setAttribute('r',15-6*Math.min(1,ad));
        ribo.setAttribute('fill',ad<0.35?'var(--bad)':ad<0.6?'var(--warn)':'var(--org)');}
      if(fold){ drawFold(ease(Math.min(1,t*1.15))); }   // chain collapses to native fold
      if(t<1){raf=requestAnimationFrame(step);}
      else{ if(xl){track.setAttribute('opacity','0');} if(c.organelle==='bench'){c._done=true;}
        i=idx; if(cb)cb(); }
    }
    raf=requestAnimationFrame(step);
  }
  function playFrom(idx){ if(idx>=DATA.checkpoints.length){finish();return;}
    animateTo(idx,function(){ setTimeout(function(){ if(playing) playFrom(idx+1); }, 650/speed); }); }
  function finish(){playing=false;document.getElementById('btnPlay').innerHTML='&#9654; Replay';}
  function goto(idx){cancelAnimationFrame(raf);playing=false;
    document.getElementById('btnPlay').innerHTML='&#9654; Play';
    benchCps.forEach(function(bc){bc._done=false;});
    // jump: set state up to idx
    var prev=idx>0?idx-1:idx; i=prev; animateTo(idx);}
  function reset(){cancelAnimationFrame(raf);i=0;playing=false;
    document.getElementById('btnPlay').innerHTML='&#9654; Play';
    benchCps.forEach(function(bc){bc._done=false;});
    setDot(ORG.nucleus[0],ORG.nucleus[1]); track.setAttribute('opacity','0'); hideFold();
    chips.forEach(function(ch){ch.className='st';}); showParams(DATA.checkpoints[0],0);
    lightBench(null);}

  document.getElementById('btnPlay').onclick=function(){
    if(i>=DATA.checkpoints.length-1 && !playing){reset();}
    playing=!playing; this.innerHTML=playing?'&#10073;&#10073; Pause':'&#9654; Play';
    if(playing) playFrom(Math.max(0,i));
  };
  document.getElementById('btnNext').onclick=function(){goto(Math.min(DATA.checkpoints.length-1,i+1));};
  document.getElementById('btnPrev').onclick=function(){goto(Math.max(0,i-1));};
  document.getElementById('btnRestart').onclick=reset;
  reset();

  // report our content height to a hosting page so an embedding iframe can self-size.
  // Harmless when not embedded (posts to self; no listener).
  var lastH=0;
  function reportHeight(){
    var h=Math.max(document.body.scrollHeight, document.documentElement.scrollHeight);
    if(Math.abs(h-lastH)<2) return; lastH=h;
    try{ parent.postMessage({synatvisHeight:h}, '*'); }catch(e){}
  }
  window.addEventListener('load', reportHeight);
  if(window.ResizeObserver){ try{ new ResizeObserver(reportHeight).observe(document.body); }catch(e){} }
  setTimeout(reportHeight, 150);
})();
</script>
"""


def render_cell_html(journey: Dict) -> str:
    safe = (journey.get("name") or "construct").replace("<", "&lt;").replace(">", "&gt;")
    html = _TEMPLATE.replace("__NAME__", safe).replace("__SVG__", _CELL_SVG)
    return html.replace("/*__DATA__*/", "var DATA = " + json.dumps(journey) + ";")


def render_cell_document(journey: Dict) -> str:
    body = render_cell_html(journey)
    name = journey.get("name") or "construct"
    return ("<!doctype html><html><head><meta charset='utf-8'>"
            "<meta name='viewport' content='width=device-width,initial-scale=1'>"
            f"<title>SynAT.Vis - journey: {name}</title>"
            "<style>body{margin:0;padding:20px;background:var(--bg,#f6f8fb)}</style>"
            "</head><body>" + body + "</body></html>")
