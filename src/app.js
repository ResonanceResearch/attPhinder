/* attP / attB Finder
   Static browser app. No sequence leaves the browser.
*/
const state = {
  phageFiles: [],
  bactFiles: [],
  phages: [],
  bacteria: [],
  results: []
};

const $ = (id) => document.getElementById(id);
const logBox = $('log');

function log(msg) {
  const stamp = new Date().toLocaleTimeString();
  logBox.textContent += `\n[${stamp}] ${msg}`;
  logBox.scrollTop = logBox.scrollHeight;
}
function resetLog(msg = 'Waiting for files…') { logBox.textContent = msg; }

function rc(seq) {
  const map = {A:'T', T:'A', G:'C', C:'G', N:'N', a:'t', t:'a', g:'c', c:'g', n:'n'};
  return seq.split('').reverse().map(b => map[b] || 'N').join('').toUpperCase();
}
function cleanSeq(s) { return (s || '').replace(/[^A-Za-z]/g, '').toUpperCase().replace(/U/g, 'T').replace(/[^ACGTN]/g, 'N'); }
function clamp(n, lo, hi) { return Math.max(lo, Math.min(hi, n)); }
function fmt(n) { return Number.isFinite(n) ? n.toLocaleString() : ''; }
function safeName(name) { return name.replace(/[^a-z0-9._-]+/gi, '_'); }

function parseFasta(text, fileName) {
  const records = [];
  let header = fileName, seq = '';
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith('>')) {
      if (seq) records.push({ name: header, sequence: cleanSeq(seq), features: [], format: 'FASTA' });
      header = line.slice(1).trim() || fileName;
      seq = '';
    } else {
      seq += line.trim();
    }
  }
  if (seq) records.push({ name: header, sequence: cleanSeq(seq), features: [], format: 'FASTA' });
  if (!records.length) records.push({ name: fileName, sequence: cleanSeq(text), features: [], format: 'FASTA' });
  return records;
}

function splitGenBankRecords(text) {
  const parts = text.split(/^\/\/\s*$/m).map(x => x.trim()).filter(Boolean);
  return parts.length ? parts : [text];
}

function parseGenBank(text, fileName) {
  const records = [];
  for (const rec of splitGenBankRecords(text)) {
    const locusMatch = rec.match(/^LOCUS\s+(\S+)/m);
    const defMatch = rec.match(/^DEFINITION\s+([\s\S]*?)(?=^ACCESSION\s|^VERSION\s|^KEYWORDS\s|^SOURCE\s)/m);
    const name = locusMatch?.[1] || defMatch?.[1]?.replace(/\s+/g, ' ').trim() || fileName;
    const origin = rec.match(/^ORIGIN([\s\S]*)/m);
    const sequence = origin ? cleanSeq(origin[1]) : '';
    const featuresBlock = rec.match(/^FEATURES\s+Location\/Qualifiers([\s\S]*?)(?=^ORIGIN|^CONTIG|^BASE COUNT|$)/m)?.[1] || '';
    const features = parseFeatures(featuresBlock);
    records.push({ name, fileName, sequence, features, format: 'GenBank' });
  }
  return records.filter(r => r.sequence.length > 0);
}

function parseFeatures(block) {
  const lines = block.split(/\r?\n/);
  const features = [];
  let current = null;
  for (const line of lines) {
    const start = line.match(/^\s{5}(\S+)\s+(.+)/);
    if (start) {
      if (current) features.push(finalizeFeature(current));
      current = { type: start[1], locationRaw: start[2].trim(), qualifiersRaw: [] };
    } else if (current) {
      const cont = line.slice(21).trim();
      if (cont.startsWith('/')) current.qualifiersRaw.push(cont);
      else if (current.qualifiersRaw.length && !cont.match(/^[a-zA-Z_]+\(/) && !cont.match(/^[<>]?\d/)) current.qualifiersRaw[current.qualifiersRaw.length - 1] += ' ' + cont;
      else current.locationRaw += cont;
    }
  }
  if (current) features.push(finalizeFeature(current));
  return features;
}

function finalizeFeature(f) {
  f.qualifiers = {};
  for (const q of f.qualifiersRaw) {
    const m = q.match(/^\/(\S+?)(?:=(.*))?$/);
    if (!m) continue;
    let val = m[2] || true;
    if (typeof val === 'string') val = val.replace(/^"|"$/g, '').replace(/"\s+"/g, ' ');
    if (f.qualifiers[m[1]]) {
      if (!Array.isArray(f.qualifiers[m[1]])) f.qualifiers[m[1]] = [f.qualifiers[m[1]]];
      f.qualifiers[m[1]].push(val);
    } else f.qualifiers[m[1]] = val;
  }
  const loc = parseLocation(f.locationRaw);
  return { ...f, ...loc };
}

function parseLocation(raw) {
  let s = raw.replace(/\s+/g, '');
  let strand = '+';
  if (s.startsWith('complement(') && s.endsWith(')')) {
    strand = '-';
    s = s.slice(11, -1);
  }
  if (s.startsWith('join(') && s.endsWith(')')) s = s.slice(5, -1);
  s = s.replace(/complement\(|join\(|\)/g, '');
  const nums = [...s.matchAll(/<?(\d+)\.\.>?(\d+)|<?(\d+)/g)].map(m => [Number(m[1] || m[3]), Number(m[2] || m[3])]);
  if (!nums.length) return { start: null, end: null, strand, parts: [] };
  const start = Math.min(...nums.map(p => Math.min(p[0], p[1])));
  const end = Math.max(...nums.map(p => Math.max(p[0], p[1])));
  return { start, end, strand, parts: nums };
}

function qualText(feature) {
  const q = feature.qualifiers || {};
  return ['gene', 'product', 'note', 'function', 'locus_tag'].map(k => Array.isArray(q[k]) ? q[k].join(' ') : (q[k] || '')).join(' ');
}

function findIntegrases(record) {
  const candidates = record.features.filter(f => f.type === 'CDS' && f.start && f.end).map(f => {
    const t = qualText(f).toLowerCase();
    let score = 0;
    if (/\bintegrase\b|\bint\b/.test(t)) score += 60;
    if (/recombinase|tyrosine|serine|site-specific/.test(t)) score += 25;
    if (/transposase|terminase|portal|capsid|tail|hypothetical/.test(t)) score -= /hypothetical/.test(t) ? 5 : 25;
    return { feature: f, text: qualText(f), score };
  }).filter(c => c.score > 0).sort((a,b) => b.score - a.score);

  // If no annotated CDS is found, return a soft placeholder so users still get feedback.
  return candidates;
}

function extractWindow(seq, start1, end1, flank, circular) {
  const n = seq.length;
  const start0 = start1 - 1;
  const end0 = end1; // exclusive
  let from = start0 - flank;
  let to = end0 + flank;
  if (!circular) {
    from = clamp(from, 0, n);
    to = clamp(to, 0, n);
    return { seq: seq.slice(from, to), start: from + 1, end: to, wraps: false };
  }
  const len = to - from;
  let out = '';
  for (let i = 0; i < len; i++) out += seq[((from + i) % n + n) % n];
  return { seq: out, start: ((from % n + n) % n) + 1, end: ((to - 1) % n + n) % n + 1, wraps: from < 0 || to > n };
}

function findTRNAs(record) {
  return record.features.filter(f => /^(tRNA|tmRNA)$/i.test(f.type) && f.start && f.end).map((f, idx) => {
    const product = f.qualifiers?.product || f.qualifiers?.gene || f.qualifiers?.note || `tRNA_${idx+1}`;
    return { feature: f, product: Array.isArray(product) ? product.join('; ') : product };
  });
}

function kmerIndex(seq, k) {
  const idx = new Map();
  for (let i = 0; i <= seq.length - k; i++) {
    const kmer = seq.slice(i, i + k);
    if (kmer.includes('N')) continue;
    if (!idx.has(kmer)) idx.set(kmer, []);
    const arr = idx.get(kmer);
    if (arr.length < 250) arr.push(i); // avoid pathological repeats
  }
  return idx;
}

function bestExactCore(querySeq, targetSeq, k, minCore) {
  if (querySeq.length < k || targetSeq.length < k) return null;
  const idx = kmerIndex(querySeq, k);
  let best = null;
  let seedHits = 0;
  for (let j = 0; j <= targetSeq.length - k; j++) {
    const kmer = targetSeq.slice(j, j + k);
    const qs = idx.get(kmer);
    if (!qs) continue;
    for (const i of qs) {
      seedHits++;
      let l = 0;
      while (i - l - 1 >= 0 && j - l - 1 >= 0 && querySeq[i-l-1] === targetSeq[j-l-1]) l++;
      let r = k;
      while (i + r < querySeq.length && j + r < targetSeq.length && querySeq[i+r] === targetSeq[j+r]) r++;
      const len = l + r;
      if (!best || len > best.len) {
        best = { len, qStart: i - l, qEnd: i + r, tStart: j - l, tEnd: j + r, seq: querySeq.slice(i - l, i + r), seedHits };
      }
    }
  }
  if (!best || best.len < minCore) return null;
  best.seedHits = seedHits;
  return best;
}

function regionCandidatesByKmers(querySeq, genomeSeq, k, maxCandidates) {
  const qIndex = kmerIndex(querySeq, k);
  const bins = new Map();
  const binSize = 1000;
  for (let j = 0; j <= genomeSeq.length - k; j++) {
    const kmer = genomeSeq.slice(j, j + k);
    if (!qIndex.has(kmer)) continue;
    const bin = Math.floor(j / binSize);
    bins.set(bin, (bins.get(bin) || 0) + 1);
  }
  return [...bins.entries()]
    .sort((a,b) => b[1] - a[1])
    .slice(0, maxCandidates)
    .map(([bin, count]) => ({ start0: Math.max(0, bin*binSize - 1000), end0: Math.min(genomeSeq.length, (bin+1)*binSize + 1000), count }));
}

function scoreHit(core, opts) {
  const coreScore = Math.min(70, core.len * 2.2);
  const seedScore = Math.min(18, Math.log2(core.seedHits + 1) * 4);
  const trnaScore = opts.searchType === 'tRNA' ? 10 : 0;
  const repeatPenalty = core.seedHits > 300 ? -10 : 0;
  return Math.max(0, Math.min(100, Math.round(coreScore + seedScore + trnaScore + repeatPenalty)));
}

function coordFromWindow(window, offset0, len, genomeLength) {
  const s0 = (window.absoluteStart0 + offset0) % genomeLength;
  const e0 = (s0 + len - 1) % genomeLength;
  return { start: s0 + 1, end: e0 + 1, wraps: s0 + len > genomeLength };
}

function makeCandidateTarget(record, start1, end1, flank, circular, label, feature=null) {
  const w = extractWindow(record.sequence, start1, end1, flank, circular);
  const start0 = w.start - 1;
  return { ...w, absoluteStart0: start0, label, feature };
}

async function readFiles(fileList) {
  const parsed = [];
  for (const file of fileList) {
    const text = await file.text();
    const isGB = /^LOCUS\s/m.test(text) || /\nFEATURES\s+Location\/Qualifiers/.test(text);
    const recs = isGB ? parseGenBank(text, file.name) : parseFasta(text, file.name);
    for (const r of recs) parsed.push({ ...r, fileName: file.name });
  }
  return parsed;
}

function analyzePair(phage, bacterium, settings) {
  const outputs = [];
  const integrases = findIntegrases(phage);
  if (!integrases.length) {
    log(`No integrase-like CDS found in phage ${phage.name}; skipping this phage.`);
    return outputs;
  }
  const trnas = findTRNAs(bacterium);
  if (!trnas.length && bacterium.format === 'GenBank') log(`No tRNA features found in ${bacterium.name}; global search only.`);
  if (bacterium.format === 'FASTA') log(`${bacterium.name} is FASTA; tRNA-first search unavailable.`);

  for (const intHit of integrases.slice(0, 3)) {
    const intF = intHit.feature;
    const pWin = extractWindow(phage.sequence, intF.start, intF.end, settings.intFlank, settings.circular);
    const pSeqPlus = pWin.seq;
    const pSeqMinus = rc(pWin.seq);
    const intName = intHit.text || `${intF.start}..${intF.end}`;

    for (const trna of trnas) {
      const tF = trna.feature;
      const target = makeCandidateTarget(bacterium, tF.start, tF.end, settings.trnaFlank, settings.circular, `near ${trna.product}`, tF);
      evaluateTarget({ outputs, phage, bacterium, intF, intName, pWin, pSeqPlus, pSeqMinus, target, searchType: 'tRNA', product: trna.product, settings });
    }

    if (settings.globalSearch) {
      const globalRegions = regionCandidatesByKmers(pSeqPlus, bacterium.sequence, settings.k, settings.maxGlobal);
      for (const region of globalRegions) {
        const target = { seq: bacterium.sequence.slice(region.start0, region.end0), absoluteStart0: region.start0, start: region.start0 + 1, end: region.end0, wraps: false, label: `global bin support ${region.count}`, feature: null };
        evaluateTarget({ outputs, phage, bacterium, intF, intName, pWin, pSeqPlus, pSeqMinus, target, searchType: 'global', product: '', settings });
      }
    }
  }
  return outputs;
}

function evaluateTarget(ctx) {
  const { outputs, phage, bacterium, intF, intName, pWin, pSeqPlus, pSeqMinus, target, searchType, product, settings } = ctx;
  const tests = [
    { strand: '+', qSeq: pSeqPlus, targetSeq: target.seq },
    { strand: '-', qSeq: pSeqPlus, targetSeq: rc(target.seq) },
    { strand: '+ vs phageRC', qSeq: pSeqMinus, targetSeq: target.seq }
  ];
  let best = null;
  for (const t of tests) {
    const core = bestExactCore(t.qSeq, t.targetSeq, settings.k, settings.minCore);
    if (core && (!best || core.len > best.core.len)) best = { core, ...t };
  }
  if (!best) return;

  const score = scoreHit(best.core, { searchType });
  const hostStart0 = best.strand === '-' ? (target.seq.length - best.core.tEnd) : best.core.tStart;
  const hostCoord = coordFromWindow(target, hostStart0, best.core.len, bacterium.sequence.length);
  const phageCoord = coordFromWindow({ absoluteStart0: pWin.start - 1 }, best.core.qStart, best.core.len, phage.sequence.length);
  outputs.push({
    score,
    phage: phage.name,
    phageFile: phage.fileName,
    host: bacterium.name,
    hostFile: bacterium.fileName,
    integrase: intName,
    integraseCoords: `${intF.start}..${intF.end} (${intF.strand})`,
    searchType,
    hostLocus: target.label,
    strand: best.strand,
    product,
    coreLength: best.core.len,
    coreSequence: best.core.seq,
    seedHits: best.core.seedHits,
    phageCore: `${phageCoord.start}..${phageCoord.end}${phageCoord.wraps ? ' (wraps)' : ''}`,
    hostCore: `${hostCoord.start}..${hostCoord.end}${hostCoord.wraps ? ' (wraps)' : ''}`,
    phageWindow: `${pWin.start}..${pWin.end}${pWin.wraps ? ' (wraps)' : ''}`,
    hostWindow: `${target.start}..${target.end}${target.wraps ? ' (wraps)' : ''}`
  });
}

function dedupeResults(results) {
  const seen = new Map();
  for (const r of results) {
    const key = [r.phage, r.host, r.coreSequence, r.hostCore, r.searchType].join('|');
    if (!seen.has(key) || r.score > seen.get(key).score) seen.set(key, r);
  }
  return [...seen.values()].sort((a,b) => b.score - a.score || b.coreLength - a.coreLength);
}

function renderFiles() {
  const el = $('fileList');
  el.innerHTML = `<strong>Phage files:</strong> ${state.phageFiles.map(f => f.name).join(', ') || 'none'}<br><strong>Bacterial files:</strong> ${state.bactFiles.map(f => f.name).join(', ') || 'none'}`;
}

function renderResults() {
  const tbody = $('resultsTable').querySelector('tbody');
  tbody.innerHTML = '';
  const max = Number($('maxPredictions').value || 150);
  state.results.slice(0, max).forEach((r, i) => {
    const tr = document.createElement('tr');
    tr.innerHTML = `<td>${i+1}</td><td><strong>${r.score}</strong></td><td>${escapeHtml(r.phage)}</td><td>${escapeHtml(r.integrase)}<br><small>${r.integraseCoords}</small></td><td>${escapeHtml(r.host)}</td><td>${r.searchType}</td><td>${escapeHtml(r.hostLocus)}</td><td>${r.strand}</td><td>${escapeHtml(r.product || '')}</td><td>${r.coreLength}</td><td><code>${r.coreSequence}</code><br><small>${r.seedHits} seed hits</small></td><td>phage core ${r.phageCore}<br>host core ${r.hostCore}<br><small>phage window ${r.phageWindow}; host window ${r.hostWindow}</small></td>`;
    tbody.appendChild(tr);
  });
  $('summary').textContent = state.results.length ? `${state.results.length} candidate(s) found. Showing top ${Math.min(max, state.results.length)}.` : 'No candidates passed the current thresholds.';
  const enabled = state.results.length > 0;
  $('downloadTsv').disabled = !enabled; $('downloadJson').disabled = !enabled; $('downloadFasta').disabled = !enabled;
}

function escapeHtml(s) {
  return String(s ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
}

function download(filename, text, type='text/plain') {
  const blob = new Blob([text], { type });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = filename; document.body.appendChild(a); a.click(); a.remove();
  URL.revokeObjectURL(url);
}

function toTSV(rows) {
  const cols = ['score','phage','phageFile','host','hostFile','integrase','integraseCoords','searchType','hostLocus','strand','product','coreLength','coreSequence','seedHits','phageCore','hostCore','phageWindow','hostWindow'];
  const esc = v => String(v ?? '').replace(/\t/g, ' ').replace(/\r?\n/g, ' ');
  return [cols.join('\t'), ...rows.map(r => cols.map(c => esc(r[c])).join('\t'))].join('\n');
}

function toFasta(rows) {
  return rows.map((r, i) => `>candidate_${i+1}|score=${r.score}|${r.phage}|${r.host}|hostCore=${r.hostCore}|phageCore=${r.phageCore}\n${r.coreSequence}`).join('\n');
}

async function runAnalysis() {
  if (!state.phageFiles.length || !state.bactFiles.length) { resetLog('Please add at least one phage file and one bacterial genome file.'); return; }
  resetLog('Parsing files…');
  state.phages = await readFiles(state.phageFiles);
  state.bacteria = await readFiles(state.bactFiles);
  log(`Parsed ${state.phages.length} phage record(s) and ${state.bacteria.length} bacterial record(s).`);

  const settings = {
    intFlank: Number($('intFlank').value || 750),
    trnaFlank: Number($('trnaFlank').value || 750),
    k: Number($('kmerSize').value || 12),
    minCore: Number($('minCore').value || 18),
    maxGlobal: Number($('maxGlobal').value || 25),
    circular: $('assumeCircular').checked,
    globalSearch: $('globalSearch').checked
  };
  const all = [];
  for (const p of state.phages) {
    log(`Phage ${p.name}: ${fmt(p.sequence.length)} bp, ${p.features.length} features.`);
    for (const b of state.bacteria) {
      log(`Comparing ${p.name} to ${b.name} (${fmt(b.sequence.length)} bp)…`);
      all.push(...analyzePair(p, b, settings));
    }
  }
  state.results = dedupeResults(all);
  log(`Finished. ${state.results.length} deduplicated candidate predictions found.`);
  renderResults();
}

function installDrop(zoneId, inputId, targetKey) {
  const zone = $(zoneId); const input = $(inputId);
  input.addEventListener('change', () => { state[targetKey] = [...input.files]; renderFiles(); });
  ['dragenter','dragover'].forEach(evt => zone.addEventListener(evt, e => { e.preventDefault(); zone.classList.add('dragover'); }));
  ['dragleave','drop'].forEach(evt => zone.addEventListener(evt, e => { e.preventDefault(); zone.classList.remove('dragover'); }));
  zone.addEventListener('drop', e => { state[targetKey] = [...e.dataTransfer.files]; renderFiles(); });
}

function makeFile(name, text) { return new File([text], name, { type: 'text/plain' }); }
function loadDemo() {
  const phage = `LOCUS       DemoPhage        5000 bp    DNA     linear   PHG 01-JAN-2026\nFEATURES             Location/Qualifiers\n     CDS             1950..2400\n                     /gene="int"\n                     /product="tyrosine integrase"\nORIGIN\n        1 ${'A'.repeat(1880)}GATCCGTTACGATCGACTGATGACCTGAACTGACCGGTA${'C'.repeat(350)}${'G'.repeat(2727)}\n//`;
  const host = `LOCUS       DemoHost        8000 bp    DNA     circular BCT 01-JAN-2026\nFEATURES             Location/Qualifiers\n     tRNA            3900..3975\n                     /product="tRNA-Lys"\nORIGIN\n        1 ${'T'.repeat(3850)}GATCCGTTACGATCGACTGATGACCTGAACTGACCGGTA${'A'.repeat(4108)}\n//`;
  state.phageFiles = [makeFile('demo_phage.gbk', phage)];
  state.bactFiles = [makeFile('demo_host.gbk', host)];
  renderFiles(); resetLog('Tiny demo loaded. Click Run attP/attB prediction.');
}

installDrop('phageDrop', 'phageFiles', 'phageFiles');
installDrop('bactDrop', 'bactFiles', 'bactFiles');
$('runBtn').addEventListener('click', runAnalysis);
$('demoBtn').addEventListener('click', loadDemo);
$('clearBtn').addEventListener('click', () => { state.phageFiles=[]; state.bactFiles=[]; state.results=[]; renderFiles(); renderResults(); resetLog(); });
$('downloadTsv').addEventListener('click', () => download('attp_attb_predictions.tsv', toTSV(state.results), 'text/tab-separated-values'));
$('downloadJson').addEventListener('click', () => download('attp_attb_predictions.json', JSON.stringify(state.results, null, 2), 'application/json'));
$('downloadFasta').addEventListener('click', () => download('attp_attb_candidate_cores.fasta', toFasta(state.results), 'text/plain'));
renderFiles();
