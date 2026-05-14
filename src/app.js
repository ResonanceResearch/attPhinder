/* attP / attB Finder UI
   Heavy parsing and genome comparison are run in a Web Worker so large bacterial
   GenBank files do not freeze the browser tab.
*/
const state = {
  phageFiles: [],
  bactFiles: [],
  results: [],
  currentPhageIndex: 0,
  worker: null
};

const $ = (id) => document.getElementById(id);
const logBox = $('log');

function log(msg) {
  const stamp = new Date().toLocaleTimeString();
  logBox.textContent += `\n[${stamp}] ${msg}`;
  logBox.scrollTop = logBox.scrollHeight;
}
function resetLog(msg = 'Waiting for files…') { logBox.textContent = msg; }
function fmtBytes(n) {
  if (!Number.isFinite(n)) return '';
  if (n < 1024) return `${n} B`;
  if (n < 1024*1024) return `${(n/1024).toFixed(1)} KB`;
  return `${(n/1024/1024).toFixed(1)} MB`;
}
function safeName(name) { return name.replace(/[^a-z0-9._-]+/gi, '_'); }
function compactValue(v) { return String(v ?? '').replace(/[|\t\r\n]+/g, ' ').trim(); }
function fastaValue(v) { return compactValue(v).replace(/\s+/g, '_') || 'NA'; }
function candidateId(r, i) {
  const phage = safeName(r.phage || r.phageFile || 'phage');
  const host = safeName(r.host || r.hostFile || 'host');
  const coord = safeName(`${r.hostCore || 'host'}_${r.phageCore || 'phage'}`);
  return `candidate_${i + 1}_${phage}_${host}_${coord}`;
}
function featureSummary(r) {
  const prox = r.featureProximity || 'unknown';
  const label = r.nearestFeatureLabel || r.product || '';
  const type = r.nearestFeatureType || '';
  const dist = r.distanceToFeatureBp === null || r.distanceToFeatureBp === undefined ? 'NA' : `${r.distanceToFeatureBp} bp`;
  if (!type && !label) return prox;
  return `${prox}; ${type || 'feature'} ${label ? '(' + label + ')' : ''}; distance ${dist}`;
}


function renderFiles() {
  const el = $('fileList');
  const phage = state.phageFiles.map(f => `${f.name} (${fmtBytes(f.size)})`).join(', ') || 'none';
  const bact = state.bactFiles.map(f => `${f.name} (${fmtBytes(f.size)})`).join(', ') || 'none';
  el.innerHTML = `<strong>Phage files:</strong> ${escapeHtml(phage)}<br><strong>Bacterial files:</strong> ${escapeHtml(bact)}`;

  const biggest = Math.max(0, ...state.bactFiles.map(f => f.size || 0));
  const warn = $('largeFileWarning');
  if (biggest > 15 * 1024 * 1024) {
    warn.hidden = false;
    warn.textContent = `Large bacterial file detected (${fmtBytes(biggest)}). For large host genomes, use a prebuilt host_index.json from tools/build_host_index_interactive.py. Direct GenBank parsing remains available for smaller files but can still be fragile in browsers.`;
  } else {
    warn.hidden = true;
  }
}

function getPhageGroups() {
  const map = new Map();
  for (const r of state.results) {
    const key = `${r.phageFile || ''}||${r.phage || ''}`;
    if (!map.has(key)) map.set(key, { key, label: r.phage || r.phageFile || 'unknown phage', file: r.phageFile || '', rows: [] });
    map.get(key).rows.push(r);
  }
  return [...map.values()].sort((a, b) => a.label.localeCompare(b.label));
}

function getCurrentPhageGroup() {
  const groups = getPhageGroups();
  if (!groups.length) return null;
  state.currentPhageIndex = Math.max(0, Math.min(state.currentPhageIndex, groups.length - 1));
  return groups[state.currentPhageIndex];
}

function ensurePhageNav() {
  let nav = $('phageNav');
  if (nav) return nav;
  const summary = $('summary');
  nav = document.createElement('div');
  nav.id = 'phageNav';
  nav.className = 'phage-nav';
  nav.innerHTML = `
    <button id="prevPhage" class="secondary" disabled>Previous phage</button>
    <span id="phageNavLabel">No phage selected</span>
    <button id="nextPhage" class="secondary" disabled>Next phage</button>
    <button id="downloadCurrentTsv" class="secondary" disabled>Download current phage TSV</button>
    <button id="downloadCurrentJson" class="secondary" disabled>Download current phage JSON</button>
    <button id="downloadCurrentFasta" class="secondary" disabled>Download current phage FASTA</button>
  `;
  summary.insertAdjacentElement('afterend', nav);
  $('prevPhage').addEventListener('click', () => { state.currentPhageIndex--; renderResults(); });
  $('nextPhage').addEventListener('click', () => { state.currentPhageIndex++; renderResults(); });
  $('downloadCurrentTsv').addEventListener('click', () => {
    const group = getCurrentPhageGroup();
    if (group) download(`attp_attb_${safeName(group.label)}.tsv`, toTSV(group.rows), 'text/tab-separated-values');
  });
  $('downloadCurrentJson').addEventListener('click', () => {
    const group = getCurrentPhageGroup();
    if (group) download(`attp_attb_${safeName(group.label)}.json`, JSON.stringify(group.rows, null, 2), 'application/json');
  });
  $('downloadCurrentFasta').addEventListener('click', () => {
    const group = getCurrentPhageGroup();
    if (group) download(`attp_attb_${safeName(group.label)}_cores.fasta`, toFasta(group.rows), 'text/plain');
  });
  return nav;
}

function renderResults() {
  const tbody = $('resultsTable').querySelector('tbody');
  tbody.innerHTML = '';
  const max = Number($('maxPredictions').value || 150);
  const nav = ensurePhageNav();
  const groups = getPhageGroups();
  const group = getCurrentPhageGroup();
  const visibleRows = group ? group.rows : [];

  visibleRows.slice(0, max).forEach((r, i) => {
    const tr = document.createElement('tr');
    const queryLabel = r.queryMode === 'whole_fasta_phage' ? 'Whole FASTA query' : 'Integrase window';
    tr.innerHTML = `<td>${i+1}</td><td><strong>${r.score}</strong></td><td>${escapeHtml(r.phage)}<br><small>${escapeHtml(queryLabel)}</small></td><td>${escapeHtml(r.integrase)}<br><small>${escapeHtml(r.integraseCoords)}</small></td><td>${escapeHtml(r.host)}${r.hostContig ? `<br><small>contig: ${escapeHtml(r.hostContig)}</small>` : ''}</td><td>${escapeHtml(r.searchType)}${r.duplicateCount && r.duplicateCount > 1 ? `<br><small>${r.duplicateCount} duplicate paths collapsed</small>` : ''}</td><td>${escapeHtml(r.hostLocus)}</td><td>${escapeHtml(r.strand)}</td><td>${escapeHtml(featureSummary(r))}<br><small>${escapeHtml(r.nearestFeatureCoords || '')}</small></td><td>${r.coreLength}</td><td><code>${escapeHtml(r.coreSequence)}</code><br><small>${r.seedHits} seed hits</small></td><td>attP/phage core ${escapeHtml(r.phageCore)}<br>attB/host core ${escapeHtml(r.hostCore)}<br><small>phage window ${escapeHtml(r.phageWindow)}; host window ${escapeHtml(r.hostWindow)}</small></td>`;
    tbody.appendChild(tr);
  });

  if (state.results.length && group) {
    $('summary').textContent = `${state.results.length} total deduplicated candidate(s) across ${groups.length} phage record(s). Showing ${Math.min(max, visibleRows.length)} of ${visibleRows.length} for current phage.`;
    $('phageNavLabel').textContent = `Phage ${state.currentPhageIndex + 1} of ${groups.length}: ${group.label}`;
  } else {
    $('summary').textContent = 'No candidates passed the current thresholds.';
    $('phageNavLabel').textContent = 'No phage selected';
  }

  const enabled = state.results.length > 0;
  $('downloadTsv').disabled = !enabled; $('downloadJson').disabled = !enabled; $('downloadFasta').disabled = !enabled;
  $('downloadCurrentTsv').disabled = !group;
  $('downloadCurrentJson').disabled = !group;
  $('downloadCurrentFasta').disabled = !group;
  $('prevPhage').disabled = !group || state.currentPhageIndex <= 0;
  $('nextPhage').disabled = !group || state.currentPhageIndex >= groups.length - 1;
  nav.hidden = !enabled;
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
  const cols = [
    'candidateId','score','phage','phageFile','phageLength','host','hostFile','hostContig','hostLength',
    'integrase','integraseCoords','queryMode','searchType','searchSources','duplicateCount',
    'hostLocus','strand','product','nearestFeatureType','nearestFeatureLabel','nearestFeatureCoords',
    'distanceToFeatureBp','featureProximity','coreLength','coreSequence','seedHits',
    'phageCore','hostCore','phageWindow','hostWindow'
  ];
  const esc = v => String(v ?? '').replace(/\t/g, ' ').replace(/\r?\n/g, ' ');
  const value = (r, c, i) => c === 'candidateId' ? candidateId(r, i) : (Array.isArray(r[c]) ? r[c].join(';') : r[c]);
  return [cols.join('\t'), ...rows.map((r, i) => cols.map(c => esc(value(r, c, i))).join('\t'))].join('\n');
}

function toFasta(rows) {
  return rows.map((r, i) => {
    const id = candidateId(r, i);
    const fields = [
      `score=${r.score}`,
      `phage=${fastaValue(r.phage)}`,
      `phageFile=${fastaValue(r.phageFile)}`,
      `queryMode=${fastaValue(r.queryMode)}`,
      `integraseCoords=${fastaValue(r.integraseCoords)}`,
      `host=${fastaValue(r.host)}`,
      `hostFile=${fastaValue(r.hostFile)}`,
      `hostContig=${fastaValue(r.hostContig)}`,
      `attP_phageCore=${fastaValue(r.phageCore)}`,
      `attB_hostCore=${fastaValue(r.hostCore)}`,
      `phageWindow=${fastaValue(r.phageWindow)}`,
      `hostWindow=${fastaValue(r.hostWindow)}`,
      `coreLength=${r.coreLength}`,
      `strand=${fastaValue(r.strand)}`,
      `searchType=${fastaValue(r.searchType)}`,
      `featureProximity=${fastaValue(r.featureProximity)}`,
      `nearestFeatureType=${fastaValue(r.nearestFeatureType)}`,
      `nearestFeatureLabel=${fastaValue(r.nearestFeatureLabel || r.product)}`,
      `nearestFeatureCoords=${fastaValue(r.nearestFeatureCoords)}`,
      `distanceToFeatureBp=${r.distanceToFeatureBp ?? 'NA'}`,
      `seedHits=${r.seedHits}`,
      `duplicatePaths=${r.duplicateCount || 1}`
    ];
    return `>${id}|${fields.join('|')}\n${r.coreSequence}`;
  }).join('\n');
}

function collectSettings() {
  return {
    intFlank: Number($('intFlank').value || 750),
    trnaFlank: Number($('trnaFlank').value || 750),
    k: Number($('kmerSize').value || 12),
    minCore: Number($('minCore').value || 18),
    maxGlobal: Number($('maxGlobal').value || 25),
    circular: $('assumeCircular').checked,
    globalSearch: $('globalSearch').checked,
    maxPredictions: Number($('maxPredictions').value || 150)
  };
}

function terminateWorker() {
  if (state.worker) {
    state.worker.terminate();
    state.worker = null;
  }
}

async function runAnalysis() {
  if (!state.phageFiles.length || !state.bactFiles.length) {
    resetLog('Please add at least one phage file and one bacterial genome file.');
    return;
  }
  terminateWorker();
  state.results = [];
  state.currentPhageIndex = 0;
  renderResults();
  resetLog('Starting worker…');
  $('runBtn').disabled = true;
  $('stopBtn').disabled = false;

  const worker = new Worker('src/worker.js');
  state.worker = worker;
  worker.onmessage = (event) => {
    const msg = event.data || {};
    if (msg.type === 'log') log(msg.message);
    if (msg.type === 'results') {
      state.results = msg.results || [];
      state.currentPhageIndex = 0;
      renderResults();
    }
    if (msg.type === 'done') {
      log(`Finished. ${state.results.length} deduplicated candidate predictions found.`);
      $('runBtn').disabled = false;
      $('stopBtn').disabled = true;
      terminateWorker();
    }
    if (msg.type === 'error') {
      log(`ERROR: ${msg.message}`);
      $('runBtn').disabled = false;
      $('stopBtn').disabled = true;
      terminateWorker();
    }
  };
  worker.onerror = (err) => {
    log(`ERROR: ${err.message || 'Worker failed.'}`);
    $('runBtn').disabled = false;
    $('stopBtn').disabled = true;
    terminateWorker();
  };
  worker.postMessage({ type: 'run', phageFiles: state.phageFiles, bactFiles: state.bactFiles, settings: collectSettings() });
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
  const sharedCore = 'GATCCGTTACGATCGACTGATGACCTGAACTGACCGGTA';
  const phage = `LOCUS       DemoPhage        5000 bp    DNA     linear   PHG 01-JAN-2026\nFEATURES             Location/Qualifiers\n     CDS             1950..2400\n                     /gene="int"\n                     /product="tyrosine integrase"\nORIGIN\n        1 ${'A'.repeat(1880)}${sharedCore}${'C'.repeat(350)}${'G'.repeat(2727)}\n//`;
  const hostIndex = {
    schema: 'att-site-host-index-v1',
    created_utc: new Date().toISOString(),
    metadata: { name: 'DemoHost host_index.json', organism: 'Demo host' },
    parameters: { flank: 500, feature_types: ['tRNA'], global_chunks_included: false },
    summary: { contig_count: 1, total_bp: 8000, feature_neighborhood_count: 1, global_chunk_count: 0 },
    contigs: [{ id: 'DemoHost', description: 'synthetic demo', length: 8000, topology: 'circular' }],
    feature_neighborhoods: [{
      id: 'DemoHost:tRNA:3900-3975:1',
      contig_id: 'DemoHost',
      feature_type: 'tRNA',
      feature_label: 'tRNA-Lys',
      feature_start: 3900,
      feature_end: 3975,
      feature_strand: 1,
      window_start_approx: 3400,
      window_end_approx: 4475,
      window_is_circular_wrapped: false,
      sequence: `${'T'.repeat(450)}${sharedCore}${'A'.repeat(581)}`,
      length: 1076
    }],
    global_chunks: []
  };
  state.phageFiles = [makeFile('demo_phage.gbk', phage)];
  state.bactFiles = [makeFile('demo_host_index.json', JSON.stringify(hostIndex))];
  renderFiles(); resetLog('Tiny host_index.json demo loaded. Click Run attP/attB prediction.');
}


installDrop('phageDrop', 'phageFiles', 'phageFiles');
installDrop('bactDrop', 'bactFiles', 'bactFiles');
$('runBtn').addEventListener('click', runAnalysis);
$('stopBtn').addEventListener('click', () => { terminateWorker(); $('runBtn').disabled = false; $('stopBtn').disabled = true; log('Stopped by user.'); });
$('demoBtn').addEventListener('click', loadDemo);
$('clearBtn').addEventListener('click', () => { terminateWorker(); state.phageFiles=[]; state.bactFiles=[]; state.results=[]; state.currentPhageIndex=0; renderFiles(); renderResults(); resetLog(); $('runBtn').disabled=false; $('stopBtn').disabled=true; });
$('downloadTsv').addEventListener('click', () => download('attp_attb_predictions.tsv', toTSV(state.results), 'text/tab-separated-values'));
$('downloadJson').addEventListener('click', () => download('attp_attb_predictions.json', JSON.stringify(state.results, null, 2), 'application/json'));
$('downloadFasta').addEventListener('click', () => download('attp_attb_candidate_cores.fasta', toFasta(state.results), 'text/plain'));
renderFiles();
