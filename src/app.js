/* attP / attB Finder UI
   Heavy parsing and genome comparison are run in a Web Worker so large bacterial
   GenBank files do not freeze the browser tab.
*/
const state = {
  phageFiles: [],
  bactFiles: [],
  results: [],
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

function renderResults() {
  const tbody = $('resultsTable').querySelector('tbody');
  tbody.innerHTML = '';
  const max = Number($('maxPredictions').value || 150);
  state.results.slice(0, max).forEach((r, i) => {
    const tr = document.createElement('tr');
    tr.innerHTML = `<td>${i+1}</td><td><strong>${r.score}</strong></td><td>${escapeHtml(r.phage)}</td><td>${escapeHtml(r.integrase)}<br><small>${escapeHtml(r.integraseCoords)}</small></td><td>${escapeHtml(r.host)}</td><td>${escapeHtml(r.searchType)}${r.duplicateCount && r.duplicateCount > 1 ? `<br><small>${r.duplicateCount} duplicate paths collapsed</small>` : ''}</td><td>${escapeHtml(r.hostLocus)}</td><td>${escapeHtml(r.strand)}</td><td>${escapeHtml(r.product || '')}</td><td>${r.coreLength}</td><td><code>${escapeHtml(r.coreSequence)}</code><br><small>${r.seedHits} seed hits</small></td><td>phage core ${escapeHtml(r.phageCore)}<br>host core ${escapeHtml(r.hostCore)}<br><small>phage window ${escapeHtml(r.phageWindow)}; host window ${escapeHtml(r.hostWindow)}</small></td>`;
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
  const cols = ['score','phage','phageFile','host','hostFile','integrase','integraseCoords','searchType','searchSources','duplicateCount','hostLocus','strand','product','coreLength','coreSequence','seedHits','phageCore','hostCore','phageWindow','hostWindow'];
  const esc = v => String(v ?? '').replace(/\t/g, ' ').replace(/\r?\n/g, ' ');
  return [cols.join('\t'), ...rows.map(r => cols.map(c => esc(r[c])).join('\t'))].join('\n');
}

function toFasta(rows) {
  return rows.map((r, i) => `>candidate_${i+1}|score=${r.score}|${r.phage}|${r.host}|hostCore=${r.hostCore}|phageCore=${r.phageCore}\n${r.coreSequence}`).join('\n');
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
$('clearBtn').addEventListener('click', () => { terminateWorker(); state.phageFiles=[]; state.bactFiles=[]; state.results=[]; renderFiles(); renderResults(); resetLog(); $('runBtn').disabled=false; $('stopBtn').disabled=true; });
$('downloadTsv').addEventListener('click', () => download('attp_attb_predictions.tsv', toTSV(state.results), 'text/tab-separated-values'));
$('downloadJson').addEventListener('click', () => download('attp_attb_predictions.json', JSON.stringify(state.results, null, 2), 'application/json'));
$('downloadFasta').addEventListener('click', () => download('attp_attb_candidate_cores.fasta', toFasta(state.results), 'text/plain'));
renderFiles();
