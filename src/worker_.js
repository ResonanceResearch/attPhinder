/* attP / attB Finder Worker
   Optimized for large host GenBank files by parsing only host tRNA/tmRNA features.
*/
function postLog(message) { self.postMessage({ type: 'log', message }); }
function rc(seq) {
  const out = new Array(seq.length);
  for (let i = seq.length - 1, j = 0; i >= 0; i--, j++) {
    const b = seq[i];
    out[j] = b === 'A' ? 'T' : b === 'T' ? 'A' : b === 'G' ? 'C' : b === 'C' ? 'G' : 'N';
  }
  return out.join('');
}
function cleanSeq(s) { return (s || '').replace(/[^A-Za-z]/g, '').toUpperCase().replace(/U/g, 'T').replace(/[^ACGTN]/g, 'N'); }
function clamp(n, lo, hi) { return Math.max(lo, Math.min(hi, n)); }
function fmt(n) { return Number.isFinite(n) ? n.toLocaleString() : ''; }
function fmtBytes(n) {
  if (!Number.isFinite(n)) return '';
  if (n < 1024) return `${n} B`;
  if (n < 1024*1024) return `${(n/1024).toFixed(1)} KB`;
  return `${(n/1024/1024).toFixed(1)} MB`;
}

function parseFasta(text, fileName) {
  const records = [];
  let header = fileName, seqParts = [];
  const lines = text.split(/\r?\n/);
  for (const line of lines) {
    if (line.startsWith('>')) {
      if (seqParts.length) records.push({ name: header, fileName, sequence: cleanSeq(seqParts.join('')), features: [], format: 'FASTA' });
      header = line.slice(1).trim() || fileName;
      seqParts = [];
    } else {
      seqParts.push(line.trim());
    }
  }
  if (seqParts.length) records.push({ name: header, fileName, sequence: cleanSeq(seqParts.join('')), features: [], format: 'FASTA' });
  if (!records.length) records.push({ name: fileName, fileName, sequence: cleanSeq(text), features: [], format: 'FASTA' });
  return records;
}

function splitGenBankRecords(text) {
  const parts = text.split(/^\/\/\s*$/m).map(x => x.trim()).filter(Boolean);
  return parts.length ? parts : [text];
}

function getNameFromRecord(rec, fileName) {
  const locusMatch = rec.match(/^LOCUS\s+(\S+)/m);
  if (locusMatch) return locusMatch[1];
  const defMatch = rec.match(/^DEFINITION\s+([\s\S]*?)(?=^ACCESSION\s|^VERSION\s|^KEYWORDS\s|^SOURCE\s)/m);
  return defMatch?.[1]?.replace(/\s+/g, ' ').trim() || fileName;
}

function extractOrigin(rec) {
  const originMatch = rec.match(/^ORIGIN\b/m);
  if (!originMatch) return '';
  const originStart = originMatch.index + originMatch[0].length;
  let originText = rec.slice(originStart);
  const end = originText.search(/^\/\//m);
  if (end >= 0) originText = originText.slice(0, end);
  return cleanSeq(originText);
}

function extractFeaturesBlock(rec) {
  const featuresMatch = rec.match(/^FEATURES\s+Location\/Qualifiers/m);
  if (!featuresMatch) return '';
  const start = featuresMatch.index + featuresMatch[0].length;
  const after = rec.slice(start);
  const endMatch = after.search(/^ORIGIN|^CONTIG|^BASE COUNT|^REFERENCE|^COMMENT|^SOURCE|^\/\//m);
  return endMatch >= 0 ? after.slice(0, endMatch) : after;
}

function parseGenBank(text, fileName, role) {
  const records = [];
  const keepTypes = role === 'host' ? new Set(['tRNA', 'tmRNA']) : new Set(['CDS', 'tRNA', 'tmRNA']);
  for (const rec of splitGenBankRecords(text)) {
    const name = getNameFromRecord(rec, fileName);
    const sequence = extractOrigin(rec);
    if (!sequence) continue;
    const block = extractFeaturesBlock(rec);
    const features = parseFeaturesFiltered(block, keepTypes);
    records.push({ name, fileName, sequence, features, format: 'GenBank' });
  }
  return records;
}

function parseFeaturesFiltered(block, keepTypes) {
  const lines = block.split(/\r?\n/);
  const features = [];
  let current = null;
  let keep = false;
  for (const line of lines) {
    const start = line.match(/^\s{5}(\S+)\s+(.+)/);
    if (start) {
      if (current && keep) features.push(finalizeFeature(current));
      keep = keepTypes.has(start[1]);
      current = keep ? { type: start[1], locationRaw: start[2].trim(), qualifiersRaw: [] } : null;
      continue;
    }
    if (!current || !keep) continue;
    const cont = line.length > 21 ? line.slice(21).trim() : '';
    if (!cont) continue;
    if (cont.startsWith('/')) current.qualifiersRaw.push(cont);
    else if (current.qualifiersRaw.length && !cont.match(/^[a-zA-Z_]+\(/) && !cont.match(/^[<>]?\d/)) current.qualifiersRaw[current.qualifiersRaw.length - 1] += ' ' + cont;
    else current.locationRaw += cont;
  }
  if (current && keep) features.push(finalizeFeature(current));
  return features;
}

function finalizeFeature(f) {
  const qualifiers = {};
  for (const q of f.qualifiersRaw) {
    const m = q.match(/^\/(\S+?)(?:=(.*))?$/);
    if (!m) continue;
    let val = m[2] || true;
    if (typeof val === 'string') val = val.replace(/^"|"$/g, '').replace(/"\s+"/g, ' ');
    if (qualifiers[m[1]]) {
      if (!Array.isArray(qualifiers[m[1]])) qualifiers[m[1]] = [qualifiers[m[1]]];
      qualifiers[m[1]].push(val);
    } else qualifiers[m[1]] = val;
  }
  const loc = parseLocation(f.locationRaw);
  return { type: f.type, locationRaw: f.locationRaw, qualifiers, ...loc };
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
  return record.features.filter(f => f.type === 'CDS' && f.start && f.end).map(f => {
    const t = qualText(f).toLowerCase();
    let score = 0;
    if (/\bintegrase\b|\bint\b/.test(t)) score += 60;
    if (/recombinase|tyrosine|serine|site-specific/.test(t)) score += 25;
    if (/transposase|terminase|portal|capsid|tail|hypothetical/.test(t)) score -= /hypothetical/.test(t) ? 5 : 25;
    return { feature: f, text: qualText(f), score };
  }).filter(c => c.score > 0).sort((a,b) => b.score - a.score);
}

function extractWindow(seq, start1, end1, flank, circular) {
  const n = seq.length;
  const start0 = start1 - 1;
  const end0 = end1;
  let from = start0 - flank;
  let to = end0 + flank;
  if (!circular) {
    from = clamp(from, 0, n);
    to = clamp(to, 0, n);
    return { seq: seq.slice(from, to), start: from + 1, end: to, wraps: false };
  }
  const len = to - from;
  let out = '';
  const chunks = [];
  for (let i = 0; i < len; i += 10000) {
    const partLen = Math.min(10000, len - i);
    let part = '';
    for (let x = 0; x < partLen; x++) part += seq[((from + i + x) % n + n) % n];
    chunks.push(part);
  }
  out = chunks.join('');
  return { seq: out, start: ((from % n + n) % n) + 1, end: ((to - 1) % n + n) % n + 1, wraps: from < 0 || to > n };
}

function findTRNAs(record) {
  return record.features.filter(f => /^(tRNA|tmRNA)$/i.test(f.type) && f.start && f.end).map((f, idx) => {
    const product = f.qualifiers?.product || f.qualifiers?.gene || f.qualifiers?.note || `${f.type}_${idx+1}`;
    return { feature: f, product: Array.isArray(product) ? product.join('; ') : product };
  });
}

function kmerIndex(seq, k) {
  const idx = new Map();
  for (let i = 0; i <= seq.length - k; i++) {
    const kmer = seq.slice(i, i + k);
    if (kmer.includes('N')) continue;
    let arr = idx.get(kmer);
    if (!arr) { arr = []; idx.set(kmer, arr); }
    if (arr.length < 250) arr.push(i);
  }
  return idx;
}

function bestExactCoreFromIndex(querySeq, qIndex, targetSeq, k, minCore) {
  if (querySeq.length < k || targetSeq.length < k) return null;
  let best = null;
  let seedHits = 0;
  for (let j = 0; j <= targetSeq.length - k; j++) {
    const kmer = targetSeq.slice(j, j + k);
    const qs = qIndex.get(kmer);
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

function regionCandidatesByKmers(qIndex, genomeSeq, k, maxCandidates) {
  if (maxCandidates <= 0) return [];
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

async function readFiles(fileList, role) {
  const parsed = [];
  for (const file of fileList) {
    postLog(`Reading ${role} file ${file.name} (${fmtBytes(file.size)})…`);
    const text = await file.text();
    const looksLikeHostIndex = role === 'host' && (file.name.toLowerCase().endsWith('.json') || text.trim().startsWith('{'));
    let recs;
    if (looksLikeHostIndex) {
      recs = [parseHostIndex(text, file.name)];
    } else {
      const isGB = /^LOCUS\s/m.test(text) || /\nFEATURES\s+Location\/Qualifiers/.test(text);
      recs = isGB ? parseGenBank(text, file.name, role) : parseFasta(text, file.name);
    }
    for (const r of recs) parsed.push(r);
    postLog(`Parsed ${recs.length} record(s) from ${file.name}.`);
  }
  return parsed;
}

function parseHostIndex(text, fileName) {
  let idx;
  try { idx = JSON.parse(text); }
  catch (e) { throw new Error(`${fileName} is not valid JSON: ${e.message}`); }
  if (!idx || idx.schema !== 'att-site-host-index-v1') {
    throw new Error(`${fileName} does not look like an att-site-host-index-v1 file. Rebuild it with the included tools/build_host_index_interactive.py script.`);
  }
  const metadata = idx.metadata || {};
  const summary = idx.summary || {};
  const contigs = new Map((idx.contigs || []).map(c => [c.id, c]));
  const featureNeighborhoods = (idx.feature_neighborhoods || []).filter(n => n.sequence).map((n) => {
    const contig = contigs.get(n.contig_id) || {};
    const start = Number(n.window_start_approx || 1);
    const end = Number(n.window_end_approx || start + (n.sequence || '').length - 1);
    return {
      seq: cleanSeq(n.sequence),
      absoluteStart0: Math.max(0, start - 1),
      start,
      end,
      wraps: !!n.window_is_circular_wrapped,
      label: `${n.feature_type || 'feature'} ${n.feature_start || ''}..${n.feature_end || ''} on ${n.contig_id || 'contig'} (${n.feature_label || 'unlabeled'})`,
      feature: null,
      product: n.feature_label || '',
      contigId: n.contig_id || '',
      contigLength: Number(contig.length || summary.total_bp || (n.sequence || '').length)
    };
  });
  const globalChunks = (idx.global_chunks || []).filter(c => c.sequence).map((c) => {
    const contig = contigs.get(c.contig_id) || {};
    return {
      seq: cleanSeq(c.sequence),
      absoluteStart0: Math.max(0, Number(c.start || 1) - 1),
      start: Number(c.start || 1),
      end: Number(c.end || (c.sequence || '').length),
      wraps: false,
      label: `indexed global chunk ${c.start || '?'}..${c.end || '?'} on ${c.contig_id || 'contig'}`,
      feature: null,
      product: '',
      contigId: c.contig_id || '',
      contigLength: Number(contig.length || summary.total_bp || (c.sequence || '').length)
    };
  });
  if (!featureNeighborhoods.length) postLog(`${fileName}: host index contains no sequence-bearing feature neighborhoods.`);
  if ((idx.global_chunks || []).length && !globalChunks.length) postLog(`${fileName}: global chunks are present but lack sequences, so broad search cannot use this index. Rebuild without --omit-global-sequences if broad browser search is desired.`);
  return {
    name: metadata.name || metadata.organism || fileName,
    fileName,
    sequence: '',
    sequenceLength: Number(summary.total_bp || 0),
    features: [],
    format: 'HostIndex',
    hostIndex: { featureNeighborhoods, globalChunks, parameters: idx.parameters || {}, contigs: idx.contigs || [], summary }
  };
}

function countKmerSupport(qIndex, targetSeq, k) {
  let count = 0;
  for (let j = 0; j <= targetSeq.length - k; j++) {
    const kmer = targetSeq.slice(j, j + k);
    if (qIndex.has(kmer)) count++;
  }
  return count;
}


function analyzePair(phage, bacterium, settings) {
  const outputs = [];
  const integrases = findIntegrases(phage);
  if (!integrases.length) {
    postLog(`No integrase-like CDS found in phage ${phage.name}; skipping this phage.`);
    return outputs;
  }

  const usingHostIndex = bacterium.format === 'HostIndex';
  const trnas = usingHostIndex ? [] : findTRNAs(bacterium);
  if (usingHostIndex) {
    const nFeat = bacterium.hostIndex.featureNeighborhoods.length;
    const nGlob = bacterium.hostIndex.globalChunks.length;
    postLog(`${bacterium.name} is a host_index.json: ${nFeat} indexed tRNA/tmRNA/rRNA neighborhood(s), ${nGlob} sequence-bearing global chunk(s).`);
  } else {
    if (!trnas.length && bacterium.format === 'GenBank') postLog(`No tRNA/tmRNA features found in ${bacterium.name}; global search only.`);
    if (bacterium.format === 'FASTA') postLog(`${bacterium.name} is FASTA; tRNA-first search unavailable.`);
  }

  for (const intHit of integrases.slice(0, 3)) {
    const intF = intHit.feature;
    const pWin = extractWindow(phage.sequence, intF.start, intF.end, settings.intFlank, settings.circular);
    const pSeqPlus = pWin.seq;
    const pSeqMinus = rc(pWin.seq);
    const pPlusIndex = kmerIndex(pSeqPlus, settings.k);
    const pMinusIndex = kmerIndex(pSeqMinus, settings.k);
    const intName = intHit.text || `${intF.start}..${intF.end}`;

    if (usingHostIndex) {
      const targets = bacterium.hostIndex.featureNeighborhoods;
      postLog(`Testing ${targets.length} indexed host feature neighborhood(s) for ${bacterium.name} against ${phage.name}.`);
      for (const target of targets) {
        evaluateTarget({ outputs, phage, bacterium, intF, intName, pWin, pSeqPlus, pSeqMinus, pPlusIndex, pMinusIndex, target, searchType: 'indexed_feature', product: target.product || '', settings });
      }
      if (settings.globalSearch) {
        const chunks = bacterium.hostIndex.globalChunks;
        if (!chunks.length) {
          postLog(`No sequence-bearing global chunks available in ${bacterium.fileName}; skipping broad search. Rebuild host_index.json with broader chunks enabled if needed.`);
        } else {
          postLog(`Ranking ${chunks.length} indexed global chunk(s); testing top ${settings.maxGlobal}.`);
          const ranked = chunks.map(ch => ({ target: ch, support: countKmerSupport(pPlusIndex, ch.seq, settings.k) }))
            .filter(x => x.support > 0)
            .sort((a,b) => b.support - a.support)
            .slice(0, settings.maxGlobal);
          for (const { target, support } of ranked) {
            const t = { ...target, label: `${target.label}; k-mer support ${support}` };
            evaluateTarget({ outputs, phage, bacterium, intF, intName, pWin, pSeqPlus, pSeqMinus, pPlusIndex, pMinusIndex, target: t, searchType: 'indexed_global', product: '', settings });
          }
        }
      }
      continue;
    }

    postLog(`Testing ${trnas.length} tRNA/tmRNA region(s) for ${bacterium.name} against ${phage.name}.`);
    for (const trna of trnas) {
      const tF = trna.feature;
      const target = makeCandidateTarget(bacterium, tF.start, tF.end, settings.trnaFlank, settings.circular, `near ${trna.product}`, tF);
      target.product = trna.product;
      target.contigLength = bacterium.sequence.length;
      evaluateTarget({ outputs, phage, bacterium, intF, intName, pWin, pSeqPlus, pSeqMinus, pPlusIndex, pMinusIndex, target, searchType: 'tRNA', product: trna.product, settings });
    }

    if (settings.globalSearch) {
      postLog(`Running broader k-mer scan for ${bacterium.name} (${fmt(bacterium.sequence.length)} bp); max ${settings.maxGlobal} candidate bins.`);
      const globalRegions = regionCandidatesByKmers(pPlusIndex, bacterium.sequence, settings.k, settings.maxGlobal);
      for (const region of globalRegions) {
        const target = { seq: bacterium.sequence.slice(region.start0, region.end0), absoluteStart0: region.start0, start: region.start0 + 1, end: region.end0, wraps: false, label: `global bin support ${region.count}`, feature: null, contigLength: bacterium.sequence.length };
        evaluateTarget({ outputs, phage, bacterium, intF, intName, pWin, pSeqPlus, pSeqMinus, pPlusIndex, pMinusIndex, target, searchType: 'global', product: '', settings });
      }
    }
  }
  return outputs;
}


function evaluateTarget(ctx) {
  const { outputs, phage, bacterium, intF, intName, pWin, pSeqPlus, pSeqMinus, pPlusIndex, pMinusIndex, target, searchType, product, settings } = ctx;
  const rcTarget = rc(target.seq);
  const tests = [
    { strand: '+', qSeq: pSeqPlus, qIndex: pPlusIndex, targetSeq: target.seq },
    { strand: '-', qSeq: pSeqPlus, qIndex: pPlusIndex, targetSeq: rcTarget },
    { strand: '+ vs phageRC', qSeq: pSeqMinus, qIndex: pMinusIndex, targetSeq: target.seq }
  ];
  let best = null;
  for (const t of tests) {
    const core = bestExactCoreFromIndex(t.qSeq, t.qIndex, t.targetSeq, settings.k, settings.minCore);
    if (core && (!best || core.len > best.core.len)) best = { core, ...t };
  }
  if (!best) return;

  const score = scoreHit(best.core, { searchType });
  const hostStart0 = best.strand === '-' ? (target.seq.length - best.core.tEnd) : best.core.tStart;
  const hostLength = Number(target.contigLength || bacterium.sequence?.length || bacterium.sequenceLength || target.seq.length);
  const hostCoord = coordFromWindow(target, hostStart0, best.core.len, hostLength);
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

async function run(payload) {
  const settings = payload.settings;
  postLog('Parsing phage files in full annotation mode…');
  const phages = await readFiles(payload.phageFiles, 'phage');
  postLog('Parsing bacterial inputs. host_index.json files are used directly; GenBank files are parsed in host-optimized mode…');
  const bacteria = await readFiles(payload.bactFiles, 'host');
  postLog(`Parsed ${phages.length} phage record(s) and ${bacteria.length} bacterial record(s).`);

  const all = [];
  for (const p of phages) {
    postLog(`Phage ${p.name}: ${fmt(p.sequence.length)} bp, ${p.features.length} parsed feature(s).`);
    for (const b of bacteria) {
      if (b.format === 'HostIndex') postLog(`Comparing ${p.name} to ${b.name}: indexed host, ${fmt(b.sequenceLength || 0)} total bp.`);
      else postLog(`Comparing ${p.name} to ${b.name}: ${fmt(b.sequence.length)} bp, ${findTRNAs(b).length} parsed tRNA/tmRNA feature(s).`);
      all.push(...analyzePair(p, b, settings));
    }
  }
  const results = dedupeResults(all);
  self.postMessage({ type: 'results', results });
  self.postMessage({ type: 'done' });
}

self.onmessage = (event) => {
  const msg = event.data || {};
  if (msg.type !== 'run') return;
  run(msg).catch(err => self.postMessage({ type: 'error', message: err?.stack || err?.message || String(err) }));
};
