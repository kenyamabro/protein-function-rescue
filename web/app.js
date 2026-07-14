/* Protein Function Rescue — static viewer.
 * Reads window.RESCUE_DATA (from data.js), renders a sortable/filterable table,
 * and shows each AlphaFold structure in 3D on selection via 3Dmol.js.
 */
(function () {
  "use strict";

  const data = window.RESCUE_DATA || { candidates: [], organism: "(no data)" };
  const candidates = data.candidates || [];

  // Flatten the nested records into a table-friendly shape.
  const rows = candidates.map((c) => ({
    accession: c.accession,
    name: c.name,
    matchDescription: c.best_structural_match.description,
    tm: c.best_structural_match.tm_score,
    plddt: c.mean_plddt,
    drug: c.pocket.top_druggability,
    composite: c.scores.composite,
    _raw: c,
  }));

  const tbody = document.querySelector("#candidates tbody");
  const searchEl = document.getElementById("search");
  const countEl = document.getElementById("count");
  const metaEl = document.getElementById("meta");
  const detailEl = document.getElementById("detail");

  metaEl.textContent = `${data.organism} · ${data.n_candidates ?? rows.length} rescued candidates`;

  const bannerEl = document.getElementById("sample-banner");
  if (data.sample && bannerEl) {
    bannerEl.textContent = "⚠ " + (data.sample_note || "Sample data — matches are illustrative.");
    bannerEl.hidden = false;
  }

  let sortKey = "composite";
  let sortAsc = false;
  let filter = "";
  let selectedAcc = null;

  let viewer = null;
  function initViewer() {
    const el = document.getElementById("viewer");
    if (!window.$3Dmol) return null;
    if (!viewer) {
      viewer = $3Dmol.createViewer(el, { backgroundColor: "0x0b0f18" });
    }
    return viewer;
  }

  function fmt(v, digits) {
    if (v === null || v === undefined || v === "") return "—";
    if (typeof v === "number") return v.toFixed(digits);
    return v;
  }

  function visibleRows() {
    let r = rows;
    if (filter) {
      const f = filter.toLowerCase();
      r = r.filter(
        (row) =>
          row.accession.toLowerCase().includes(f) ||
          (row.name || "").toLowerCase().includes(f) ||
          (row.matchDescription || "").toLowerCase().includes(f)
      );
    }
    const dir = sortAsc ? 1 : -1;
    r = r.slice().sort((a, b) => {
      let av = a[sortKey];
      let bv = b[sortKey];
      if (av === null || av === undefined) av = -Infinity;
      if (bv === null || bv === undefined) bv = -Infinity;
      if (typeof av === "string" || typeof bv === "string") {
        return dir * String(av).localeCompare(String(bv));
      }
      return dir * (av - bv);
    });
    return r;
  }

  function renderTable() {
    const r = visibleRows();
    countEl.textContent = `${r.length} shown`;
    tbody.innerHTML = "";
    const frag = document.createDocumentFragment();

    for (const row of r) {
      const tr = document.createElement("tr");
      if (row.accession === selectedAcc) tr.classList.add("selected");
      tr.innerHTML = `
        <td><span class="acc">${row.accession}</span></td>
        <td>${escapeHtml(row.name)}<span class="badge">unannotated</span></td>
        <td>${escapeHtml(row.matchDescription)}</td>
        <td class="num">${fmt(row.tm, 3)}</td>
        <td class="num">${fmt(row.plddt, 1)}</td>
        <td class="num">${fmt(row.drug, 2)}</td>
        <td class="num"><strong>${fmt(row.composite, 3)}</strong></td>`;
      tr.addEventListener("click", () => select(row.accession));
      frag.appendChild(tr);
    }
    tbody.appendChild(frag);

    document.querySelectorAll("thead th").forEach((th) => {
      th.classList.toggle("sorted", th.dataset.key === sortKey);
      th.classList.toggle("asc", th.dataset.key === sortKey && sortAsc);
    });
  }

  function select(acc) {
    selectedAcc = acc;
    renderTable();
    const c = candidates.find((x) => x.accession === acc);
    if (c) {
      renderDetail(c);
      loadStructure(c);
    }
  }

  function scoreMeter(label, value) {
    const pct = value === null || value === undefined ? 0 : Math.round(value * 100);
    const shown = value === null || value === undefined ? "—" : value.toFixed(2);
    return `<div class="lbl">${label}</div>
      <div><div class="meter"><span style="width:${pct}%"></span></div></div>
      <div class="lbl">&nbsp;</div><div class="num" style="text-align:right">${shown}</div>`;
  }

  function renderDetail(c) {
    const m = c.best_structural_match;
    const p = c.pocket;
    const s = c.scores;
    // fpocket reports a trained "druggability"; the geometry backend a "cavity score".
    const pocketLabel = p.method === "geometry" ? "cavity score" : "druggability";
    const evaluePart = m.evalue != null ? ` · E-value ${m.evalue}` : "";
    detailEl.innerHTML = `
      <h2>${escapeHtml(c.name)} <span class="acc">${c.accession}</span></h2>
      <div class="row"><span class="k">Length / mean pLDDT</span>
        <span class="v">${c.length} aa · ${c.mean_plddt}</span></div>
      <div class="row"><span class="k">Best structural match</span>
        <span class="v">${escapeHtml(m.description)} <span class="pill">${escapeHtml(m.target)}</span></span></div>
      <div class="row"><span class="k">Match confidence</span>
        <span class="v">TM-score ${m.tm_score}${
          m.tm_score_query != null && m.tm_score_reference != null
            ? ` (query ${m.tm_score_query}; reference ${m.tm_score_reference})`
            : ""
        }${evaluePart}</span></div>
      <div class="row"><span class="k">Binding pocket</span>
        <span class="v">${p.n_pockets} pocket(s)${
      p.top_druggability != null ? ` · top ${pocketLabel} ${p.top_druggability}` : ""
    }</span></div>
      <div class="scoregrid">
        ${scoreMeter("Structural similarity", s.structural_similarity)}
        ${scoreMeter("Model confidence", s.model_confidence)}
        ${scoreMeter("Pocket quality", s.pocket_quality)}
        <div class="lbl"><strong>Composite</strong></div>
        <div><div class="meter"><span style="width:${Math.round(
          (s.composite || 0) * 100
        )}%"></span></div></div>
        <div class="lbl">&nbsp;</div><div class="num" style="text-align:right"><strong>${
          s.composite
        }</strong></div>
      </div>
      <div class="row" style="margin-top:10px">
        <a href="${c.afdb_entry_url}" target="_blank" rel="noopener">Open in AlphaFold DB ↗</a>
      </div>`;
  }

  function loadStructure(c) {
    const v = initViewer();
    const el = document.getElementById("viewer");
    if (!v) {
      showViewerError("3Dmol.js failed to load (offline?).");
      return;
    }
    v.clear();
    // Try PDB first, fall back to CIF.
    fetch(c.structure_pdb_url)
      .then((r) => {
        if (!r.ok) throw new Error(r.status);
        return r.text();
      })
      .then((pdb) => {
        v.addModel(pdb, "pdb");
        colourByConfidence(v);
        v.zoomTo();
        v.render();
      })
      .catch(() => {
        showViewerError(
          `Could not fetch the structure directly.<br/>` +
            `<a href="${c.afdb_entry_url}" target="_blank" rel="noopener">View it on AlphaFold DB ↗</a>` +
            `<br/><small>(Serving this folder over http, e.g. <code>python -m http.server</code>, avoids browser CORS limits.)</small>`
        );
      });
  }

  // AlphaFold convention: colour by pLDDT (stored in the B-factor column).
  function colourByConfidence(v) {
    v.setStyle({}, {
      cartoon: {
        colorfunc: (atom) => {
          const b = atom.b;
          if (b > 90) return "0x0053d6"; // very high
          if (b > 70) return "0x65cbf3"; // confident
          if (b > 50) return "0xffdb13"; // low
          return "0xff7d45"; // very low
        },
      },
    });
  }

  function showViewerError(html) {
    const el = document.getElementById("viewer");
    let err = el.querySelector(".viewer-error");
    if (!err) {
      err = document.createElement("div");
      err.className = "viewer-error";
      el.appendChild(err);
    }
    err.innerHTML = html;
  }

  function escapeHtml(s) {
    return String(s == null ? "" : s).replace(/[&<>"']/g, (ch) => ({
      "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
    }[ch]));
  }

  // Wire up sorting + search.
  document.querySelectorAll("thead th").forEach((th) => {
    th.addEventListener("click", () => {
      const key = th.dataset.key;
      if (key === sortKey) sortAsc = !sortAsc;
      else {
        sortKey = key;
        sortAsc = false;
      }
      renderTable();
    });
  });
  searchEl.addEventListener("input", (e) => {
    filter = e.target.value;
    renderTable();
  });

  renderTable();
  if (rows.length) select(visibleRows()[0].accession);
})();
