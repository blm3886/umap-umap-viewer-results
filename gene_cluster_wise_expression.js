import * as THREE from "three";
import Plotly from 'plotly.js-dist-min';
import Papa from "papaparse";
import { confusionData } from "./per_macrogene_confusion_matrix_with_og.js";
import { geneList } from "./geneList.js";

// const geneList = ["AT1G01580", "AT1G04310", "AT1G04360"];
const sidebar = document.getElementById("gene-sidebar");
const plotDiv = document.getElementById("gene-plot");
const toggleBtn = document.getElementById("toggle-gene-list");
const geneListDiv = document.getElementById("gene-list");

// let geneToFile = {
//   "AT1G01580": "AT1G01580_MG_1807_gene_cluster_contrib.json",
//   "AT1G04310": "AT1G04310_MG_1372_gene_cluster_contrib.json",
//   "AT1G04360": "AT1G04360_MG_0014_gene_cluster_contrib.json"
// }

let geneToMG = {
"AT1G01580": "MG_1807",
"AT1G04310" : "MG_1372",
"AT1G04360" : "MG_0014",
}


Object.keys(geneList).forEach(gene => {
  const record = geneList[gene][0];
  const mgId = record.macrogenes;
  const fileName = `${gene}_${mgId}_gene_cluster_contrib.json`; 
  const gene_name = record["Gene Name"]

  const btn = document.createElement("button");
  btn.textContent = `${gene} – ${gene_name}`;
  btn.onclick = () => {
    console.log("Look for file", fileName)
    // const fileName = geneToFile[gene];
    if (!fileName) {
      console.error("No file found for gene:", gene);
      return;
    }
    fetch(`/umap-umap-viewer-site/Expression_per_Cluster/${fileName}`)
      .then(r => r.json())
      .then(data => {
        renderHeatmap(gene, data);
        renderConfusionMatrices(mgId);
    });

    const umapFile = `${mgId}_umap3d.json`;
    console.log("UMAP FILE TO BE READ",umapFile)
    const testumapFIle = `MG_0450_umap3d.json`;
    fetch(`/umap-umap-viewer-site/perMacrogene_UMAPP_CSV/${testumapFIle}`)
        .then(r => r.json())
        .then(umapData => {renderUMAPP(mgId, umapData), renderClusterFrequenciesAccordion(umapData)});
    };

  geneListDiv.appendChild(btn);
});


function renderHeatmap(gene, data) {
  const genes = Object.keys(data);
  const clusters = Object.keys(data[genes[0]]);
  const matrix = genes.map(g => clusters.map(c => data[g][c]));

  const trace = {
    z: matrix,
    x: clusters,
    y: genes,
    type: 'heatmap',
    colorscale: 'Viridis',
    hovertemplate:
      "Cluster: %{x}<br>" +
      "Gene: %{y}<br>" +
      "Fraction: %{z:.3f}<extra></extra>"
  };

  const tickvals = genes;
  const ticktext = genes.map(g => {
  const shortId = g.split("_")[1];
  return shortId === gene
    ? `<span style="color:red;">${g}</span>`
    : g;
});

  const layout = {
    margin: { t: 40, l: 150, b: 80 },
    xaxis: {
      title: "Clusters",
      tickmode: "linear",
      dtick: 1,
      automargin: true
    },
    yaxis: {
      title: "Genes",
      automargin: true,
      tickmode: "array",
      tickvals: tickvals,
      ticktext: ticktext
    },
    width: 800,
    height: genes.length * 20
  };

  Plotly.newPlot('gene-plot', [trace], layout);
}

function renderConfusionMatrices(mgId) {
  const container = document.getElementById("confusion-container");
  container.innerHTML = ""; // clear old

  if (!confusionData[mgId]) {
    container.innerHTML = `<p>No data for ${mgId}</p>`;
    return;
  }

  confusionData[mgId].forEach(row => {
    const { orthogroup, TP, FP, FN, TN, precision, recall } = row;

    // Create table
    const table = document.createElement("table");
    table.style.borderCollapse = "collapse";
    table.style.marginBottom = "20px";
    table.style.width = "100%";
    table.style.textAlign = "center";
    table.style.color = "white";
    table.style.border = "1px solid #444";
    table.style.background = "#111";

    table.innerHTML = `
    <caption style="caption-side: top; font-weight:bold; margin-bottom:5px;">
        ${mgId} – ${orthogroup} 
        (P=${(+precision).toFixed(2)}, R=${(+recall).toFixed(2)})
    </caption>
    <tr>
        <td></td>
        <td style="padding:4px;border:1px solid #444;">Macrogene +</td>
        <td style="padding:4px;border:1px solid #444;">Macrogene -</td>
    </tr>
    <tr>
        <td style="padding:4px;border:1px solid #444;">Orthogroup +</td>
        <td style="padding:4px;border:1px solid #444;">${TP}</td>
        <td style="padding:4px;border:1px solid #444;">${FP}</td>
       
    </tr>
    <tr>
        <td style="padding:4px;border:1px solid #444;">Orthogroup -</td>
        <td style="padding:4px;border:1px solid #444;">${TN}</td>
        <td style="padding:4px;border:1px solid #444;">${FN}</td>
    </tr>
    `;
    container.appendChild(table);
  });
}

function renderUMAPP(macrogene_id, umapData) {
  const container = document.getElementById("umap-container");
  const plotDiv = document.getElementById("umap-plot");

  const skipCols = new Set([
    "UMAP1", "UMAP2", "UMAP3",
    "cell_id", "species", "cluster", "cell_type"
  ]);

  const allGenes = Object.keys(umapData[0]).filter(k => !skipCols.has(k));

  // max expression per cell
  const expr = umapData.map(row =>
    Math.max(...allGenes.map(g => +row[g]))
  );

  const x = umapData.map(r => +r.UMAP1);
  const y = umapData.map(r => +r.UMAP2);
  const z = umapData.map(r => +r.UMAP3);
  const species = umapData.map(r => r.species);

  const speciesColors = {
    Arabidopsis: "red",
    Pennycress: "green",
    Soybean: "blue"
  };

  // find smallest non-zero
  const minNonZero = Math.min(...expr.filter(v => v > 0));
  const maxExpr = Math.max(...expr);

  function makeTrace(threshold) {
    const colors = expr.map((v, i) =>
      v > 0 && v >= threshold ? (speciesColors[species[i]] || "black")
                              : "rgba(80, 80, 80, 0.3)"
    );

    return {
      x, y, z,
      mode: "markers",
      type: "scatter3d",
      marker: { size: 3, color: colors },
    };
  }

  const layout = {
    title: `3D UMAP – ${macrogene_id}`,
    scene: { xaxis: {title:"UMAP1"}, yaxis: {title:"UMAP2"}, zaxis: {title:"UMAP3"} },
    height: 600, width: 800
  };

  // start with minNonZero threshold
  Plotly.newPlot(plotDiv, [makeTrace(minNonZero)], layout);

  // === Controls ===
  let controls = document.getElementById("umap-controls");
  if (!controls) {
    controls = document.createElement("div");
    controls.id = "umap-controls";
    controls.style.marginBottom = "10px";
    container.insertBefore(controls, plotDiv);
  }
  controls.innerHTML = "";

  const thresholdLabel = document.createElement("div");
  thresholdLabel.style.color = "white";
  thresholdLabel.textContent = `Threshold: ${minNonZero.toFixed(2)}`;

  const slider = document.createElement("input");
  slider.type = "range";
  slider.min = minNonZero;
  slider.max = maxExpr;
  slider.step = 0.01;
  slider.value = minNonZero;
  slider.style.width = "80%";

  const geneBox = document.createElement("div");
  geneBox.id = "gene-box";
  geneBox.style.marginTop = "10px";
  geneBox.style.color = "white";
  geneBox.style.fontSize = "12px";
  geneBox.style.maxHeight = "200px";
  geneBox.style.overflowY = "auto";
  geneBox.style.border = "1px solid #555";
  geneBox.style.padding = "5px";

  function updateContributingGenes(threshold) {
    const contributing = new Set();
    umapData.forEach(row => {
      allGenes.forEach(g => {
        if (+row[g] >= threshold) contributing.add(g);
      });
    });
    geneBox.innerHTML = "<b>Contributing genes:</b><br>" + Array.from(contributing).join(", ");
  }

  slider.oninput = e => {
    const t = +e.target.value;
    thresholdLabel.textContent = `Threshold: ${t.toFixed(2)}`;
    Plotly.react(plotDiv, [makeTrace(t)], layout);
    updateContributingGenes(t);
  };

  // init
  controls.appendChild(thresholdLabel);
  controls.appendChild(slider);
  controls.appendChild(geneBox);
  updateContributingGenes(minNonZero);
}


function renderClusterFrequenciesAccordion(umapData, targetDivId = "cluster-summary") {
  // Group counts: cluster → species → cell_type counts
  const counts = {};
  umapData.forEach(r => {
    if (!counts[r.cluster]) counts[r.cluster] = {};
    if (!counts[r.cluster][r.species]) counts[r.cluster][r.species] = {};
    counts[r.cluster][r.species][r.cell_type] =
      (counts[r.cluster][r.species][r.cell_type] || 0) + 1;
  });

  const container = document.getElementById(targetDivId);
  container.innerHTML = "";

  Object.entries(counts).forEach(([cluster, speciesMap]) => {
    // === Collapsible header ===
    const header = document.createElement("button");
    header.textContent = `Cluster ${cluster}`;
    header.style.width = "100%";
    header.style.textAlign = "left";
    header.style.padding = "10px";
    header.style.margin = "5px 0";
    header.style.border = "1px solid #555";
    header.style.background = "#222";
    header.style.color = "white";
    header.style.cursor = "pointer";
    header.style.fontWeight = "bold";

    // === Content div ===
    const content = document.createElement("div");
    content.style.display = "none";
    content.style.padding = "15px";
    content.style.background = "#111";
    content.style.border = "1px solid #333";
    content.style.marginBottom = "20px";

    // Per-species plots
    ["Arabidopsis", "Pennycress", "Soybean"].forEach(species => {
      // Species title
      const speciesTitle = document.createElement("h3");
      speciesTitle.textContent = species;
      speciesTitle.style.color = "white";
      speciesTitle.style.margin = "10px 0";
      speciesTitle.style.fontSize = "16px";
      content.appendChild(speciesTitle);

      // Plot div
      const plotDiv = document.createElement("div");
      plotDiv.style.width = "100%";
      plotDiv.style.height = "300px";
      plotDiv.style.marginBottom = "200px"; // space between plots
      content.appendChild(plotDiv);

      if (speciesMap[species]) {
        const ctMap = speciesMap[species];
        const total = Object.values(ctMap).reduce((a, b) => a + b, 0);
        const cellTypes = Object.keys(ctMap);
        const values = cellTypes.map(ct => (ctMap[ct] / total) * 100);

        const trace = {
          x: cellTypes,
          y: values,
          type: "bar",
          marker: { color: "rgba(100,149,237,0.8)" },
          text: values.map(v => v.toFixed(1) + "%"),
          textposition: "auto",
          hovertemplate:
            `Species: ${species}<br>Cluster: ${cluster}<br>Cell type: %{x}<br>Percentage: %{y:.1f}%<extra></extra>`
        };

        const layout = {
          xaxis: { tickangle: -45, automargin: true },
          yaxis: { title: "Percentage", range: [0, 100], fixedrange: true },
          margin: { t: 10, b: 80, l: 60, r: 20 },
        };

        Plotly.newPlot(plotDiv, [trace], layout, { displayModeBar: false });
      } else {
        plotDiv.innerHTML = `<p style="color:gray;">No data for ${species}</p>`;
      }
    });

    // Toggle behavior
    header.onclick = () => {
      content.style.display = (content.style.display === "none" ? "block" : "none");
    };

    // Append to container
    container.appendChild(header);
    container.appendChild(content);
  });
}





