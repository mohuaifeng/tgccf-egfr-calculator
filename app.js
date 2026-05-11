const MODEL = {
  intercept: 4.283243,
  betaLnHeightScr: 0.184304,
  betaNsAge1: 0.046885,
  betaNsAge2: 0.336329,
  betaNsAge3: 0.046227,
  boundaryLow: 0.33,
  boundaryHigh: 20,
  knot1: 5,
  knot2: 10,
  tau2: 0.01007562,
  sigma2: 0.0345,
  wSmall: 0.03,
  wLarge: 0.08,
  ageMin: 0,
  ageMax: 18,
  heightMin: 63,
  heightMax: 189,
  scrUmolMin: 9,
  scrUmolMax: 273.1,
  scrConversionFactor: 88.4
};

const state = {
  base: null,
  calibration: null,
  batchResults: []
};

function formatNumber(value, digits = 2) {
  if (!Number.isFinite(value)) return "--";
  return Number(value).toFixed(digits);
}

function parseNumber(value) {
  const out = Number(String(value).trim());
  return Number.isFinite(out) ? out : NaN;
}

function scrToUmol(value, unit) {
  if (!Number.isFinite(value)) return NaN;
  return unit === "mgdl" ? value * MODEL.scrConversionFactor : value;
}

function scrRangeText(unit) {
  if (unit === "mgdl") {
    return `Eligible range: ${formatNumber(MODEL.scrUmolMin / MODEL.scrConversionFactor, 2)}–${formatNumber(MODEL.scrUmolMax / MODEL.scrConversionFactor, 2)} mg/dL.`;
  }
  return `Eligible range: ${formatNumber(MODEL.scrUmolMin, 1)}–${formatNumber(MODEL.scrUmolMax, 1)} μmol/L.`;
}

function cubicPositivePart(x) {
  return x > 0 ? x * x * x : 0;
}

function bsplineBasis(i, degree, x, knots) {
  if (degree === 0) {
    const isLastInterval = i === knots.length - 2 && x === knots[knots.length - 1];
    return (knots[i] <= x && x < knots[i + 1]) || isLastInterval ? 1 : 0;
  }

  let left = 0;
  const leftDen = knots[i + degree] - knots[i];
  if (leftDen !== 0) {
    left = ((x - knots[i]) / leftDen) * bsplineBasis(i, degree - 1, x, knots);
  }

  let right = 0;
  const rightDen = knots[i + degree + 1] - knots[i + 1];
  if (rightDen !== 0) {
    right = ((knots[i + degree + 1] - x) / rightDen) * bsplineBasis(i + 1, degree - 1, x, knots);
  }

  return left + right;
}

function bsplineDerivative(i, degree, deriv, x, knots) {
  if (deriv === 0) return bsplineBasis(i, degree, x, knots);
  if (degree === 0) return 0;

  let left = 0;
  const leftDen = knots[i + degree] - knots[i];
  if (leftDen !== 0) {
    left = (degree / leftDen) * bsplineDerivative(i, degree - 1, deriv - 1, x, knots);
  }

  let right = 0;
  const rightDen = knots[i + degree + 1] - knots[i + 1];
  if (rightDen !== 0) {
    right = (degree / rightDen) * bsplineDerivative(i + 1, degree - 1, deriv - 1, x, knots);
  }

  return left - right;
}

function rawBsplineRow(age, deriv = 0) {
  const lower = MODEL.boundaryLow;
  const upper = MODEL.boundaryHigh;
  const knots = [
    lower, lower, lower, lower,
    MODEL.knot1, MODEL.knot2,
    upper, upper, upper, upper
  ];
  const degree = 3;
  const nBasis = knots.length - degree - 1;

  function rowAt(x, derivativeOrder) {
    const xx = x === upper ? upper - 1e-12 : x;
    return Array.from({ length: nBasis }, (_, i) => bsplineDerivative(i, degree, derivativeOrder, xx, knots));
  }

  if (deriv !== 0) return rowAt(age, deriv);

  if (age < lower) {
    const basis = rowAt(lower, 0);
    const slope = rowAt(lower, 1);
    return basis.map((v, j) => v + (age - lower) * slope[j]);
  }

  if (age > upper) {
    const basis = rowAt(upper, 0);
    const slope = rowAt(upper, 1);
    return basis.map((v, j) => v + (age - upper) * slope[j]);
  }

  return rowAt(age, 0);
}

function dot(a, b) {
  return a.reduce((sum, v, i) => sum + v * b[i], 0);
}

function norm(a) {
  return Math.sqrt(dot(a, a));
}

function subtract(a, b) {
  return a.map((v, i) => v - b[i]);
}

function scalarMultiply(a, scalar) {
  return a.map(v => v * scalar);
}

function gramSchmidtComplete(initialColumns, dimension) {
  const columns = [];

  for (const col of initialColumns) {
    let v = col.slice();
    for (const q of columns) {
      v = subtract(v, scalarMultiply(q, dot(v, q)));
    }
    const n = norm(v);
    if (n > 1e-12) columns.push(scalarMultiply(v, 1 / n));
  }

  for (let i = 0; i < dimension; i++) {
    let v = Array.from({ length: dimension }, (_, j) => (i === j ? 1 : 0));
    for (const q of columns) {
      v = subtract(v, scalarMultiply(q, dot(v, q)));
    }
    const n = norm(v);
    if (n > 1e-12) columns.push(scalarMultiply(v, 1 / n));
    if (columns.length === dimension) break;
  }

  return columns;
}

function transformToRNsBasis(rawValues) {
  const transform = [
    [-0.324442851, 0.811107103, -0.486664271],
    [0.945905310, 0.278207455, -0.166924469],
    [0.000000006, 0.514495745, 0.857492939]
  ];

  return [0, 1, 2].map(j =>
    rawValues[0] * transform[0][j] +
    rawValues[1] * transform[1][j] +
    rawValues[2] * transform[2][j]
  );
}

function rCompatibleNsAge(age) {
  try {
    const basisRaw = rawBsplineRow(age, 0).slice(1);
    const constLow = rawBsplineRow(MODEL.boundaryLow, 2).slice(1);
    const constHigh = rawBsplineRow(MODEL.boundaryHigh, 2).slice(1);

    const qColumns = gramSchmidtComplete([constLow, constHigh], basisRaw.length);
    const naturalColumns = qColumns.slice(2);
    const values = naturalColumns.map(col => dot(basisRaw, col));

    if (values.length !== 3 || values.some(v => !Number.isFinite(v))) {
      throw new Error("Natural spline basis failed.");
    }

    const rValues = transformToRNsBasis(values);
    if (rValues.some(v => !Number.isFinite(v))) {
      throw new Error("R-compatible natural spline transformation failed.");
    }
    return rValues;
  } catch (error) {
    console.error(error);
    return [NaN, NaN, NaN];
  }
}

function calculateLogBase(age, height, scrUmol) {
  const ns = rCompatibleNsAge(age);
  const lnHeightScr = Math.log(height / scrUmol);
  return (
    MODEL.intercept +
    MODEL.betaLnHeightScr * lnHeightScr +
    MODEL.betaNsAge1 * ns[0] +
    MODEL.betaNsAge2 * ns[1] +
    MODEL.betaNsAge3 * ns[2]
  );
}

function calculateBaseEgfr(age, height, scrUmol) {
  const logBase = calculateLogBase(age, height, scrUmol);
  return {
    logBase,
    egfr: Math.exp(logBase),
    ns: rCompatibleNsAge(age)
  };
}

function validatePatientInputs(age, height, scrUmol) {
  const errors = [];
  const warnings = [];

  if (!Number.isFinite(age) || age < MODEL.ageMin || age > MODEL.ageMax) {
    errors.push(`Age must be between ${MODEL.ageMin} and ${MODEL.ageMax} years.`);
  } else if (age < MODEL.boundaryLow) {
    warnings.push("Age is below 0.33 years, the lower spline boundary used by the model. Interpret cautiously.");
  }

  if (!Number.isFinite(height) || height < MODEL.heightMin || height > MODEL.heightMax) {
    errors.push(`Height must be between ${MODEL.heightMin} and ${MODEL.heightMax} cm.`);
  }

  if (!Number.isFinite(scrUmol) || scrUmol < MODEL.scrUmolMin || scrUmol > MODEL.scrUmolMax) {
    errors.push(`Serum creatinine must be between ${MODEL.scrUmolMin} and ${MODEL.scrUmolMax} μmol/L, equivalent to ${formatNumber(MODEL.scrUmolMin / MODEL.scrConversionFactor, 2)}–${formatNumber(MODEL.scrUmolMax / MODEL.scrConversionFactor, 2)} mg/dL.`);
  }

  return { errors, warnings };
}

function updateBaseResult(result, warnings = []) {
  document.getElementById("baseEgfrValue").textContent = formatNumber(result.egfr, 1);

  let thresholdNote = "Base estimate calculated.";
  if (result.egfr < 85) {
    thresholdNote = "Base eGFR is below 85 mL/min/1.73 m².";
  } else if (result.egfr < 90) {
    thresholdNote = "Base eGFR is below 90 mL/min/1.73 m².";
  } else {
    thresholdNote = "Base eGFR is ≥90 mL/min/1.73 m².";
  }

  if (warnings.length > 0) thresholdNote += " " + warnings.join(" ");
  document.getElementById("baseEgfrNote").textContent = thresholdNote;
}

function showInputMessage(errors, warnings) {
  const box = document.getElementById("inputWarning");
  const messages = [...errors, ...warnings];
  if (messages.length === 0) {
    box.classList.add("hidden");
    box.textContent = "";
    return;
  }
  box.classList.remove("hidden");
  box.textContent = messages.join(" ");
}

function parseCsv(text) {
  const lines = text.split(/\r?\n/).map(line => line.trim()).filter(Boolean);
  if (lines.length < 2) return [];

  const headers = lines[0].split(",").map(h => h.trim());
  return lines.slice(1).map(line => {
    const values = line.split(",").map(v => v.trim());
    const row = {};
    headers.forEach((header, index) => {
      row[header] = values[index];
    });
    return row;
  });
}

function normalizeCalibrationRows(rows, scrUnit) {
  return rows.map((row, index) => {
    const age = parseNumber(row.age ?? row.Age);
    const height = parseNumber(row.height ?? row.Height);
    const scrRaw = parseNumber(row.Scr ?? row.scr ?? row.SCR);
    const scrUmol = scrToUmol(scrRaw, scrUnit);
    const mgfr = parseNumber(row.mGFR ?? row.mgfr ?? row.MGFR ?? row.measured_gfr);
    return { index: index + 1, age, height, scrRaw, scrUmol, scrUnit, mgfr };
  });
}

function validateCalibrationRows(rows) {
  const valid = [];
  const invalid = [];

  for (const row of rows) {
    const check = validatePatientInputs(row.age, row.height, row.scrUmol);
    const reasons = [...check.errors];

    if (!Number.isFinite(row.mgfr) || row.mgfr <= 0) {
      reasons.push("mGFR must be a positive value.");
    }

    if (reasons.length > 0) {
      invalid.push({ index: row.index, reasons });
    } else {
      valid.push(row);
    }
  }

  return { valid, invalid };
}

function estimateCalibration(rows) {
  const n = rows.length;
  if (n < 10) {
    throw new Error("Calibration requires at least 10 valid local samples in this public prototype.");
  }

  const residuals = rows.map(row => {
    const pred = calculateBaseEgfr(row.age, row.height, row.scrUmol);
    return Math.log(row.mgfr) - pred.logBase;
  });

  const wNaive = residuals.reduce((sum, v) => sum + v, 0) / n;
  const useEb = n < 20;
  const lambda = useEb ? MODEL.tau2 / (MODEL.tau2 + MODEL.sigma2 / n) : 1;
  const wUsed = useEb ? lambda * wNaive : wNaive;
  const factor = Math.exp(wUsed);

  return {
    n,
    residuals,
    wNaive,
    lambda,
    useEb,
    wUsed,
    factor,
    route: useEb ? "EB intercept-only" : "Naive intercept-only"
  };
}

function recommendationFromCalibration(calibration) {
  const absW = Math.abs(calibration.wNaive);
  if (absW < MODEL.wSmall) {
    return {
      type: "bad",
      title: "Calibration not recommended",
      text: "The estimated center-level intercept shift is small. The base eGFR estimate should generally be reported without center calibration."
    };
  }

  if (absW < MODEL.wLarge) {
    return {
      type: "review",
      title: "Review before routine calibration",
      text: "The estimated center-level intercept shift is moderate. Additional local samples or independent validation are recommended before routine use."
    };
  }

  return {
    type: "good",
    title: "Center calibration may be considered",
    text: "The local calibration sample suggests a meaningful center-level intercept shift. Intercept-only center calibration may be considered, subject to local validation."
  };
}

function resetCalibrationOutputs(message = "A recommendation requires at least 10 valid local mGFR calibration samples.") {
  state.calibration = null;
  document.getElementById("wNaiveValue").textContent = "--";
  document.getElementById("wUsedValue").textContent = "--";
  document.getElementById("calibrationFactorValue").textContent = "--";
  document.getElementById("calibrationRouteValue").textContent = "--";
  document.getElementById("calibratedEgfrValue").textContent = "--";
  document.getElementById("calibratedEgfrNote").textContent = "Upload at least 10 valid local mGFR samples to enable center calibration.";
  const card = document.getElementById("recommendationCard");
  card.className = "recommendation-card error";
  document.getElementById("recommendationTitle").textContent = "Calibration unavailable";
  document.getElementById("recommendationText").textContent = message;
}

function updateCalibrationDisplay(calibration) {
  state.calibration = calibration;

  document.getElementById("wNaiveValue").textContent = formatNumber(calibration.wNaive, 4);
  document.getElementById("wUsedValue").textContent = formatNumber(calibration.wUsed, 4);
  document.getElementById("calibrationFactorValue").textContent = formatNumber(calibration.factor, 3);
  document.getElementById("calibrationRouteValue").textContent = calibration.route;

  const rec = recommendationFromCalibration(calibration);
  const card = document.getElementById("recommendationCard");
  card.className = `recommendation-card ${rec.type}`;
  document.getElementById("recommendationTitle").textContent = rec.title;
  document.getElementById("recommendationText").textContent = rec.text;

  if (state.base) {
    const calibrated = state.base.egfr * calibration.factor;
    document.getElementById("calibratedEgfrValue").textContent = formatNumber(calibrated, 1);
    document.getElementById("calibratedEgfrNote").textContent = `Applied ${calibration.route} using local n = ${calibration.n}.`;
  }
}

function renderCsvPreview(rows) {
  const container = document.getElementById("csvPreview");
  if (!rows.length) {
    container.classList.add("hidden");
    container.innerHTML = "";
    return;
  }

  const previewRows = rows.slice(0, 5);
  const html = `
    <table>
      <thead>
        <tr><th>#</th><th>Age</th><th>Height</th><th>Scr, μmol/L</th><th>mGFR</th></tr>
      </thead>
      <tbody>
        ${previewRows.map(row => `
          <tr>
            <td>${row.index}</td>
            <td>${formatNumber(row.age, 2)}</td>
            <td>${formatNumber(row.height, 1)}</td>
            <td>${formatNumber(row.scrUmol, 2)}</td>
            <td>${formatNumber(row.mgfr, 1)}</td>
          </tr>
        `).join("")}
      </tbody>
    </table>
  `;
  container.innerHTML = html;
  container.classList.remove("hidden");
}

function renderExcludedRows(invalid) {
  const box = document.getElementById("excludedRowsBox");
  if (!invalid.length) {
    box.classList.add("hidden");
    box.innerHTML = "";
    return;
  }

  const maxShow = 8;
  const shown = invalid.slice(0, maxShow);
  const extra = invalid.length > maxShow ? `<p>${invalid.length - maxShow} additional excluded row(s) not shown.</p>` : "";
  box.innerHTML = `
    <strong>Excluded row details</strong>
    <ul>
      ${shown.map(item => `<li>Row ${item.index}: ${item.reasons.join(" ")}</li>`).join("")}
    </ul>
    ${extra}
  `;
  box.classList.remove("hidden");
}

function updateRowCounts(uploaded, valid, invalid) {
  document.getElementById("uploadedRowsValue").textContent = Number.isFinite(uploaded) ? uploaded : "--";
  document.getElementById("validRowsValue").textContent = Number.isFinite(valid) ? valid : "--";
  document.getElementById("excludedRowsValue").textContent = Number.isFinite(invalid) ? invalid : "--";
}


function normalizeBatchRows(rows, scrUnit) {
  return rows.map((row, index) => {
    const id = row.id ?? row.ID ?? row.Id ?? String(index + 1);
    const age = parseNumber(row.age ?? row.Age);
    const height = parseNumber(row.height ?? row.Height);
    const scrRaw = parseNumber(row.Scr ?? row.scr ?? row.SCR);
    const scrUmol = scrToUmol(scrRaw, scrUnit);
    return { index: index + 1, id, age, height, scrRaw, scrUmol, scrUnit };
  });
}

function validateBatchRows(rows) {
  const valid = [];
  const invalid = [];

  for (const row of rows) {
    const check = validatePatientInputs(row.age, row.height, row.scrUmol);
    const reasons = [...check.errors];
    if (reasons.length > 0) {
      invalid.push({ index: row.index, reasons });
    } else {
      valid.push(row);
    }
  }

  return { valid, invalid };
}

function calculateBatchBaseEgfr(rows) {
  return rows.map(row => {
    const pred = calculateBaseEgfr(row.age, row.height, row.scrUmol);
    return {
      id: row.id,
      age: row.age,
      height: row.height,
      Scr_umol_L: row.scrUmol,
      log_base_eGFR: pred.logBase,
      base_eGFR: pred.egfr
    };
  });
}

function renderBatchExcludedRows(invalid) {
  const box = document.getElementById("batchExcludedRowsBox");
  if (!invalid.length) {
    box.classList.add("hidden");
    box.innerHTML = "";
    return;
  }

  const maxShow = 8;
  const shown = invalid.slice(0, maxShow);
  const extra = invalid.length > maxShow ? `<p>${invalid.length - maxShow} additional excluded row(s) not shown.</p>` : "";
  box.innerHTML = `
    <strong>Excluded row details</strong>
    <ul>
      ${shown.map(item => `<li>Row ${item.index}: ${item.reasons.join(" ")}</li>`).join("")}
    </ul>
    ${extra}
  `;
  box.classList.remove("hidden");
}

function renderBatchPreview(results) {
  const container = document.getElementById("batchPreview");
  const downloadButton = document.getElementById("downloadBatchResults");

  if (!results.length) {
    container.classList.add("hidden");
    container.innerHTML = "";
    downloadButton.classList.add("hidden");
    return;
  }

  const previewRows = results.slice(0, 10);
  const html = `
    <table>
      <thead>
        <tr><th>ID</th><th>Age</th><th>Height</th><th>Scr, μmol/L</th><th>Base eGFR</th></tr>
      </thead>
      <tbody>
        ${previewRows.map(row => `
          <tr>
            <td>${String(row.id)}</td>
            <td>${formatNumber(row.age, 2)}</td>
            <td>${formatNumber(row.height, 1)}</td>
            <td>${formatNumber(row.Scr_umol_L, 2)}</td>
            <td>${formatNumber(row.base_eGFR, 1)}</td>
          </tr>
        `).join("")}
      </tbody>
    </table>
  `;
  container.innerHTML = html;
  container.classList.remove("hidden");
  downloadButton.classList.remove("hidden");
}

function csvEscape(value) {
  const text = String(value ?? "");
  if (/[",\n\r]/.test(text)) return `"${text.replace(/"/g, '""')}"`;
  return text;
}

function resultsToCsv(results) {
  const headers = ["id", "age", "height", "Scr_umol_L", "log_base_eGFR", "base_eGFR"];
  const lines = [headers.join(",")];
  for (const row of results) {
    lines.push(headers.map(header => csvEscape(row[header])).join(","));
  }
  return lines.join("\n");
}

function downloadBatchResults() {
  if (!state.batchResults.length) return;
  const blob = new Blob([resultsToCsv(state.batchResults)], { type: "text/csv;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = "tgccf_base_egfr_results.csv";
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
}

function resetBatchOutputs(message = "No batch file loaded.") {
  state.batchResults = [];
  document.getElementById("batchStatus").textContent = message;
  document.getElementById("batchPreview").classList.add("hidden");
  document.getElementById("batchPreview").innerHTML = "";
  document.getElementById("batchExcludedRowsBox").classList.add("hidden");
  document.getElementById("batchExcludedRowsBox").innerHTML = "";
  document.getElementById("downloadBatchResults").classList.add("hidden");
}

function handleBatchCsvUpload(event) {
  const file = event.target.files[0];
  if (!file) {
    resetBatchOutputs("No batch file loaded.");
    return;
  }

  const reader = new FileReader();
  reader.onload = function (e) {
    try {
      const rawRows = parseCsv(String(e.target.result));
      const scrUnit = document.getElementById("batchScrUnit").value;
      const rows = normalizeBatchRows(rawRows, scrUnit);
      const { valid, invalid } = validateBatchRows(rows);

      renderBatchExcludedRows(invalid);

      if (valid.length === 0) {
        state.batchResults = [];
        renderBatchPreview([]);
        document.getElementById("batchStatus").textContent = `Uploaded rows: ${rows.length}. Valid rows: 0. Excluded rows: ${invalid.length}. No valid rows were available for base eGFR calculation.`;
        return;
      }

      const results = calculateBatchBaseEgfr(valid);
      state.batchResults = results;
      renderBatchPreview(results);
      document.getElementById("batchStatus").textContent = `Uploaded rows: ${rows.length}. Valid rows: ${valid.length}. Excluded rows: ${invalid.length}. Base eGFR was calculated for valid rows.`;
    } catch (error) {
      resetBatchOutputs(error.message);
    }
  };
  reader.readAsText(file);
}

function calculateCurrentPatient() {
  const age = parseNumber(document.getElementById("ageInput").value);
  const height = parseNumber(document.getElementById("heightInput").value);
  const scrRaw = parseNumber(document.getElementById("scrInput").value);
  const scrUnit = document.getElementById("scrUnitInput").value;
  const scrUmol = scrToUmol(scrRaw, scrUnit);

  const validation = validatePatientInputs(age, height, scrUmol);
  showInputMessage(validation.errors, validation.warnings);

  if (validation.errors.length > 0) return;

  const result = calculateBaseEgfr(age, height, scrUmol);
  state.base = result;
  updateBaseResult(result, validation.warnings);

  if (state.calibration) {
    const calibrated = result.egfr * state.calibration.factor;
    document.getElementById("calibratedEgfrValue").textContent = formatNumber(calibrated, 1);
    document.getElementById("calibratedEgfrNote").textContent = `Applied ${state.calibration.route} using local n = ${state.calibration.n}.`;
  }
}

function handleCsvUpload(event) {
  const file = event.target.files[0];
  const status = document.getElementById("csvStatus");

  if (!file) {
    status.textContent = "No local calibration file loaded.";
    updateRowCounts(NaN, NaN, NaN);
    return;
  }

  const reader = new FileReader();
  reader.onload = function (e) {
    try {
      const rawRows = parseCsv(String(e.target.result));
      const scrUnit = document.getElementById("csvScrUnit").value;
      const rows = normalizeCalibrationRows(rawRows, scrUnit);
      const { valid, invalid } = validateCalibrationRows(rows);

      updateRowCounts(rows.length, valid.length, invalid.length);
      renderCsvPreview(valid);
      renderExcludedRows(invalid);

      if (valid.length === 0) {
        throw new Error("No valid calibration rows were found. Required columns are age, height, Scr, and mGFR.");
      }

      status.textContent = `Uploaded rows: ${rows.length}. Valid rows: ${valid.length}. Excluded rows: ${invalid.length}.`;

      if (valid.length < 10) {
        resetCalibrationOutputs("Calibration was not applied because this public prototype requires at least 10 valid local samples.");
        status.textContent += " Calibration not applied: at least 10 valid samples are required.";
        return;
      }

      const calibration = estimateCalibration(valid);
      updateCalibrationDisplay(calibration);
      status.textContent += ` Calibration applied using ${calibration.route}.`;
    } catch (error) {
      status.textContent = error.message;
      document.getElementById("csvPreview").classList.add("hidden");
      resetCalibrationOutputs(error.message);
    }
  };
  reader.readAsText(file);
}

function updateScrRangeHint() {
  const unit = document.getElementById("scrUnitInput").value;
  document.getElementById("scrRangeHint").textContent = scrRangeText(unit);
}

document.getElementById("calculateButton").addEventListener("click", calculateCurrentPatient);
document.getElementById("csvInput").addEventListener("change", handleCsvUpload);
document.getElementById("batchCsvInput").addEventListener("change", handleBatchCsvUpload);
document.getElementById("downloadBatchResults").addEventListener("click", downloadBatchResults);
document.getElementById("scrUnitInput").addEventListener("change", updateScrRangeHint);

document.getElementById("csvScrUnit").addEventListener("change", () => {
  const input = document.getElementById("csvInput");
  if (input.files && input.files.length > 0) {
    handleCsvUpload({ target: input });
  }
});

document.getElementById("batchScrUnit").addEventListener("change", () => {
  const input = document.getElementById("batchCsvInput");
  if (input.files && input.files.length > 0) {
    handleBatchCsvUpload({ target: input });
  }
});

["ageInput", "heightInput", "scrInput"].forEach(id => {
  document.getElementById(id).addEventListener("keydown", event => {
    if (event.key === "Enter") calculateCurrentPatient();
  });
});

updateScrRangeHint();

window.TGCCF = {
  calculateBaseEgfr,
  rCompatibleNsAge,
  estimateCalibration,
  calculateBatchBaseEgfr,
  scrToUmol,
  MODEL
};
