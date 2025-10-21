<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>CanEpiRisk R Package</title>
  <style>
    body {
      font-family: Arial, sans-serif;
      line-height: 1.6;
      max-width: 800px;
      margin: auto;
      padding: 2em;
    }
    code {
      background-color: #f4f4f4;
      padding: 2px 4px;
      border-radius: 4px;
      font-family: Consolas, monospace;
    }
    pre {
      background-color: #f4f4f4;
      padding: 1em;
      border-radius: 5px;
      overflow-x: auto;
    }
    h1, h2, h3 {
      color: #2c3e50;
    }
  </style>
</head>
<body>

<h1>📦 CanEpiRisk</h1>

<p><strong>CanEpiRisk</strong> is an R package for estimating lifetime attributable cancer risk due to radiation exposure. It combines epidemiological models with demographic data to support individual and population-level risk assessments.</p>

<h2>🚀 Installation</h2>

<p>You can install the development version of CanEpiRisk from GitHub using:</p>

<pre><code># install.packages("devtools")
devtools::install_github("yourusername/CanEpiRisk")
</code></pre>

<h2>📘 Example Usage</h2>

<pre><code>library(CanEpiRisk)

# Prepare input data
ex_data &lt;- data.frame(
  sex = "female",
  birth = 1965,
  exposure = 1985,
  site = "thyroid",
  exposure_rate = "acute",
  dosedist = "fixedvalue",
  dose1 = 20,
  dose2 = NA,
  dose3 = NA
)

# Load built-in life table and cancer incidence data
data("life_table")
data("incidence_table")

# Calculate lifetime attributable risk
result &lt;- LAR(
  data = ex_data,
  basedata = list(life_table, incidence_table),
  sim = 500,
  seed = 123,
  current = 2025,
  ci = 0.95,
  DDREF = TRUE
)

# View results
print(result)
plot(result)
</code></pre>

<h2>🔧 Key Functions</h2>
<ul>
  <li><code>LAR()</code>: Calculates lifetime attributable risk for individuals</li>
  <li><code>LAR_group()</code>: Estimates risk across population subgroups</li>
  <li><code>life_table</code>, <code>incidence_table</code>: Built-in demographic datasets</li>
</ul>

<h2>📄 License</h2>
<p>This package is released under the MIT License.</p>

<h2>🙋‍♀️ Contributing</h2>
<p>Contributions are welcome! Feel free to open issues or submit pull requests to improve the package.</p>

</body>
</html>
