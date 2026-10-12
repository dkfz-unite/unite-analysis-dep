# Differential Protein Expression Analysis Service (Limma)

Shared architecture, development workflows, and agent instructions:
[UNITE project documentation](https://github.com/dkfz-unite/.github/blob/main/README.md).

## General
R [limma](https://bioconductor.org/packages/limma) based differential protein expression analysis application wrapped with web API.

The R script ([src/run/run.R](src/run/run.R)) is run by the [commands](https://github.com/dkfz-unite/unite-commands) web service, which runs a configured command for each request. The service is used by the [analysis service](https://github.com/dkfz-unite/unite-analysis) for [DEP analysis](https://github.com/dkfz-unite/unite-analysis/blob/main/Docs/api-analysis-dep.md): the analysis service writes the input files to the shared analysis data directory and calls the service with the analysis id as process key.

The script performs the following steps:
1. Reads the input data, metadata and options files.
2. Preprocesses the data (filters features, log2 transforms and normalizes data, and imputes missing values), see [Options](#options) below.
3. Performs differential expression analysis between the two conditions using the `lmFit`, `eBayes` and `topTable` functions of `limma`.
4. Saves the results to a file.

R packages: `limma`, `edgeR`, `sva`, `readr`, `dplyr`, `jsonlite` (installed by [src/install/install.R](src/install/install.R) on top of the base image `ghcr.io/dkfz-unite/docker-rdev-biocmanager`), `tibble`.

## Configuration
To configure the application, change environment variables as required in [commands](https://github.com/dkfz-unite/unite-commands/blob/main/README.md#configuration) web service. The [Dockerfile](Dockerfile) sets the following values:
- `UNITE_COMMAND` - command to run the analysis script (`Rscript`).
- `UNITE_COMMAND_ARGUMENTS` - command arguments (`run.R {data}/{proc}`).
  - `{data}` is replaced with `UNITE_DATA_PATH`, `{proc}` with the process key (analysis id).
- `UNITE_SOURCE_PATH` - location of the script in docker container, the working directory of the command (`/src`).
- `UNITE_DATA_PATH` - location of the data in docker container (`/mnt/data`).
- `UNITE_PROCESS_LIMIT` - maximum number of concurrent processes (`1` - process is heavy and uses a lot of CPU). If the limit is reached, the service returns `501`.

## Usage
- Place the files `data.tsv`, `metadata.tsv` and `options.json` to the `{proc}` subdirectory of the data directory (`/mnt/data/{proc}` in the container).
- Send a `POST` request to `/api/run?key={proc}`.
- The service runs `Rscript run.R /mnt/data/{proc}`, which reads the input files from this directory and writes `results.tsv` to it. The service returns `200` with the script output, or `500` if the script fails.

### Data format
Data file `data.tsv` should be in the following format:
```tsv
feature	sample1	sample2	sample3	sample4
protein1	10	20	30	40
protein2	15	25	35	45
protein3	20	30	40	50
protein4	25	35	45	55
```

Where:
- `feature` - identifier of the feature (protein). Should be first column. The analysis service uses Ensembl protein ids.
- `sample1`, `sample2`, `sample3`, `sample4` - names of the samples. The analysis service uses internal sample ids.
- Values in the table are raw protein intensity values (not normalised, not filtered). Missing values and `0` are treated as missing.

### Metadata format
Metadata file `metadata.tsv` should be in the following format:
```tsv
sample	condition	batch
sample1	A	1
sample2	A	2
sample3	B	1
sample4	B	2
```

Where:
- `sample` - name of the sample. Should be first column.
- `condition` - condition of the sample (e.g. control or treatment). The condition of the last row (in the example above `B`) becomes the reference category of the analysis. The analysis service uses dataset names as conditions, so the last dataset is the reference.
- `batch` - optional batch variable. If there is no batch variable, this column should be present but empty. If all samples have a batch and there are at least two batches, the batch is added as a fixed effect to the linear model.
- There should be exactly two conditions, with at least two samples for each condition.

The analysis service also writes `donor`, `specimen` and `specimen_type` columns; they are not used by the script.

### Options
`options.json` configures the preprocessing. All options are required.

- `normalization_log_offset`: [number > 0] log2 transformation applies an offset `log2(x + offset)`.
- `normalization_method`: [`"median"`, `"quantile"`] after log2 transformation, data are either median centered per sample (`"median"`) or quantile normalized with `normalizeQuantiles` of `limma` (`"quantile"`).
- `imputation_method`: [`"mindet"`, `"minprob"`] missing values are imputed using either:
  - `"mindet"` - each missing value is imputed with the "minimum" (actually 1st percentile) of intensity values of its sample.
  - `"minprob"` - imputed values are sampled from a Gaussian distribution centered on the "minimum" (see above). The standard deviation of this distribution is set to the median of the standard deviations of all proteins.
- `stratify_imputation_by_batch`: [`true`, `false`] if `true`, the imputation is done separately for each batch. If any protein has less than three non-missing values in any batch, the imputation falls back to non-stratified processing.
- `batch_correction_method`: [`"combat"`, `"limma"`, `null`] ignored by this analysis (the script sets it to `null`); see `batch` above.
- `min_non_missing_fraction`: [number, 0 < x <= 1] the minimum proportion of non-missing values of a protein required for it to be retained for analysis.
- `require_min_fraction_one_class`: [`true`, `false`] if `true`, a protein is retained if it reaches `min_non_missing_fraction` in at least one condition; if `false`, across all samples.

```json
{
  "normalization_log_offset": 0.1,
  "normalization_method": "median",
  "imputation_method": "mindet",
  "stratify_imputation_by_batch": false,
  "batch_correction_method": null,
  "min_non_missing_fraction": 0.7,
  "require_min_fraction_one_class": true
}
```

### Results format
`results.tsv` - `limma` `topTable` results for the condition coefficient, sorted by significance:
- `feature` - feature identifier.
- `logFC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`, `B` - `topTable` columns.
- `contrast` - compared conditions as `{test} - {reference}`.

```tsv
feature	logFC	AveExpr	t	P.Value	adj.P.Val	B	contrast
protein7	-2.0586	0.6670	-10.9903	4.69e-6	6.56e-5	4.7893	B - A
```

### Tests
Tests are in [src/tests](src/tests) with test data in [src/tests/test_data](src/tests/test_data).

## Installation

### Docker Compose
The easiest way to install the application is to use docker-compose:
- [Environment configuration and installation scripts](https://github.com/dkfz-unite/unite-env)
- [DEP analysis service configuration and installation scripts](https://github.com/dkfz-unite/unite-env/tree/main/applications/unite-analysis-dep)

### Docker
The image of the service is available in our [registry](https://github.com/dkfz-unite/unite-analysis-dep/pkgs/container/unite-analysis-dep).

The image is built by the [docker-images](.github/workflows/docker-images.yml) workflow on pushes to `main`: it checks out [unite-commands](https://github.com/dkfz-unite/unite-commands), publishes it with `publish-linux-x64.sh` and `publish-linux-arm64.sh`, copies the result (`commands` executable) to `./app`, builds amd64 and arm64 images with the [Dockerfile](Dockerfile) and pushes them under one `latest` manifest.

To build an image locally, publish unite-commands for the target platform to `./app` first:
```bash
# in the unite-commands repository
sh publish-linux-x64.sh # or publish-linux-arm64.sh

# in this repository
rm -rf app
cp -r ../unite-commands/publish/linux-x64 app
chmod a+rx app/commands
docker build -t unite.analysis.dep:latest .
```

All application components should run in the same docker network.
To create common docker network if not yet available run the following command:
```bash
docker network create unite
```

To run application in docker run the following command:
```bash
docker run \
--name unite.analysis.dep \
--restart unless-stopped \
--net unite \
--net-alias dep.analysis.unite.net \
-p 127.0.0.1:5310:80 \
-e ASPNETCORE_ENVIRONMENT=Release \
-v unite.analysis.data:/mnt/data:rw \
-d \
unite.analysis.dep:latest
```

The volume `unite.analysis.data` is the analysis data directory shared with the analysis service. For standalone use, a host directory can be mounted instead (e.g. `-v ./data:/mnt/data:rw`).
