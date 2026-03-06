# Differential Protein Expression Analysis Service (Limma)

## General
R Limma based differential protein expression analysis application wrapped with web API.

The service performs the following steps:
1. Reads the input data and metadata files.
2. Filters the data to remove missing values.
3. Normalises the data using the `normalizeBetweenArrays` function from the `limma` package.
4. Performs differential expression analysis using the `lmFit` and `eBayes` functions from the `limma` package.
5. Saves the results to a file.

## Configuration
To configure the application, change environment variables as required in [commands](https://github.com/dkfz-unite/unite-commands/blob/main/README.md#configuration) web service:
- `UNITE_COMMAND` - command to run the analysis package (`Rscript`).
- `UNITE_COMMAND_ARGUMENTS` - command arguments (`run.R {data}/{proc}/data.tsv {data}/{proc}/metadata.tsv {data}/{proc}/results.tsv`).
- `UNITE_SOURCE_PATH` - location of the source code in docker container (`/src`).
- `UNITE_DATA_PATH` - location of the data in docker container (`/mnt/data`).
- `UNITE_PROCESS_LIMIT` - maximum number of concurrent jobs (`1` - process is heavy and uses a lot of CPU).

## Installation

### Docker Compose
The easiest way to install the application is to use docker-compose:
- Environment configuration and installation scripts: https://github.com/dkfz-unite/unite-environment
- Analysis service configuration and installation scripts: https://github.com/dkfz-unite/unite-environment/tree/main/applications/unite-analysis-dpe

### Docker
[Dockerfile](Dockerfile) is used to build an image of the application.
To build an image run the following command:
```
docker build -t unite.analysis.dpe:latest .
```

All application components should run in the same docker network.
To create common docker network if not yet available run the following command:
```bash
docker network create unite
```

To run application in docker run the following command:
```bash
docker run \
--name unite.analysis.dpe \
--restart unless-stopped \
--net unite \
--net-alias dpe.analysis.unite.net \
-p 127.0.0.1:5310:80 \
-e ASPNETCORE_ENVIRONMENT=Release \
-e UNITE_COMMAND=Rscript \
-v ./data:/mnt/data:rw \
-d \
unite.analysis.dpe:latest
```

## Usage
- Place the data files `data.tsv` and `metadata.tsv` in the `./data/{proc}` directory on the host machine.
- Send a POST request to the `localhost:5310/api/run?key=[key]` endpoint, where `[key]` is the process key.
- Analysis will run the command `Rscript` with the arguments `run.R {data}/{proc}/data.tsv {data}/{proc}/metadata.tsv {data}/{proc}/results.tsv` where `{proc}` is the process key.
  - All entries of `{data}` will be replaced with the path to the data location in docker container (In the example `./data` on the host machine will be mounted to `/mnt/data` in container).
  - All entries of `{proc}` will be replaced with the process key.
- Analysis will try to find the files `data.tsv` and `metadata.tsv` in the `{proc}` subfolder of the data location and use them as input.
- Analysis will save the results to the file `results.tsv` in the same subfolder.

### Data format
Data file `{proc}/data.tsv` should be in the following format:
```tsv
feature sample1 sample2 sample3 sample4
protein1 10 20 30 40
protein2 15 25 35 45
protein3 20 30 40 50
protein4 25 35 45 55
```  

Where:
- `feature` - identifier of the feature (protein). Should be first column.
- `sample1`, `sample2`, `sample3`, `sample4` - names of the samples.
- Values in the table are raw protein intensity values (not normalissed, not filtered).

### Metadata format
Metadata file `{proc}/metadata.tsv` should be in the following format:
```tsv
sample condition
sample1 A
sample2 A
sample3 B
sample4 B
```

Where:
- `sample` - name of the sample. Should be first column.
- `condition` - condition of the sample (e.g. control or treatment). Should be second column.
- Shoule be at least two samples for each condition.
- Should be exactly two conditions.
