# circRNA features 與 predictCS workflow

這是可獨立搬移的 Snakemake release。每次分析都會對指定的 circRNA
完整計算 annotation、alternative splicing、conservation、MaxEntScan 與
predictCS，不讀取或寫入跨分析 run 共用的 feature cache。

固定 reference 所衍生的 annotation DB、AS site index、liftOver 結果、predictCS
reference DB 與 score rank indexes 會存放在 `results/references/`，並由 Snakemake
依輸入檔案時間戳記判斷是否需要重建。部分程式執行期間會使用暫存 SQLite 或
記憶體 LRU 控制記憶體與提高速度；這些不是可跨 run 重用的分析 cache。

## 輸入格式

`input_circrnas` 是四欄 tab-separated table：

```text
chrom	start	end	strand
```

- 座標使用 1-based closed convention。
- `strand` 必須是 `+` 或 `-`。
- 第一列可以沒有 header。
- 若有 header，chromosome 欄可使用 `chr`、`chrom`、`chromosome` 或
  `chrm`；座標欄可使用 `start`/`pos1` 與 `end`/`pos2`。
- chromosome 名稱會標準化為 UCSC style，例如 `1` 轉為 `chr1`、`MT`
  轉為 `chrM`。
- 不接受重複座標。

## Reference resources

執行前請在 `config.yaml` 設定以下檔案：

- GRCh38 genome FASTA 與 gene annotation GTF。
- hg38 phyloP 與 phastCons bigWig。
- human 與 mouse database-presence tables。
- mm10-to-hg38 chain file。
- circFL-seq hg19 table與 hg19-to-hg38 chain file。
- 包含四組 predictCS score 基準分布的 `scores_and_ranks.tsv`。
- 執行 checkAA 時另需一個 non-genomic reference FASTA。

相對路徑均以 release 資料夾為基準。大型 reference 檔不包含在這份 release
中。

## 執行主要分析

可直接修改 `config.yaml` 的 `input_circrnas` 與 `run_name`，也可由指令列覆蓋：

```bash
cd release
snakemake --use-conda --cores 8 \
  --config input_circrnas=/absolute/path/circRNAs.tsv run_name=example
```

主要輸出位於 `results/runs/<run_name>/`：

- `circRNAs.features.tsv`：predictCS 使用的 feature table。
- `circRNAs.features.detailed.tsv`：包含較完整 annotation 與 AS 欄位。
- `circRNAs.predictCS_scores.tsv`：四組 predictCS scores。
- `circRNAs.predictCS_score_ranks.tsv`：四組 predictCS percentile ranks。
- `circRNAs.results.tsv`：輸入座標、predictCS features、scores 與 ranks。
- `circRNAs.results.detailed.tsv`：詳細 features、scores 與 ranks。

`7 DB`、`FL-circAS or circFL_seq` 與 `mouse_conserved` 查無對應資料時為
`0`。若 conservation window 沒有完整 bigWig coverage，或 MaxEntScan window
超出 genome 範圍，相關 feature 保持空值。

四組 percentile rank 使用固定基準分布計算。若新 score 等於任一 reference
score，直接採用該 reference 同分群組的既有 average rank：

```text
((count(reference_score < score) + 0.5 * (count(reference_score == score) + 1)) / N) * 100
```

若 reference 中沒有相等的 score，則使用：

```text
(count(reference_score < score) / N) * 100
```

結果輸出至小數第 10 位。
predictCS score 為 `.` 時，對應 rank 亦為 `.`。基準分布會先建立可重用的排序
索引並存放於 `results/references/predictCS/score_rank_reference/`。

## 執行 checkAA

先設定 `config.yaml` 的 `checkAA.other_references`，再執行獨立 target：

```bash
snakemake --use-conda --cores 8 checkAA \
  --config input_circrnas=/absolute/path/circRNAs.tsv run_name=example
```

checkAA 使用 pblat，且不讀寫 checkAA cache。輸出包括：

- `circRNAs.checkAA.tsv`
- `circRNAs.results.with_checkAA.tsv`
- `circRNAs.results.detailed.with_checkAA.tsv`

## 重新執行

同一個 `run_name` 會指向相同輸出路徑，Snakemake 會依既有 output 與輸入時間戳記
決定是否重跑。要建立完全獨立的分析結果，請使用新的 `run_name`。
