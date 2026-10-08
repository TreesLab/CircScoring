# checkAA 內含原始碼

`checkAA_reads.py` 與 `checkaa_cache.py` 複製自
`/Users/chiangtw/Projects/Tools_dev/checkAA` 的 commit
`aa8c34c47be3be5ffb09a83640ff1b23ae805bdf`。

Snakemake 整合固定使用 `pblat` backend。`checkAA_reads.py` 仍保留舊有 BLAT
classes，讓內含的核心程式與上述版本完全一致；但 checkAA 環境不會安裝或呼叫
BLAT、samtools 或 `mp_blat.py`。
