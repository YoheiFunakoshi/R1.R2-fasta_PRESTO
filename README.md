# R1/R2 FASTQ merge app (PRESTO + cutadapt)

このリポジトリは、BCRペアリード（R1/R2 FASTQ）をマージするための
CLI/GUIアプリ一式です。cutadaptでのトリムと、pRESTOのAssemblePairsによる
マージを自動化します。Windowsでの再現性を重視しています。

## できること

- 3つのフローを選択可能
  - `merge-only`: マージのみ（トリムなし）
  - `trim-merge`: トリム→マージ
  - `merge-trim`: マージ→トリム
- 既定の固定配列（アダプター/プライマー）を使いつつ、サンプルごとに変更可能
- 出力は同一フォルダに集約（merged fastq / fasta / AP_align.log）

## 収録ファイル

- `bcr_merge_cli.py`: CLI本体
- `bcr_merge_gui.pyw`: GUI本体（ダブルクリック起動）
- `AIRR_outfmt19_PRESTO_AssemblePairs_R1R2_to_mergedFASTA_beta_1_ver1.0_plus_code.docx`
- `KKF340_R1R2_trimming_cutadapt_for_PRESTO_merge_ver1.0_plus_code.docx`

## 前提環境（推奨）

- Windows 10/11
- Miniforge または Anaconda/Miniconda
- Python 3.11
- `presto` と `cutadapt` が同一環境にインストール済み

環境の作り方（初回のみ）:

```bash
conda create -n presto_env -c conda-forge python=3.11 numpy scipy pandas biopython
conda activate presto_env
pip install presto cutadapt
python -c "import presto; print(presto.__version__)"
cutadapt --version
```

## GUIの使い方（おすすめ）

1. `bcr_merge_gui.pyw` をダブルクリックで起動
2. R1/R2を選択
3. Output base folderを確認
4. Modeを選択（初回は `merge-only` または `merge-trim` 推奨）
5. `Run` を押す

GUIが開かない場合（ショートカット作成）:

```
C:\miniforge3\envs\presto_env\pythonw.exe "C:\Users\Yohei Funakoshi\Desktop\R1.R2-fasta_PRESTO\bcr_merge_gui.pyw"
```

GUIの `Run inside conda env` をONにする場合:

- Conda exe: `C:\miniforge3\Scripts\conda.exe`
- Env name or path: `C:\miniforge3\envs\presto_env`

## CLIの使い方

### 1) merge-only（トリムなし）

```bash
python bcr_merge_cli.py --mode merge-only --r1 C:\path\R1.fastq --r2 C:\path\R2.fastq
```

### 2) trim-merge（トリム→マージ）

```bash
python bcr_merge_cli.py --mode trim-merge --r1 C:\path\R1.fastq --r2 C:\path\R2.fastq
```

### 3) merge-trim（マージ→トリム）

```bash
python bcr_merge_cli.py --mode merge-trim --r1 C:\path\R1.fastq --r2 C:\path\R2.fastq
```

conda環境を明示して実行する場合:

```bash
C:\miniforge3\Scripts\conda.exe run -p C:\miniforge3\envs\presto_env ^
  python bcr_merge_cli.py --mode merge-only --r1 C:\path\R1.fastq --r2 C:\path\R2.fastq
```

## 出力ファイルとフォルダ

出力は `Output base folder` の下に、モード名と日時付きのサブフォルダとして作成されます。

例:

```
C:\path\KKF340_withoutTrim_20260110_153000\
C:\path\KKF340_trimThenMerge_20260110_153500\
C:\path\KKF340_mergeThenTrim_20260110_154200\
```

モードにより以下のいずれかが作成されます（ファイル名にモード名と日時が入ります）。
すべて同一の実行フォルダ内に出力されます。

- `*_<mode>_<timestamp>_presto_assemble-pass.fastq`（AssemblePairsのマージ結果）
- `*_<mode>_<timestamp>_presto_assemble-pass.fasta`（FASTA変換）
- `*_<mode>_<timestamp>_AP_align.log`（AssemblePairsログ）
- `*_<mode>_<timestamp>_trim_R1.fastq`, `*_<mode>_<timestamp>_trim_R2.fastq`（trim-merge時）
- `*_<mode>_<timestamp>.fastq`, `*_<mode>_<timestamp>.fasta`（merge-trim時の最終出力）

## 固定配列（デフォルト）

必要に応じてGUIまたはCLIで上書きできます。

- R1 5' : `TGAGTTCCACGACACCGTCA`
- R2 5' : `CTAATACGACTCCGAATTCC`
- R1 3' : `GGAATTCGGAGTCGTATTAG`
- R2 3' : `TGACGGTGTCGTGGAACTCA`

## どのモードを使うべきか

- `trim-merge` は、トリム後にオーバーラップが短くなりPASSが減ることがあります。
- その場合は `merge-trim` が有効です（PASSを維持しやすい）。

実例（KKF340）:

- merge-only: PASS 75,248
- trim-merge: PASS 29,395
- merge-trim: PASS 75,248 → トリム後 75,020

## 注意点・トラブルシュート

- `cutadapt` の WARNING（adapter incomplete）は、プライマー/固定配列では
  よく出るため基本的に無視してOKです。
- Windowsで `.py` がVS Codeで開かれる問題があるため、必ず `python` 経由で実行します。
- `conda run -n presto_env` が失敗する場合は、環境のフルパス（`-p`）を使ってください。
- FASTAヘッダは空白以降を除去して出力します（IDの安定化のため）。
- R1/R2の向きが合わない場合は `--swap-r1r2` を試してください（GUIではチェックボックス）。
- cutadapt の警告「adapter is preceded by 'G' extremely often」は、
  プライマー/固定配列のトリムではよく出る注意喚起です。結果に大きな影響はないため、
  基本は無視してOKです。もしアダプター残存が疑われる場合のみ、3'アダプターの5'側に
  1塩基追加して再実行してください。
- cutadapt の「Indexing could take long」は速度/メモリの注意で、結果の正しさには影響しません。
- merge-trim では 5' と 3' の両方を確実に落とすため、single-end の cutadapt に `--times 2` を使用します。

## 再現性のために残すもの（推奨）

- 実行日時、作業フォルダ
- 入力FASTQ名（R1/R2）とサイズ
- 実行コマンド全文
- `AP_align.log`
- 出力FASTQ/FASTA名

## Gitについて

大きなFASTQや生成物は `.gitignore` により除外しています。
必要なら別途保存してください。
