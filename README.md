# PKD Real Data 程式說明

貝氏約束群組變數選擇於有序 Probit 迴歸模型，套用於 DAT-SPECT Parkinson's disease 資料集。

## 檔案需求

執行前需在工作目錄放置：

- `train_pkd.csv`
- `test_pkd.csv`

兩檔案欄位格式：

- 目標欄位 `Y`：1/2/3
- 主效應欄位：`Sex`, `Age`, 及各影像特徵
- 交互作用欄位：以底線分隔，如 `Sex_S.R`、`Age_AP.L`

## 相依套件

```r
install.packages(c("MASS", "mvtnorm", "truncnorm", "coda", "mcmcse"))
```

## 程式結構

| 函數 | 功能 |
|---|---|
| `compute_gamma()` | 依 η、φ（分組）、δ（階層約束）、ξ（反階層約束）遞迴計算 γ（各變數是否納入模型） |
| `calc_marginal_likelihood_ordinal()` | 計算給定變數子集合下 Y* 的邊際概似（用於 η 更新時的模型比較） |
| `run_single_chain()` | 單條 MCMC chain：依序更新 β（聯合抽樣）、η（逐一 collapsed Gibbs／邊際概似比較）、τ（MH 步）、Y*（截斷常態） |
| `select_median_model()` | 依後驗選入機率 > 0.5 決定 median probability model |
| `predict_ordinal_probit()` | 依 β、τ 估計值對新資料做類別預測 |
| `load_pkd_data()` | 讀取 CSV、依訓練集相關係數做主效應分組、整理交互作用欄位命名 |
| `build_constraint_matrices_pkd()` | 建構 φ（分組）、δ（交互作用對父變數的階層約束）、ξ（反階層，本設定為全 0）矩陣 |
| `cv_select_tau1sq()` | k-fold CV 選取最佳先驗變異數 `b_sq` |
| `find_burnin_realdata()` | 多條 chain 平行跑，以 τ 的 mpsrf（`coda::gelman.diag`）與 γ 的 MCSE（`mcmcse::mcse.mat`）判斷收斂，動態決定 burn-in 步數 |
| `run_pkd_analysis()` | 主流程：讀資料 → 建約束矩陣 → CV 選 `b_sq` → 找 burn-in → burn-in 抽樣 → 後驗抽樣 → 模型選擇 → 預測與混淆矩陣 |

## 主要參數（`run_pkd_analysis()`）

| 參數 | 說明 |
|---|---|
| `cor_threshold` | 主效應分組相關係數門檻 |
| `group_interactions` | 是否讓交互作用項繼承父影像變數的分組（`TRUE` = 將交互作用項分群；`FALSE` = 交互作用項各自獨立） |
| `b_sq` / `b_sq_grid` | 先驗變異數 b²；設 `NULL` 則透過 `b_sq_grid` 做 CV 選取 |
| `cv_folds`, `cv_burnin`, `cv_post` | CV 的 fold 數與每個 fold 的 burn-in／後驗步數 |
| `w` | η 的先驗選入機率 |
| `a_tau` | τ 之 MH 步的先驗標準差 |
| `n_chains`, `max_burnin`, `check_every` | 收斂偵測用的 chain 數、最大 burn-in、每次檢查間隔 |
| `grd_threshold`, `mcse_threshold` | mpsrf(τ) 與 max(MCSE(γ)) 之收斂門檻 |
| `n_post_burnin` | burn-in 後的正式後驗抽樣步數 |
| `preset_burnin` | 若指定正整數，則跳過動態收斂搜尋，直接使用固定 burn-in |
| `data_seed` | 控制 CV 切分、chain 初始化與正式抽樣的隨機種子 |

## 執行

執行後會依序印出：CV 各 `b_sq` 之平均 ACC → 選定的 `b_sq` → 收斂後的 burn-in 步數 → 選中變數（M 個）→ 各變數後驗選入機率表（僅列出 > 0.05 者）→ 訓練/測試集 ACC 與混淆矩陣。

3. 若需重新評估先驗變異數，將 `b_sq` 設回 `NULL` 並調整 `b_sq_grid`
本程式碼實作貝氏受限群組變數選擇之有序 Probit 迴歸模型，應用於 PKD（DAT-SPECT 巴金森氏症）資料集
