# scSeuratBuilder

> 智能单细胞 Seurat 对象构建工具
> 支持 metadata 验证、多格式数据读取、基因名同步、10X 文件验证
>
> `scSeuratBuilder` 提供两种构建模式：
> - **`create_seurat_obj()`** - 经典批量处理模式
> - **`create_seurat_obj_v2()`** - 增强版单样本循环模式

---

## 功能特点

### 核心功能

- **三种数据结构支持**：自动检测并处理不同的数据组织方式
- **强制 metadata 验证**：确保样本信息与实验设计一致
- **多格式支持**：H5、MTX（三兄弟）、RDS、RDATA
- **基因名同步**：集成 geneSync，消除 Cell Ranger 版本差异
- **Seurat V4/V5 兼容**：自动检测对象类型并适配

### V2 新增功能

- **10X 文件验证**：barcodes 格式一致性、features 列数检查、三文件对齐验证
- **断点续传**：支持从上次中断处继续处理
- **中间结果保存**：单样本处理完成后可保存 RDS 文件
- **细粒度错误处理**：单个样本失败不影响其他样本
- **节省存储空间**：默认跳过 scale.data 生成

---

## 安装

```r
# 从 GitHub 安装
remotes::install_github("xiaoqqjun/scSeuratBuilder")
# 或
devtools::install_github("xiaoqqjun/scSeuratBuilder")

# 安装依赖包 geneSync（可选，用于基因名同步）
remotes::install_github("xiaoqqjun/geneSync")
```

---

## 依赖包

| 包名 | 说明 | 必需 |
|------|------|------|
| Seurat | 单细胞分析核心包 | ✅ |
| Matrix | 稀疏矩阵处理 | ✅ |
| data.table | 高性能数据读取 | ✅ |
| readxl | Excel 文件读取 | ✅ |
| geneSync | 基因名同步 | 可选 |
| hdf5r | H5 文件支持 | 可选 |

---

# `create_seurat_obj()` - 经典批量处理模式

经典模式，适合快速构建小规模数据集。

## 基本用法

```r
library(scSeuratBuilder)

# 创建 Seurat 对象
obj <- create_seurat_obj(
  data_dir = "path/to/your/data",
  species = "mmu",              # 小鼠，可选 "homo"(人), "rat"(大鼠)
  enable_gene_sync = TRUE,       # 启用基因名同步
  output_dir = "output"          # 保存目录，NULL 则不保存
)
```

## 参数说明

```r
create_seurat_obj(
  data_dir,                     # 数据目录路径（必需）
  metadata_file = NULL,         # metadata 文件路径
  species = NULL,               # 物种："homo"/"human"/"hsa", "mus"/"mouse"/"mmu", "rat"/"rattus"/"rno"
  enable_gene_sync = FALSE,     # 是否启用基因名同步
  min_cells = 3,                # 最小细胞数阈值
  min_features = 100,           # 最小特征数阈值
  verbose = TRUE,               # 是否显示详细信息
  output_dir = NULL             # 输出目录，NULL 则不保存到文件
)
```

## 支持的数据结构

### Type A：子文件夹结构（推荐）

每个样本一个子文件夹：

```
Data/
├── metadata.xlsx
├── Sample1/
│   └── filtered_feature_bc_matrix.h5
├── Sample2/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
└── Sample3/
    └── filtered_feature_bc_matrix.h5
```

### Type B：直接文件结构

```
Data/
├── metadata.xlsx
├── Sample1.h5
├── Sample2.h5
└── Sample3.h5
```

### Type C：已合并对象

```
Data/
├── metadata.xlsx
└── merged_seurat_object.rds
```

## Metadata 文件格式

**必须**包含以下两列：

| samples | groups | (可选列) |
|---------|--------|----------|
| Sample1 | Control | condition |
| Sample2 | Treatment | batch |
| Sample3 | Control | ... |

- `samples`：样本名称，必须与数据文件/文件夹名称匹配
- `groups`：样本分组信息

## 示例

### H5 文件处理

```r
obj <- create_seurat_obj(
  data_dir = "G:/project/Data",
  species = "mmu",
  enable_gene_sync = TRUE,
  min_cells = 3,
  min_features = 100,
  output_dir = "G:/project/01.build_obj"
)
```

### 已合并对象处理

```r
obj <- create_seurat_obj(
  data_dir = "G:/project/MergedData",
  species = "homo",
  enable_gene_sync = TRUE
)
```

### 不保存到文件（仅返回对象）

```r
obj <- create_seurat_obj(
  data_dir = "path/to/Data",
  species = "human",
  enable_gene_sync = TRUE
  # output_dir = NULL (默认不保存)
)
```

## 工作流程

```
1. 验证 metadata 文件
2. 检测数据结构 (Type A/B/C)
3. 样本名匹配
4. 读取数据（批量）
5. 基因名同步（单样本级别）
6. 合并对象
7. 整合 metadata
8. 保存结果（可选）
```

---

# `create_seurat_obj_v2()` - 增强版单样本循环模式

增强版模式，适合大规模数据集、需要断点续传或更严格验证的场景。

## V2 新增特性

| 特性 | 说明 |
|------|------|
| **10X 文件验证** | 验证三文件对齐、barcodes 格式、features 列数 |
| **断点续传** | 处理中断后可从上次继续 |
| **中间保存** | 单样本结果可保存为 RDS |
| **错误恢复** | 单个样本失败不影响其他样本 |
| **进度可见** | 实时显示处理进度 |
| **节省存储** | 默认跳过 scale.data 生成 |

## 基本用法

```r
library(scSeuratBuilder)

# 基本使用
obj <- create_seurat_obj_v2(
  data_dir = "path/to/Data",
  metadata_file = "path/to/metadata.xlsx"
)

# 启用断点续传
obj <- create_seurat_obj_v2(
  data_dir = "path/to/Data",
  save_intermediate = TRUE,
  resume = TRUE
)

# 跳过问题样本
obj <- create_seurat_obj_v2(
  data_dir = "path/to/Data",
  skip_samples = c("bad_sample1", "bad_sample2"),
  on_sample_error = "skip"
)
```

## 参数说明

```r
create_seurat_obj_v2(
  # ===== 基础参数 =====
  data_dir,                     # 数据目录路径（必需）
  metadata_file = NULL,         # metadata 文件路径
  species = NULL,               # 物种
  enable_gene_sync = FALSE,     # 是否启用基因名同步
  min_cells = 3,                # 最小细胞数阈值
  min_features = 100,           # 最小特征数阈值
  verbose = TRUE,               # 是否显示详细信息
  output_dir = NULL,            # 输出目录

  # ===== V2 新增参数 =====
  save_intermediate = FALSE,    # 保存单样本 RDS 文件
  intermediate_dir = NULL,       # 中间文件目录
  resume = FALSE,               # 断点续传模式
  skip_samples = NULL,           # 要跳过的样本列表
  on_sample_error = "stop",      # 单样本错误处理: "stop" / "skip" / "warn"
  check_10x = TRUE,             # 验证 10X 文件
  skip_scale_data = TRUE        # 跳过 scale.data 生成
)
```

## V2 高级用法

### 断点续传

```r
# 第一次运行
obj <- create_seurat_obj_v2(
  data_dir = "path/to/Data",
  save_intermediate = TRUE,    # 必须启用中间保存
  intermediate_dir = "path/to/intermediate"
)

# 如果处理中断，可以继续运行
obj <- create_seurat_obj_v2(
  data_dir = "path/to/Data",
  save_intermediate = TRUE,
  resume = TRUE                 # 跳过已处理的样本
)
```

### 错误处理

```r
# 严格模式（默认）：任何样本失败都停止
obj <- create_seurat_obj_v2(
  data_dir = "path/to/Data",
  on_sample_error = "stop"
)

# 宽松模式：跳过失败的样本
obj <- create_seurat_obj_v2(
  data_dir = "path/to/Data",
  on_sample_error = "skip"
)

# 警告模式：继续但记录错误
obj <- create_seurat_obj_v2(
  data_dir = "path/to/Data",
  on_sample_error = "warn"
)
```

### 10X 文件验证

```r
# 启用验证（默认）
obj <- create_seurat_obj_v2(
  data_dir = "path/to/Data",
  check_10x = TRUE
)

# 跳过验证（数据已确认时）
obj <- create_seurat_obj_v2(
  data_dir = "path/to/Data",
  check_10x = FALSE
)
```

### 批量处理时保存中间结果

```r
obj <- create_seurat_obj_v2(
  data_dir = "path/to/Data",
  save_intermediate = TRUE,
  intermediate_dir = "path/to/intermediate",
  species = "human",
  enable_gene_sync = TRUE
)
```

## V2 工作流程

```
1. 验证 metadata 文件
2. 检测数据结构和格式
3. 验证 10X 文件（可选）
4. 样本名匹配
5. 创建中间目录（可选）
6. 单样本循环处理 ⭐
   ├── 读取数据
   ├── 创建 Seurat 对象
   ├── 基因同步（单样本级别）
   └── 保存中间文件（可选）
7. 合并对象
8. 整合 metadata
9. 返回最终对象
```

## V2 输出结构

```
data_dir/
├── metadata.xlsx
├── Sample1/
├── Sample2/
├── ...
└── .scSeuratBuilder_intermediate/    # 中间文件目录
    ├── Sample1.rds
    ├── Sample2.rds
    └── ...
```

---

# 10X 文件验证工具

包提供了独立的 10X 文件验证函数：

```r
library(scSeuratBuilder)

# 检查单个目录
result <- check_10x_triplets(
  data_dir = "path/to/10x/data",
  feature_col = 2,        # 使用第几列作为基因名（2=symbol, 1=ensembl_id）
  make_unique_names = TRUE,
  verbose = TRUE
)

if (result$valid) {
  mat <- result$matrix  # 获取验证通过的矩阵
} else {
  print(result$errors)
}

# 批量检查多个目录
dirs <- c("sample1/", "sample2/", "sample3/")
batch_result <- batch_check_10x(dirs, verbose = TRUE)

# 验证 Seurat 对象与原始文件的一致性
check_seurat_consistency(obj, "path/to/original/data/")
```

## 验证内容

| 检查项 | 说明 |
|--------|------|
| 文件存在性 | barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz |
| 维度对齐 | features行数 == matrix行数, barcodes行数 == matrix列数 |
| MTX索引 | 行/列索引不越界 |
| 名称质量 | NA值、空值、重复 |
| Barcodes格式 | 所有样本格式一致 |
| Features列数 | 检测是2列还是3列 |

---

# 基因名同步 (geneSync)

当 `enable_gene_sync = TRUE` 时，会自动调用 geneSync 包进行基因名标准化：

- 消除不同 Cell Ranger 版本之间的基因名差异
- 将基因名统一为权威符号（authority symbol）
- 过滤掉不可靠的基因名
- 支持人、小鼠、大鼠的跨物种转换

**注意**：需要先安装 geneSync 包：

```r
remotes::install_github("xiaoqqjun/geneSync")
```

---

# 支持的文件格式

| 格式 | 扩展名 | 说明 |
|------|--------|------|
| Cell Ranger H5 | .h5 | Cell Ranger v3+ 输出 |
| MTX 三兄弟 | .mtx/.tsv.gz | barcodes, features, matrix |
| RDS | .rds | R 序列化文件 |
| RDATA | .rdata, .rda | R 数据文件 |

---

# 两个函数的对比

| 特性 | create_seurat_obj | create_seurat_obj_v2 |
|------|-------------------|---------------------|
| **处理模式** | 批量处理 | 单样本循环 |
| **适用场景** | 小规模数据 | 大规模/长时运行 |
| **10X验证** | ❌ | ✅ |
| **断点续传** | ❌ | ✅ |
| **中间保存** | ❌ | ✅ |
| **进度显示** | 批量显示 | 逐样本详细显示 |
| **错误处理** | 全局停止 | 可配置跳过/继续 |
| **scale.data** | Seurat默认 | 默认跳过（节省空间） |
| **H5文件** | ✅ 支持 | ⚠️ 仅提示 |
| **基因同步** | ✅ 单样本级别 | ✅ 单样本级别 |
| **性能** | 快 | 相同（都是循环） |

---

# 输出

## create_seurat_obj 输出

当 `output_dir` 指定时：

```
output_dir/
├── seurat_object.rds    # 合并后的 Seurat 对象
└── summary.txt          # 运行摘要
```

## create_seurat_obj_v2 输出

当 `output_dir` 指定时：

```
output_dir/
├── seurat_object_v2.rds  # 合并后的 Seurat 对象
└── summary.txt           # 运行摘要
```

中间文件（当 `save_intermediate = TRUE` 时）：

```
data_dir/.scSeuratBuilder_intermediate/
├── Sample1.rds
├── Sample2.rds
└── ...
```

---

# 注意事项

- metadata 文件**必须**包含 `samples` 和 `groups` 两列
- 样本名称必须与数据文件/文件夹名称**精确匹配**
- H5 文件在 v2 中仅支持提示，建议使用 v1 处理
- 大文件保存可能需要较长时间
- 使用 `enable_gene_sync = TRUE` 需要安装 geneSync 包

---

# 常见问题

## Q: 两个函数应该如何选择？

**A:**
- 小规模数据（<10个样本）、快速分析 → `create_seurat_obj`
- 大规模数据（>10个样本）、长时间运行 → `create_seurat_obj_v2`
- 需要断点续传 → `create_seurat_obj_v2`
- 需要10X文件验证 → `create_seurat_obj_v2`

## Q: 如何启用基因名同步？

**A:**
```r
# 首先安装 geneSync
remotes::install_github("xiaoqqjun/geneSync")

# 然后设置参数
obj <- create_seurat_obj(
  data_dir = "path/to/Data",
  species = "human",        # 或 "mouse", "rat"
  enable_gene_sync = TRUE
)
```

## Q: H5 文件如何处理？

**A:**
- v1 (`create_seurat_obj`) 支持 H5 文件
- v2 (`create_seurat_obj_v2`) 检测到 H5 文件会提示不支持

## Q: 断点续传如何使用？

**A:**
```r
# 首次运行
obj <- create_seurat_obj_v2(
  data_dir = "path/to/Data",
  save_intermediate = TRUE
)

# 如果中断，再次运行相同代码即可继续
obj <- create_seurat_obj_v2(
  data_dir = "path/to/Data",
  save_intermediate = TRUE,
  resume = TRUE    # 关键参数
)
```

---

# 版本历史

- **v2.2.0** (2026-03-18)
  - 新增 `create_seurat_obj_v2()` 函数
  - 新增 10X 文件验证功能
  - 新增断点续传支持
  - 新增中间结果保存
  - 优化错误处理机制
  - 默认跳过 scale.data 生成

- **v1.0.0** (2025-02-27)
  - 初始版本
  - 基础 Seurat 对象构建功能
  - metadata 验证
  - 基因名同步支持

---

# 作者

**xiaoqqjun** (xiaoqqjun@sina.com)

- GitHub: [@xiaoqqjun](https://github.com/xiaoqqjun)
- ORCID: [0000-0003-1813-1669](https://orcid.org/0000-0003-1813-1669)
- WeChat: 博士后的小酒馆

---

# 许可证

MIT License

---

# 相关链接

- [GitHub 仓库](https://github.com/xiaoqqjun/scSeuratBuilder)
- [geneSync 包](https://github.com/xiaoqqjun/geneSync)
- [问题反馈](https://github.com/xiaoqqjun/scSeuratBuilder/issues)
