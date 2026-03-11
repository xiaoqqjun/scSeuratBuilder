# scSeuratBuilder

智能 Seurat 对象构建工具，支持强制 metadata 验证和基因名同步。

## 功能特点

- **三种数据结构支持**：自动检测并处理不同的数据组织方式
- **强制 metadata 验证**：确保样本信息与实验设计一致
- **多格式支持**：H5、MTX（三兄弟）、RDS、RDATA
- **基因名同步**：集成 geneSync，消除 Cell Ranger 版本差异
- **Seurat V4/V5 兼容**：自动检测对象类型并适配

## 安装

```r
# 从本地目录安装
devtools::install("path/to/scSeuratBuilder")
```

## 依赖包

- Seurat (>= 4.0.0)
- Matrix
- data.table
- readxl
- geneSync (可选，用于基因名同步)

## 快速开始

### 基本用法

```r
library(scSeuratBuilder)

# 创建 Seurat 对象
obj <- create_seurat_obj(
  data_dir = "path/to/your/data",
  species = "mmu",              # 小鼠，可选 "homo"(人), "rat"(大鼠)
  enable_gene_sync = TRUE       # 启用基因名同步
)
```

## 支持的数据结构

### Type A：子文件夹结构（推荐）

每个样本一个子文件夹，文件夹名与 metadata 中的样本名对应：

```
Data/
├── metadata.xlsx
├── Sample1/
│   └── filtered_feature_bc_matrix.h5
├── Sample2/
│   └── filtered_feature_bc_matrix.h5
└── Sample3/
    ├── barcodes.tsv
    ├── features.tsv
    └── matrix.mtx
```

### Type B：直接文件结构

数据文件直接放在 data_dir 中，文件名包含样本名：

```
Data/
├── metadata.xlsx
├── Sample1.h5
├── Sample2.h5
└── Sample3.h5
```

### Type C：已合并对象

单个已合并的 Seurat 对象文件：

```
Data/
├── metadata.xlsx
└── merged_seurat_object.rds
```

## metadata 文件格式

**必须**包含以下两列：

| samples | groups | ... |
|---------|--------|-----|
| Sample1 | Control | ... |
| Sample2 | Treatment | ... |
| Sample3 | Control | ... |

- `samples`：样本名称，必须与数据文件/文件夹名称匹配
- `groups`：样本分组信息

## 参数说明

```r
create_seurat_obj(
  data_dir,                     # 数据目录路径（必需）
  metadata_file = NULL,         # metadata 文件路径（默认在 data_dir 中查找）
  species = NULL,               # 物种："homo"/"human"/"hsa", "mus"/"mouse"/"mmu", "rat"/"rattus"/"rno"
  enable_gene_sync = FALSE,     # 是否启用基因名同步
  min_cells = 3,                # 最小细胞数阈值
  min_features = 100,           # 最小特征数阈值
  verbose = TRUE                # 是否显示详细信息
)
```

## 支持的文件格式

| 格式 | 扩展名 | 说明 |
|------|--------|------|
| Cell Ranger H5 | .h5 | Cell Ranger v3+ 输出 |
| MTX 三兄弟 | .mtx + .tsv | barcodes.tsv, features.tsv, matrix.mtx |
| RDS | .rds | R 序列化文件 |
| RDATA | .rdata, .rda | R 数据文件 |

## 输出

运行成功后，会在当前工作目录创建输出文件夹：

```
scSeuratBuilder_output_YYYY_MM_DD/
├── seurat_object.rds    # 合并后的 Seurat 对象
└── summary.txt          # 运行摘要
```

## 基因名同步（geneSync）

当 `enable_gene_sync = TRUE` 时，会自动调用 geneSync 包进行基因名标准化：

- 消除不同 Cell Ranger 版本之间的基因名差异
- 将基因名统一为权威符号（authority symbol）
- 过滤掉不可靠的基因名

**注意**：需要先安装 geneSync 包：

```r
devtools::install_local("path/to/geneSync")
```

## 示例

### H5 文件处理

```r
obj <- create_seurat_obj(
  data_dir = "G:/project/Data",
  species = "mmu",
  enable_gene_sync = TRUE,
  min_cells = 3,
  min_features = 100
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

## 工作流程

1. **验证 metadata 文件** - 检查必需列和格式
2. **检测数据结构** - 自动识别 Type A/B/C
3. **样本名匹配** - 验证数据与 metadata 一致性
4. **读取数据** - 创建 Seurat 对象
5. **基因名同步** - 可选，标准化基因名
6. **合并对象** - 整合 metadata 信息
7. **保存结果** - 输出到指定目录

## 注意事项

- metadata.xlsx 文件**必须**存在于数据目录中
- 样本名称必须与数据文件/文件夹名称**精确匹配**
- Type C（已合并对象）样本名不匹配时会显示警告但继续执行
- 大文件保存可能需要较长时间

## 作者

xiaoqqjun (xiaoqqjun@sina.com)

## 许可证

MIT License
