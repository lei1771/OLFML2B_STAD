# OLFML2B-STAD GitHub 上传说明

这是用于论文代码公开的干净仓库版本。

## 上传前

1. **先解压本 ZIP。不要把 ZIP 本身作为 GitHub 仓库唯一文件上传。**
2. 在 GitHub 新建公开仓库，推荐名称：`OLFML2B_STAD`。
3. 将解压后 `OLFML2B_STAD/` 目录**里面的全部文件和文件夹**上传到仓库根目录。
4. 不要额外上传你本地工程中的 `data/`、`output/`、`logs/`、RDS/RData、HDF5、PDC/TIGER 原始文件。
5. 上传后创建 Release：tag `v1.0.0`。
6. 连接 Zenodo 后归档该 Release，获得 DOI。

## 推荐首个 commit message

`Initial reproducible release for OLFML2B-STAD manuscript`

## 运行

如果缺包，先运行：

```r
source("00_INSTALL_REQUIRED_PACKAGES_ONCE.R", encoding="UTF-8", local=FALSE)
```

完整运行：

```r
source("00_RECOVER_AND_RUN_OLFML2B_FROM_ZERO_v1_0_2.R", encoding="UTF-8", local=FALSE)
```

## 手工数据

- `DATA_INBOX/PDC000614/`：放 PDC000614 必需的原始下载文件。
- `DATA_INBOX/TIGER/`：放冻结设计所需的两个 TIGER exact RDS。

实际输入已被 `.gitignore` 屏蔽，避免误上传。
