---
name: analysis-report
description: >-
  Write step-by-step analysis reports as self-contained HTML (R Markdown), in the
  style "每一步：材料方法 → 结果 → 解读", with software versions, parameter tables,
  per-step statistics tables, embedded figures, and an interpretation paragraph for
  every result. Use for 分析报告、HTML 报告、流程报告、项目结题报告、给导师/合作者看的
  结果汇报, or whenever an analysis must be documented from raw data to final figures.
  Not for manuscript prose (use nature-writing / paper-writer) and not for slide decks.
license: CC-BY-4.0
---

# 分析报告写作：材料方法 → 结果 → 解读

把一次完整分析（从原始数据到最终图表）写成一份**可直接双击打开的 HTML 报告**。
核心不是排版，而是**每一步都交代清楚：用了什么软件/参数、得到什么数字、这些数字意味着什么**。

## 何时用

- 跑完一条分析流程，需要交付给导师、合作者或存档
- 项目结题 / 中期汇报 / 方法学附录
- 需要把散落在脚本和日志里的数字整理成可核查的文档

不用于：论文正文（用 `nature-writing`）、纯数据分析（不产出文档）。

## 报告骨架

```
摘要（一张总表：步骤 → 做了什么 → 关键结果）
步骤 1：xxx（软件 + 版本）
  材料方法  ← 软件、版本、参数、输入、阈值（要能复现）
  结果      ← 统计表 + 图，只放事实
  解读      ← 这个数字说明什么、能推出什么、有什么保留
步骤 2：…
…
讨论与局限
软件与版本（汇总表）
文件清单
```

**步骤划分按流程节点，不按图表**：质控 → 比对 → 定量 → 过滤 → 降维 → 差异 →
功能注释 → 建模 → 分层 → 富集。每一步都是"上一段的输出 = 下一段的输入"。

## 三条硬规则

1. **每个步骤必须有三段**（材料方法 / 结果 / 解读）。缺"材料方法"读者无法复现，
   缺"解读"就只是数据倾倒。
2. **"结果"里只放事实，"解读"里才下判断**。数字、表、图在结果；"说明/提示/意味着"
   在解读。两者混在一起，审阅者分不清哪句是数据、哪句是推断。
3. **软件 + 版本 + 参数必须写在"材料方法"里**，且与最后的版本汇总表一致。
   常见漏写：富集分析用的 R 包、建库/比对版本、阈值来源（软件推荐值还是人为选定）。

## 写作要点

详见 `references/style-guide.md`，这里只列最常踩的：

- **数字要能对上**：正文里的数必须等于表里的数；改了一个要全文搜一遍
- **阈值要交代来源**：软件推荐值 / 文献值 / 自己扫描后选定，三者写法不同
- **图必须被引用**：每张图前面要有"如图 N 所示"或在图注里编号，不能出现孤儿图
- **图注自足**：只读图注就能知道画的是什么、颜色/大小代表什么、n 是多少
- **保留不利信息**：过滤掉多少基因、哪些注释是假象、哪些步骤是近似——写清楚比藏起来安全
- **不写服务器路径与账号**：用 `data/`、`figs/` 这种相对路径
- **中英文混排**：软件名、参数、基因 ID 用英文；正文用中文；数字与单位之间留空格

## 技术实现

- **R Markdown → 自包含 HTML**（`self_contained: true`），图片 base64 内嵌，
  单文件即可传阅
- 表格用 `knitr::kable()`，不要手写 Markdown 表（列宽和对齐会崩）
- 图用 `ggsave()` 先出 PNG/PDF，Rmd 里用 `![](figs/x.png)` 引用，避免重复计算
- 渲染脚本与模板见 `templates/`
- **pandoc 必须在 PATH 里**（否则报 "pandoc version 1.12.3 or higher is required"）：
  `export PATH=/path/to/conda/envs/xxx/bin:$PATH`

## 交付

- 文件名带项目名与日期：`项目_分析内容_YYYY-MM-DD.html`
- 同时保留 `.Rmd` 源文件（改数字/改图后一条命令重出）
- 报告与图一起归档，报告内嵌图不依赖外部文件
- 归档前跑一遍 `references/style-guide.md` 末尾的自检清单

## 参考

- `references/style-guide.md`：逐段写作模板、措辞规则、自检清单、常见错误
- `templates/report_template.Rmd`：可直接改用的骨架
- `templates/render.R`：渲染脚本（含 pandoc PATH 处理）
