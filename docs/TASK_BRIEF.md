# Task brief — Intestinal oncofetal reference signature

> Founding task description (recorded 2026-10-01, verbatim). Not yet started.

---

我认为这个方向有潜力，而且比“直接沿用 Qin/revCSC signature”更适合作为你 TWEAKR paper 前半部分的方法学基础。但需要把定义收紧：你真正要建立的不是一个泛泛的 “oncofetal signature”，而是一个有明确证据层级的 intestinal oncofetal reference signature。

最关键的一点先回答：fetal 绝对不是天然 pan-fetal。 fetal brain、fetal liver、fetal intestine、placenta 有大量 tissue-specific developmental programs。你现在研究的 TACSTD2/TROP2–ANXA1–LY6A–CLU–YAP 这一套，本质上来自 fetal intestinal epithelium / injury-induced fetal-like regeneration，不能因为叫 fetal 就默认它代表所有胎儿组织。2013 年 Mustata 的原始工作直接证明 TACSTD2/Trop2 和 GJA1/Cx43 富集于 fetal intestinal spheroid，并在 E14 intestinal epithelium 表达；同时作者还特别提醒，Sca1/Ly6a 在 spheroid culture 中很高，但在 E15 intestinal epithelium 中很低，因此它甚至可能部分是 ex vivo/regenerative artifact，而不是纯粹的 in vivo fetal marker。

这反而给你的项目提供了一个很好的切入点。

## 我建议把整个 signature framework 分成三层，而不是两层

### Level 0 — Literature reference signatures

先不要创造任何东西。系统整理 field 已经使用过的 signatures，例如：

Mustata fetal spheroid；Yui fetal/regenerative；Nusse/Clu revival stem cell；Ayyaz revSC；Pikkupeura FEnS/fetal epithelium；Chen fetal/regenerative；Qin revCSC；以及后续 CRC revCSC signatures。

这是你的 benchmarking layer。

现在已有文献实际上经常把 fetal / regenerative / revival 混用。Chen 2023 就同时拿了 regeneration、fetal spheroid、revival stem cell 和 YAP signatures 做 GSEA，并发现 TGFB1 可以同时诱导它们。 (TGFB1 induces fetal reprogramming and enhances intestinal regeneration.pdf) 这说明这些 signatures 有共同 biology，但不能证明它们是同一个 biological object。

### Level 1 — Literature-supported Core Intestinal Oncofetal Markers

这才是你说的 core。

但我建议不要简单规定：“每个 gene ≥2 high-impact papers”。这会有 citation-count bias，而且容易 circular。应该要求每个 core gene 同时满足三个彼此独立的证据轴：

| Evidence axis | Core requirement |
|---|---|
| Fetal | fetal intestinal epithelium > adult intestinal epithelium |
| Onco | CRC malignant epithelium > normal adult epithelium |
| Literature replication | ≥2 independent studies / experimental systems support fetal/regenerative/oncofetal identity |
| Cross-species | mouse + human directionally concordant |
| Cell identity | epithelial/tumor-cell intrinsic，而不是 CAF/immune contamination |
| Robustness | replicated across datasets rather than driven by one cohort |

这样你得到的才是真正意义上的：Core Intestinal Oncofetal Program，而不是 “genes frequently called fetal-like”。

你现在提到的几个 gene，我不会全部预先判定为 Core。例如：

- **TACSTD2/TROP2** 是目前非常强的 candidate。Mustata 不只是 transcriptomic association，而是直接证明 Trop2 富集于 fetal spheroid，并在 embryonic intestinal epithelium 中表达；文章同时指出 Trop2 在多种 malignant tumors、包括 CRC 中升高。 (Identification of Lgr5-Independent Spheroid-Generating Progenitors of the Mouse Fetal Intestinal Epithelium.pdf)
- **CLU** 很强，但更接近 fetal-regenerative/revival axis。损伤后的 regenerative epithelium 明确出现 Clu-high population，而且与 homeostatic OLFM4 state 相反。 (TGFB1 induces fetal reprogramming and enhances intestinal regeneration.pdf)
- **ANXA1** 也很有希望。revSC literature 已经把 ANXA1、CLU、LY6A 作为典型 revSC markers； (Plastic persisters- revival stem cells in colorectal cancer.pdf) TGFB1-induced fetal reprogramming 中 Anxa1 chromatin accessibility 也会上升。 (TGFB1 induces fetal reprogramming and enhances intestinal regeneration.pdf)
- **LY6A/SCA1** 我反而建议暂时降级。它是极经典 regenerative/fetal-like marker，但 Mustata 明确指出 Sca1 在 spheroid 中高、E15 intestinal epithelium 中却 barely detected。 (Identification of Lgr5-Independent Spheroid-Generating Progenitors of the Mouse Fetal Intestinal Epithelium.pdf)

这是你这个项目特别值得做的地方：**field-known marker ≠ automatically bona fide fetal marker.** 你的 Atlas 可以真正把这件事重新清理一次。

### Level 2 — Data-derived Extended Intestinal Oncofetal Signature

这部分就完全可以发挥你的 Atlas 优势。

不要从 Core genes 扩相关基因。那会把 YAP、EMT、stress、cell cycle 等 downstream correlates 全吸进来。而应该重新从 phenotype 定义：

- mouse fetal intestine → mouse adult intestine
- human fetal intestine → human adult intestine
- human CRC malignant epithelial → adult normal epithelial

三个 contrasts independently 做。然后取方向一致的 intersection / meta-analysis candidates。

这样 Extended signature 的逻辑是：genes independently rediscovered from the definition of the biological state，而不是：genes correlated with genes that we already believed were oncofetal. 这是非常重要的方法学区别。

## 你的“相关性过滤”想法是对的，但我会把它放到最后

得到 Extended candidates 后，再利用 CRC Atlas 的几百万 epithelial cells 做内部结构。你会得到一个 gene–gene correlation/network。

这里真正有意思的不是简单删除低 correlation genes，而是看看是否出现 Core module：TACSTD2, ANXA1, CLU, EMP1, TNFRSF12A, CTGF/CCN2, CYR61/CCN1, ……

Pikkupeura 的 fetal epithelium 数据其实已经提供很漂亮的 mechanistic bridge：fetal progenitors 显示持续高 YAP activity，并高表达直接 YAP targets ANKRD1, CTGF, CYR61, CLU, TNFRSF12A；而且这种 YAP signature 在 E16.5 epithelium 相对 adult crypt 也富集。 (Transcriptional and epigenomic profiling identifies YAP signaling as a key regulator of intestinal epithelium maturation.pdf)

这对你的 TWEAKR paper 特别关键，因为 TNFRSF12A 本身就在 independently defined fetal YAP program 里面。这比“我们的 CRC 数据发现 TWEAKR 和 oncofetal score correlation”强很多。它建立了一条完全独立的 developmental evidence chain：

fetal intestine → YAP-high → TNFRSF12A-high

然后你的 cancer 数据再建立：

CRC → oncofetal/revCSC → TNFRSF12A-high

最后 perturbation 才问：

TWEAK → TWEAKR → YAP → state transition?

逻辑会非常干净。而且 YAP 在 CRC 中确实具有独立的 tumorigenic biology；已有 colon cancer 工作显示 YAP1 对部分 CRC cell proliferation/tumorigenesis 必需。 (YAP1 and PRDM14 converge to promote cell survival and tumorigenesis.pdf)

## 我甚至建议再增加一个非常重要的维度：tissue specificity

最终不要只有 Fetal vs Adult，而应该做：

- Fetal intestine vs Adult intestine
- Fetal intestine vs fetal brain/liver/lung/kidney/etc.

这样每个 gene 可以获得两个完全不同的 annotation：developmental fetalness 和 intestinal fetal specificity。于是你最终可能发现三类 genes：

1. **Pan-developmental fetal genes** — 多个 fetal tissues 都高。
2. **Intestinal-fetal genes** — fetal intestine 特别高。
3. **Regenerative fetal-mimic genes** — 真正 fetal intestine 未必非常高，但 injury/YAP/revSC 强烈诱导，例如某些 Ly6a-like genes。

这三类 biology 根本不一样。

甚至 placenta 最好单独处理，不要和 fetal organs 混起来。placenta/cancer 本身还有一整套独立的 oncofetal/onco-placental literature，例如 placenta-specific antigens 在正常成人组织低而在多种 cancer 中重新表达。 (biomedicines-11-00316.pdf)

## 最终我建议你发布的不是一个 signature，而是一个小型 signature family

| Signature | Size | Definition | Use |
|---|---|---|---|
| Core-OF | ~10–20 genes | 极严格、高可信、cross-species、fetal-high、CRC-high、literature replicated | IHC / Xenium / functional interpretation |
| Extended-OF | ~50–150 genes | 数据驱动 | UCell / AUCell / GSEA; scRNA / bulk / spatial scoring |
| Intestinal-Fetal | — | developmental reference，不要求 cancer-high | developmental reference |
| CRC-Oncofetal | — | fetal intestine + CRC intersection | — |

这比一个 50-gene list 同时承担所有功能合理得多。

## Figure 1 analytical architecture

- **Fig 1A** — Literature landscape：不同经典 fetal/regenerative/revSC/oncofetal signatures 的来源和 overlap。
- **Fig 1B** — Development axis：mouse + human fetal intestine vs adult intestine。
- **Fig 1C** — Cancer axis：CRC malignant epithelium vs normal adult epithelium。
- **Fig 1D** — intersection → Core / Extended OF。
- **Fig 1E** — cross-species replication。
- **Fig 1F** — single-cell gene–gene coherence/network，识别稳定 module 和 outliers。
- **Fig 1G** — benchmark：你的 Core/Extended OF vs Qin/revSC/fetal spheroid 等 signatures，在 independent fetal、regeneration、CRC datasets 上测试。

## Working definition

> **Oncofetal ≠ fetal marker + cancer marker.**
> **Oncofetal = a developmental epithelial program enriched in fetal intestine, extinguished during adult maturation, and reacquired by malignant epithelium.**

这样 TWEAKR story 就不再是“我们挑了 Qin signature，然后发现 TWEAKR 跟它相关”，而会变成：

> We reconstructed the intestinal oncofetal state from first principles → TNFRSF12A emerges as a conserved component → then ask what regulates it and what it does.

这在文章逻辑上明显更强，也更不容易被认为是 signature cherry-picking。
