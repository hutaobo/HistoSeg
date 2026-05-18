import fs from "node:fs";
import path from "node:path";
import {
  Presentation,
  PresentationFile,
  row,
  column,
  grid,
  panel,
  text,
  shape,
  rule,
  fill,
  hug,
  fixed,
  wrap,
  grow,
  fr,
  auto,
} from "@oai/artifact-tool";
import { paint } from "@oai/artifact-tool/presentation-jsx";

const ROOT = process.env.DECK_ROOT ?? "D:/GitHub/HistoSeg/artifacts/histoseg_anchor_only_residual_tps";
const OUT = path.join(ROOT, "output");
const SCRATCH = path.join(ROOT, "scratch");
const PREVIEWS = path.join(SCRATCH, "previews");

fs.mkdirSync(OUT, { recursive: true });
fs.mkdirSync(PREVIEWS, { recursive: true });

const W = 1920;
const H = 1080;

const C = {
  bg: "#F6F8FB",
  ink: "#12212F",
  muted: "#536579",
  faint: "#D8E1EA",
  soft: "#EAF0F6",
  blue: "#1E63B6",
  blueSoft: "#DCEBFF",
  teal: "#0F8C79",
  tealSoft: "#DDF5EF",
  amber: "#B46A00",
  amberSoft: "#FFF1D8",
  red: "#B64040",
  redSoft: "#FBE5E5",
  purple: "#6E4BB5",
  purpleSoft: "#EEE8FF",
  green: "#1A7F4B",
  greenSoft: "#E4F5EA",
  white: "#FFFFFF",
  black: "#0B1117",
};

const font = "Microsoft YaHei";
const mono = "Consolas";

const S = {
  title: { fontFace: font, fontSize: 52, bold: true, color: C.ink },
  subtitle: { fontFace: font, fontSize: 25, color: C.muted },
  h2: { fontFace: font, fontSize: 31, bold: true, color: C.ink },
  body: { fontFace: font, fontSize: 24, color: C.ink },
  small: { fontFace: font, fontSize: 17, color: C.muted },
  label: { fontFace: font, fontSize: 18, bold: true, color: C.muted },
  code: { fontFace: mono, fontSize: 20, color: C.black },
  codeSmall: { fontFace: mono, fontSize: 16, color: C.black },
  metric: { fontFace: mono, fontSize: 42, bold: true, color: C.ink },
  coverTitle: { fontFace: font, fontSize: 80, bold: true, color: C.ink },
  coverSubtitle: { fontFace: font, fontSize: 30, color: C.muted },
};

function T(value, style = S.body, opts = {}) {
  return text(value, { width: opts.width ?? fill, height: opts.height ?? hug, style, ...opts });
}

function Box(children, opts = {}) {
  return panel(
    {
      name: opts.name,
      width: opts.width ?? fill,
      height: opts.height ?? hug,
      padding: opts.padding ?? { x: 24, y: 20 },
      fill: paint(opts.fill ?? C.white),
      line: opts.line ?? { color: opts.stroke ?? C.faint, width: 1 },
      borderRadius: opts.radius ?? 12,
      columnSpan: opts.columnSpan,
      rowSpan: opts.rowSpan,
    },
    children,
  );
}

function Tag(label, color = C.blue, bg = C.blueSoft) {
  return panel(
    {
      name: `tag-${label}`,
      width: hug,
      height: hug,
      padding: { x: 14, y: 7 },
      fill: paint(bg),
      line: { color, width: 1 },
      borderRadius: 14,
    },
    T(label, { fontFace: font, fontSize: 16, bold: true, color }, { width: hug }),
  );
}

function Footer(source) {
  return row(
    { name: "footer", width: fill, height: hug, justify: "between", align: "center" },
    [
      T(source, S.small, { width: wrap(1320), name: "source-rail" }),
      T("HistoSeg 3D alignment logic", S.small, { width: hug }),
    ],
  );
}

function Slide(pres, title, subtitle, body, source = "Local repo state and A100 validation logs, 2026-05-08.") {
  const slide = pres.slides.add();
  slide.compose(
    panel(
      { name: "slide-bg", width: fill, height: fill, padding: { x: 76, y: 56 }, fill: paint(C.bg), line: { color: C.bg, width: 0 } },
      grid(
        {
          name: "slide-root",
          width: fill,
          height: fill,
          columns: [fr(1)],
          rows: [auto, fr(1), auto],
          rowGap: 28,
        },
        [
          column({ name: "title-stack", width: fill, height: hug, gap: 10 }, [
            T(title, S.title, { name: "slide-title", width: wrap(1600) }),
            subtitle ? T(subtitle, S.subtitle, { name: "slide-subtitle", width: wrap(1460) }) : T("", S.subtitle),
          ]),
          body,
          Footer(source),
        ],
      ),
    ),
    { frame: { left: 0, top: 0, width: W, height: H }, baseUnit: 8 },
  );
  return slide;
}

function Metric(label, value, detail, color = C.blue, soft = C.blueSoft) {
  return Box(
    column({ width: fill, height: hug, gap: 8 }, [
      T(label, { fontFace: font, fontSize: 17, bold: true, color }, { width: fill }),
      T(value, S.metric, { width: fill }),
      T(detail, S.small, { width: fill }),
    ]),
    { fill: soft, stroke: color, padding: { x: 22, y: 18 } },
  );
}

function Arrow(label = "->") {
  return T(label, { fontFace: mono, fontSize: 30, bold: true, color: C.muted }, { width: hug });
}

function ProcessNode(title, lines, fillColor, lineColor) {
  return Box(
    column({ width: fill, height: hug, gap: 10 }, [
      T(title, { fontFace: font, fontSize: 25, bold: true, color: lineColor }, { width: fill }),
      ...lines.map((line) => T(line, S.small, { width: fill })),
    ]),
    { fill: fillColor, stroke: lineColor, padding: { x: 24, y: 20 }, height: fixed(170) },
  );
}

function TinyPoint(color, label) {
  return row({ width: hug, height: hug, gap: 8, align: "center" }, [
    shape({ geometry: "ellipse", width: fixed(15), height: fixed(15), fill: paint(color), line: { color, width: 0 } }),
    T(label, S.small, { width: hug }),
  ]);
}

function Formula(lines) {
  return Box(
    column({ width: fill, height: hug, gap: 16 }, lines.map((line) => T(line, S.code, { width: fill }))),
    { fill: "#F8FBFF", stroke: C.blue, padding: { x: 30, y: 24 } },
  );
}

function LegendRow(items) {
  return row({ width: fill, height: hug, gap: 24, align: "center" }, items.map(([color, label]) => TinyPoint(color, label)));
}

function CodeCell(code, label) {
  return Box(
    column({ width: fill, height: hug, gap: 6 }, [
      T(code, S.codeSmall, { width: fill }),
      T(label, S.small, { width: fill }),
    ]),
    { fill: C.white, stroke: C.faint, padding: { x: 16, y: 12 } },
  );
}

function ComparisonRow(leftTitle, leftText, rightTitle, rightText, leftColor, rightColor) {
  return grid(
    { width: fill, height: hug, columns: [fr(1), fr(1)], columnGap: 28 },
    [
      Box(
        column({ width: fill, height: hug, gap: 12 }, [
          T(leftTitle, { fontFace: font, fontSize: 28, bold: true, color: leftColor }),
          T(leftText, S.body),
        ]),
        { fill: leftColor === C.red ? C.redSoft : C.amberSoft, stroke: leftColor, padding: { x: 28, y: 24 } },
      ),
      Box(
        column({ width: fill, height: hug, gap: 12 }, [
          T(rightTitle, { fontFace: font, fontSize: 28, bold: true, color: rightColor }),
          T(rightText, S.body),
        ]),
        { fill: rightColor === C.green ? C.greenSoft : C.tealSoft, stroke: rightColor, padding: { x: 28, y: 24 } },
      ),
    ],
  );
}

function addCover(pres) {
  const slide = pres.slides.add();
  slide.compose(
    panel(
      { width: fill, height: fill, fill: paint(C.bg), line: { color: C.bg, width: 0 }, padding: { x: 78, y: 62 } },
      grid(
        { width: fill, height: fill, columns: [fr(1.05), fr(0.95)], rows: [fr(1), auto], columnGap: 56 },
        [
          column({ width: fill, height: fill, justify: "center", gap: 24 }, [
            row({ width: fill, height: hug, gap: 12 }, [
              Tag("major upgrade", C.purple, C.purpleSoft),
              Tag("partial-anchor", C.teal, C.tealSoft),
              Tag("A100 validated", C.blue, C.blueSoft),
            ]),
            T("HistoSeg 3D 对齐算法逻辑", S.coverTitle, { width: wrap(860), name: "cover-title" }),
            T("Anchor-only Residual TPS: 从“边界强制重合”转为“证据锚点驱动”的软配准。", S.coverSubtitle, {
              width: wrap(780),
              name: "cover-subtitle",
            }),
          ]),
          Box(
            column({ width: fill, height: fill, justify: "center", gap: 28 }, [
              row({ width: fill, height: hug, justify: "center", gap: 22, align: "center" }, [
                shape({ geometry: "ellipse", width: fixed(42), height: fixed(42), fill: paint(C.blue), line: { color: C.blue, width: 0 } }),
                shape({ geometry: "ellipse", width: fixed(42), height: fixed(42), fill: paint(C.purple), line: { color: C.purple, width: 0 } }),
                shape({ geometry: "ellipse", width: fixed(42), height: fixed(42), fill: paint(C.teal), line: { color: C.teal, width: 0 } }),
              ]),
              Formula(["x_hard = H(x_raw)", "r_i = fixed_i - hard_aligned_i", "x_soft = x_hard + R(x_hard)"]),
              T("只有 active anchors 产生拉力；passive/no-counterpart contours 只跟随变换。", S.body, { width: wrap(650) }),
            ]),
            { height: fill, fill: C.white, stroke: C.faint, padding: { x: 46, y: 40 } },
          ),
          Footer("Prepared from HistoSeg main branch, GitHub issue #9, and A100 validation logs."),
        ],
      ),
    ),
    { frame: { left: 0, top: 0, width: W, height: H }, baseUnit: 8 },
  );
}

function authoredTable(headers, rowsData, widths, name) {
  const children = [];
  for (const h of headers) {
    children.push(Box(T(h, { fontFace: font, fontSize: 16, bold: true, color: C.white }, { width: fill }), {
      name: `${name}-header`,
      fill: C.ink,
      stroke: C.ink,
      padding: { x: 10, y: 9 },
      radius: 4,
    }));
  }
  for (const rowData of rowsData) {
    for (const cell of rowData) {
      children.push(
        Box(T(String(cell.value ?? cell), cell.style ?? S.small, { width: fill }), {
          fill: cell.fill ?? C.white,
          stroke: C.faint,
          padding: { x: 10, y: 8 },
          radius: 4,
        }),
      );
    }
  }
  return grid({ name, width: fill, height: hug, columns: widths, rowGap: 6, columnGap: 6 }, children);
}

const presentation = Presentation.create({ slideSize: { width: W, height: H } });

addCover(presentation);

Slide(
  presentation,
  "为什么需要这次升级",
  "旧系统已经能找到 partial-anchor，但 soft TPS 仍可能把物理缺失或撕裂强行缝合。",
  column({ width: fill, height: fill, gap: 34, justify: "center" }, [
    ComparisonRow(
      "旧 semantic soft TPS",
      "假设同名结构应该整体重合，用边界到边界的最近点拟合 TPS。局部撕裂或缺失会变成“需要修正”的几何误差。",
      "新 anchor-only residual TPS",
      "只信任 hard 阶段筛出的高置信度 anchor。没有点对点证据的 contour 保持 passive，不参与拉伸。",
      C.red,
      C.green,
    ),
    row({ width: fill, height: hug, gap: 18, align: "center" }, [
      ProcessNode("Label-free hard", ["识别局部可匹配结构", "允许切片错位、缺失、撕裂"], C.tealSoft, C.teal),
      Arrow(),
      ProcessNode("旧 soft", ["又尝试全局边界吸附", "语义一致性压过物理容忍"], C.redSoft, C.red),
      Arrow(),
      ProcessNode("升级目标", ["软变形只来自证据锚点", "无证据区域不产生拉力"], C.greenSoft, C.green),
    ]),
  ]),
);

Slide(
  presentation,
  "升级后的核心不变量",
  "HistoSeg 3D stack 的对齐目标从“最大化轮廓重合”改为“最大化有证据区域的残差一致性，同时保留无证据区域的物理状态”。",
  grid({ width: fill, height: fill, columns: [fr(0.95), fr(1.05)], columnGap: 42 }, [
    column({ width: fill, height: fill, gap: 22, justify: "center" }, [
      Metric("判断对象", "anchor pairs", "不是整张切片，也不是所有同名 contour", C.purple, C.purpleSoft),
      Metric("通过标准", "residual + topology", "median/p90 residual、invalid geometry、negative Jacobian", C.teal, C.tealSoft),
      Metric("IoU 的地位", "context only", "anchor-only 模式不再用全局 IoU gain 一票否决", C.amber, C.amberSoft),
    ]),
    Box(
      column({ width: fill, height: fill, justify: "center", gap: 22 }, [
        T("这解决了原来的矛盾：", S.h2),
        T("1. Hard alignment 认为只应使用可信局部锚点。", S.body),
        T("2. 旧 semantic TPS 又默认同名结构应该整体贴合。", S.body),
        T("3. 现在 auto 模式会让 label-free-group 直接进入 anchor-only soft，不再回到 semantic TPS。", S.body),
      ]),
      { fill: C.white, stroke: C.faint, height: fill, padding: { x: 42, y: 36 } },
    ),
  ]),
);

Slide(
  presentation,
  "控制流：soft_alignment_mode=auto",
  "策略显式化后，hard 阶段的胜出方法决定 soft 阶段的语义。",
  column({ width: fill, height: fill, gap: 30, justify: "center" }, [
    row({ width: fill, height: hug, gap: 18, align: "center" }, [
      ProcessNode("pairwise hard winner", ["label-free-group", "contour-tps", "coda-image"], C.blueSoft, C.blue),
      Arrow(),
      ProcessNode("resolve soft mode", ["auto -> anchor-only for label-free-group", "auto -> semantic for legacy backends"], C.purpleSoft, C.purple),
      Arrow(),
      ProcessNode("write provenance", ["active_soft_alignment_mode", "fallback reason", "runtime and QC columns"], C.greenSoft, C.green),
    ]),
    authoredTable(
      ["User mode", "label-free-group", "contour-tps / coda-image", "Failure behavior"],
      [
        ["auto", "anchor-only", "semantic", "label-free anchor-only unsafe -> hard-only fallback"],
        ["anchor-only", "anchor-only", "anchor-only", "unsafe -> hard-only fallback"],
        ["semantic", "semantic", "semantic", "legacy acceptance logic"],
        ["none / --no-soft", "hard-only", "hard-only", "no soft artifacts"],
      ],
      [fr(0.85), fr(1.05), fr(1.25), fr(1.45)],
      "mode-table",
    ),
  ]),
  "Code map: src/histoseg/threed/multislice.py::_resolve_soft_alignment_mode and pairwise metrics writer.",
);

Slide(
  presentation,
  "Hard 阶段输出什么",
  "Anchor-only soft 不重新发明匹配。它继承 hard 阶段已经筛好的 partial-anchor evidence。",
  grid({ width: fill, height: fill, columns: [fr(1.05), fr(0.95)], columnGap: 42 }, [
    column({ width: fill, height: fill, gap: 18, justify: "center" }, [
      CodeCell("partial_matches_csv", "hard 阶段匹配表"),
      CodeCell("used_for_transform == True", "只有这些 anchor 进入 residual TPS"),
      CodeCell("aligned_moving_centroid_x/y", "hard transform 后的 moving source 点"),
      CodeCell("fixed_centroid_x/y", "fixed target 点"),
      CodeCell("match status", "active / passive / no-counterpart 的溯源基础"),
    ]),
    Box(
      column({ width: fill, height: fill, justify: "center", gap: 24 }, [
        LegendRow([
          [C.blue, "fixed contour"],
          [C.red, "moving hard-aligned"],
          [C.purple, "used anchor"],
        ]),
        grid({ width: fill, height: fixed(360), columns: [fr(1), fr(1), fr(1)], columnGap: 18 }, [
          shape({ geometry: "ellipse", width: fill, height: fixed(220), fill: paint(C.blueSoft), line: { color: C.blue, width: 3 } }),
          shape({ geometry: "ellipse", width: fill, height: fixed(220), fill: paint(C.redSoft), line: { color: C.red, width: 3 } }),
          shape({ geometry: "ellipse", width: fill, height: fixed(220), fill: paint(C.purpleSoft), line: { color: C.purple, width: 3 } }),
        ]),
        T("关键点：不是所有轮廓都能拉 TPS。只有 hard 阶段确认的 active anchors 能产生残差约束。", S.body),
      ]),
      { height: fill, fill: C.white, stroke: C.faint, padding: { x: 34, y: 32 } },
    ),
  ]),
  "Code map: src/histoseg/threed/label_free_alignment.py::run_anchor_only_residual_tps.",
);

Slide(
  presentation,
  "Residual Field TPS 的数学定义",
  "先用 hard transform 解决大位姿，再对小残差拟合平滑场，避免 TPS 直接处理大旋转/大平移。",
  grid({ width: fill, height: fill, columns: [fr(1), fr(1)], columnGap: 42 }, [
    column({ width: fill, height: fill, gap: 24, justify: "center" }, [
      Formula([
        "s_i = H(moving_centroid_i)",
        "f_i = fixed_centroid_i",
        "r_i = f_i - s_i",
        "R = RBFInterpolator(s_i -> r_i)",
        "x_soft = x_hard + R(x_hard)",
      ]),
      T("实现上用 residual field，而不是直接拟合 moving -> fixed。这样 passive geometry 在 hard 坐标系里只接受小的、锚点约束的松弛。", S.body),
    ]),
    Box(
      column({ width: fill, height: fill, gap: 18, justify: "center" }, [
        T("默认参数", S.h2),
        authoredTable(
          ["Parameter", "Default", "Purpose"],
          [
            ["bbox padding", "0.10", "union bbox 外扩，覆盖边缘 passive geometry"],
            ["identity padding", "16", "远离锚点处残差固定为 0"],
            ["RBF smoothing", "1e-4", "轻度平滑，降低过拟合"],
            ["Jacobian grid", "50 x 50", "检测 fold-over / self-crossing"],
          ],
          [fr(1), fr(0.75), fr(1.8)],
          "param-table",
        ),
      ]),
      { height: fill, fill: C.white, stroke: C.faint, padding: { x: 30, y: 26 } },
    ),
  ]),
);

Slide(
  presentation,
  "Identity padding",
  "在 fixed + hard-aligned moving 的 union bbox 外扩区域加入 residual=0 的虚拟点，限制 TPS 外推漂移。",
  grid({ width: fill, height: fill, columns: [fr(1.05), fr(0.95)], columnGap: 42 }, [
    Box(
      column({ name: "identity-padding-diagram", width: fill, height: fill, justify: "center", gap: 18 }, [
        grid({ width: fill, height: fixed(340), columns: [fr(1), fr(1), fr(1)], rows: [fr(1), fr(1), fr(1)], columnGap: 18, rowGap: 18 }, [
          shape({ geometry: "ellipse", width: fixed(24), height: fixed(24), fill: paint(C.faint), line: { color: C.faint, width: 0 } }),
          shape({ geometry: "ellipse", width: fixed(24), height: fixed(24), fill: paint(C.faint), line: { color: C.faint, width: 0 } }),
          shape({ geometry: "ellipse", width: fixed(24), height: fixed(24), fill: paint(C.faint), line: { color: C.faint, width: 0 } }),
          shape({ geometry: "ellipse", width: fixed(24), height: fixed(24), fill: paint(C.faint), line: { color: C.faint, width: 0 } }),
          column({ width: fill, height: fill, justify: "center", align: "center", gap: 18 }, [
            shape({ geometry: "ellipse", width: fixed(70), height: fixed(70), fill: paint(C.purple), line: { color: C.purple, width: 0 } }),
            T("active anchor cluster", { fontFace: font, fontSize: 20, bold: true, color: C.purple }, { width: wrap(260) }),
          ]),
          shape({ geometry: "ellipse", width: fixed(24), height: fixed(24), fill: paint(C.faint), line: { color: C.faint, width: 0 } }),
          shape({ geometry: "ellipse", width: fixed(24), height: fixed(24), fill: paint(C.faint), line: { color: C.faint, width: 0 } }),
          shape({ geometry: "ellipse", width: fixed(24), height: fixed(24), fill: paint(C.faint), line: { color: C.faint, width: 0 } }),
          shape({ geometry: "ellipse", width: fixed(24), height: fixed(24), fill: paint(C.faint), line: { color: C.faint, width: 0 } }),
        ]),
        LegendRow([
          [C.purple, "real residual anchors"],
          [C.faint, "identity padding anchors, residual=0"],
        ]),
      ]),
      { height: fill, fill: C.white, stroke: C.faint, padding: { x: 34, y: 28 } },
    ),
    column({ width: fill, height: fill, justify: "center", gap: 22 }, [
      T("为什么需要 padding", S.h2),
      T("如果真实 anchor 只集中在组织一侧，RBF/TPS 在远离 anchor 的区域会外推。identity padding 把无证据区域的残差拉回 0，避免 passive contour 被非物理地甩动。", S.body),
      T("这不是新增匹配证据；它只是数值稳定约束。", { fontFace: font, fontSize: 26, bold: true, color: C.amber }),
    ]),
  ]),
);

Slide(
  presentation,
  "Active / Passive / No-counterpart 的处理方式",
  "新的语义不是“所有 contour 都对齐”，而是“有证据者拉动，无证据者随动”。",
  column({ width: fill, height: fill, gap: 28, justify: "center" }, [
    authoredTable(
      ["Contour state", "Produces TPS force?", "Moves after soft?", "Interpretation"],
      [
        [
          { value: "Active anchor", style: { fontFace: font, fontSize: 20, bold: true, color: C.purple }, fill: C.purpleSoft },
          "Yes",
          "Yes, residual corrected",
          "hard 阶段确认的高置信度对应关系",
        ],
        [
          { value: "Passive follower", style: { fontFace: font, fontSize: 20, bold: true, color: C.teal }, fill: C.tealSoft },
          "No",
          "Yes, follows residual field",
          "可见但没有足够证据作为控制点",
        ],
        [
          { value: "No-counterpart", style: { fontFace: font, fontSize: 20, bold: true, color: C.amber }, fill: C.amberSoft },
          "No",
          "Yes, constrained by padding",
          "可能是缺失、撕裂、漂移或真实不连续",
        ],
      ],
      [fr(1.05), fr(0.85), fr(1.05), fr(1.85)],
      "state-table",
    ),
    Box(
      T("这正是截图里 partial-anchor 示例的语义：绿色 moving after soft 可以局部改善，但没有 anchor 的区域不会被强迫去覆盖蓝色 fixed。", S.body),
      { fill: C.white, stroke: C.faint, padding: { x: 28, y: 20 } },
    ),
  ]),
);

Slide(
  presentation,
  "Anchor-only 的验收逻辑",
  "通过标准从 IoU gain 切换为 anchor residual + topology QC。IoU 仍记录，但不作为 anchor-only 的主判据。",
  grid({ width: fill, height: fill, columns: [fr(0.9), fr(1.1)], columnGap: 42 }, [
    column({ width: fill, height: fill, gap: 20, justify: "center" }, [
      Metric("最低证据", "anchor count", "anchor 太少则 fallback hard-only", C.purple, C.purpleSoft),
      Metric("拟合质量", "median / p90 residual", "看 active anchors 是否被安全修正", C.blue, C.blueSoft),
      Metric("拓扑安全", "negative Jacobian", "默认阈值 0.001，50x50 grid", C.teal, C.tealSoft),
      Metric("几何有效性", "invalid geometry count", "变形后 polygon 不应大量失效", C.amber, C.amberSoft),
    ]),
    Box(
      column({ width: fill, height: fill, justify: "center", gap: 22 }, [
        T("semantic 模式仍沿用原有 IoU/soft acceptance 逻辑。", S.body),
        rule({ width: fill, stroke: C.faint, weight: 2 }),
        T("anchor-only 模式下：", S.h2),
        T("低 IoU 不自动代表失败。它可能表示系统正确保留了物理缺失、撕裂缝隙或非对应组织。", S.body),
        T("真正要检查的是：锚点残差是否降低，Jacobian 是否安全，passive/no-counterpart 是否只随动不拉扯。", S.body),
      ]),
      { height: fill, fill: C.white, stroke: C.faint, padding: { x: 36, y: 34 } },
    ),
  ]),
);

Slide(
  presentation,
  "输出和可追溯性",
  "每一对切片都会记录 active mode、anchor-only QC、fallback reason 和 review artifact，方便回看和文章统计。",
  grid({ width: fill, height: fill, columns: [fr(1), fr(1)], columnGap: 40 }, [
    column({ width: fill, height: fill, gap: 16, justify: "center" }, [
      CodeCell("pairwise_alignment_metrics.csv", "新增 active_soft_alignment_mode、anchor count、residual、Jacobian、runtime"),
      CodeCell("anchor_only_soft_tps/anchor_only_tps_summary.json", "单 pair 的 anchor-only residual TPS 摘要"),
      CodeCell("anchor_only_soft_tps/anchor_only_tps_landmarks.csv", "used_for_tps、landmark_kind、source/target/residual 坐标"),
      CodeCell("anchor_only_soft_tps/anchor_only_soft_aligned_contours.geojson", "soft 后 contour geometry"),
    ]),
    Box(
      column({ width: fill, height: fill, gap: 16, justify: "center" }, [
        T("Review HTML 需要回答的问题", S.h2),
        T("1. 哪些 contour 是 active anchors？", S.body),
        T("2. residual vectors 是否小且方向合理？", S.body),
        T("3. passive/no-counterpart 是否没有被强迫闭合？", S.body),
        T("4. negative Jacobian ratio 是否接近 0？", S.body),
        T("5. fallback 是证据不足，还是拓扑不安全？", S.body),
      ]),
      { fill: C.white, stroke: C.faint, height: fill, padding: { x: 34, y: 30 } },
    ),
  ]),
);

Slide(
  presentation,
  "与 local-z orientation 的关系",
  "本升级改变的是 3D stack contour alignment 的 soft 阶段；transcript local-z correction 的坐标公式不变。",
  column({ width: fill, height: fill, gap: 30, justify: "center" }, [
    row({ width: fill, height: hug, gap: 18, align: "center" }, [
      ProcessNode("Contour stack alignment", ["更可靠的 XY/contour 位置", "partial-anchor 语义保留"], C.greenSoft, C.green),
      Arrow(),
      ProcessNode("Contour-aware local-z scoring", ["按 aligned contour 建上下层 transcript profile", "用相邻 contour continuity 判断 preserve/reverse"], C.blueSoft, C.blue),
      Arrow(),
      ProcessNode("Transcript-only z output", ["local_z_corrected_um", "z_3d_um"], C.purpleSoft, C.purple),
    ]),
    Box(
      T("细胞/轮廓云输出不因为 transcript local-z flip 而改变；local-z 修正只作用于 transcript 层。Anchor-only TPS 让 local-z 的 contour 单元更接近真实切片物理状态。", S.body),
      { fill: C.white, stroke: C.faint, padding: { x: 30, y: 22 } },
    ),
  ]),
);

Slide(
  presentation,
  "A100 验证：已证明的部分",
  "Pairwise smoke 和 full run 前几对显示 anchor-only residual TPS 能在真实 polyp pair 上安全接受。",
  column({ width: fill, height: fill, gap: 26, justify: "center" }, [
    row({ width: fill, height: hug, gap: 18 }, [
      Metric("Smoke pair", "accepted", "anchors=12, p90=0.0299, neg_jac=0.0", C.green, C.greenSoft),
      Metric("Wall time", "17.52s", "smoke pair runtime, max RSS 397 MB", C.blue, C.blueSoft),
      Metric("Topology", "0.0", "negative Jacobian ratio on smoke", C.teal, C.tealSoft),
    ]),
    authoredTable(
      ["Real pair", "Anchors", "Post p90 residual", "Negative Jacobian", "Soft IoU"],
      [
        ["001_to_002", "22", "28.7897", "0.000416", "0.725883"],
        ["002_to_003", "23", "39.3646", "0", "0.710260"],
        ["003_to_004", "21", "36.8260", "0.000416", "0.676023"],
        ["004_to_005", "16", "88.9183", "0.000416", "0.780452"],
      ],
      [fr(1), fr(0.65), fr(1.05), fr(1.05), fr(0.75)],
      "a100-table",
    ),
  ]),
  "A100 run root: /data/taobo.hu/projects/histoseg_polyp_xeniumslide/a100_anchor_only_residual_tps_20260508.",
);

Slide(
  presentation,
  "A100 full 32 run 的当前边界",
  "截至 2026-05-08，full 32 run 不是因为 anchor-only TPS 失败，而是在后续 contour generation 的 cluster normalization 上停止。",
  grid({ width: fill, height: fill, columns: [fr(1.05), fr(0.95)], columnGap: 42 }, [
    Box(
      column({ width: fill, height: fill, gap: 18, justify: "center" }, [
        T("失败点", S.h2),
        T("slice 6 contour generation", { fontFace: mono, fontSize: 28, bold: true, color: C.red }),
        T("ValueError: requested structures could not be matched after cluster normalization: Structure 1, Structure 2, Structure 5", S.codeSmall),
        rule({ width: fill, stroke: C.faint, weight: 2 }),
        T("wall time 28:24.75, max RSS 4,892,640 KB", S.small),
      ]),
      { fill: C.redSoft, stroke: C.red, height: fill, padding: { x: 34, y: 32 } },
    ),
    column({ width: fill, height: fill, gap: 22, justify: "center" }, [
      T("解释", S.h2),
      T("这说明新的 soft alignment 可以进入真实数据并写出前几对 accepted metrics，但完整 32 张验证还需要修复或绕过 contour materialization 的 structure normalization 不一致。", S.body),
      T("下一步应把 cluster/structure normalization 与 aligned contour reuse 逻辑修好，再重复 full run 并更新 issue #9 与 RTD。", S.body),
    ]),
  ]),
);

Slide(
  presentation,
  "代码层面的主要落点",
  "算法升级不是一个孤立函数，而是贯穿配置、调度、拟合、验收、输出和文档。",
  grid({ width: fill, height: fill, columns: [fr(1), fr(1)], columnGap: 40 }, [
    column({ width: fill, height: fill, gap: 16, justify: "center" }, [
      CodeCell("src/histoseg/threed/multislice.py", "CLI/config 接入；auto 解析；pairwise metrics；semantic/anchor-only 调度"),
      CodeCell("src/histoseg/threed/label_free_alignment.py", "partial-anchor matching；anchor-only residual TPS；identity padding；Jacobian QC"),
      CodeCell("tests/test_threed_multislice.py", "auto mode、CLI、fallback、synthetic tearing regression"),
      CodeCell("docs/3d_analysis.md", "方法说明、QC 指标、解释原则"),
    ]),
    Box(
      column({ width: fill, height: fill, justify: "center", gap: 18 }, [
        T("GitHub provenance", S.h2),
        T("Issue: hutaobo/HistoSeg#9", S.body),
        T("Commit: Add anchor-only residual TPS alignment", S.body),
        T("Commit: Document anchor-only residual TPS upgrade", S.body),
        T("A100 result comment records smoke success and full-run blocker.", S.body),
      ]),
      { fill: C.white, stroke: C.faint, height: fill, padding: { x: 34, y: 30 } },
    ),
  ]),
);

Slide(
  presentation,
  "结论：之前的矛盾是否解决了？",
  "控制流层面已经解决；完整 32 张生产验证还有一个 contour generation 输入一致性问题待处理。",
  column({ width: fill, height: fill, gap: 30, justify: "center" }, [
    ComparisonRow(
      "已解决",
      "auto + label-free-group 不再回退到 semantic soft TPS。Soft 变形只由 active anchors 拟合 residual field；passive/no-counterpart 不产生拉力。",
      "仍需收尾",
      "full 32 run 在 slice 6 的 structure normalization 阶段失败。它不否定 anchor-only TPS，但阻止了完整生产级统计。",
      C.green,
      C.amber,
    ),
    Box(
      column({ width: fill, height: hug, gap: 12 }, [
        T("当前可防守的表述", S.h2),
        T("HistoSeg 已完成从 semantic boundary attraction 到 evidence-first residual TPS 的算法转换；真实 polyp pair 的 smoke/partial validation 通过，full-stack validation 正在被一个下游 contour materialization issue 阻塞。", S.body),
      ]),
      { fill: C.white, stroke: C.faint, padding: { x: 30, y: 24 } },
    ),
  ]),
);

const pptx = await PresentationFile.exportPptx(presentation);
const pptxPath = path.join(OUT, "output.pptx");
await pptx.save(pptxPath);

const previewManifest = [];
for (let i = 0; i < presentation.slides.items.length; i += 1) {
  const slide = presentation.slides.items[i];
  const blob = await slide.export({ format: "png" });
  const buffer = Buffer.from(await blob.arrayBuffer());
  const previewPath = path.join(PREVIEWS, `slide_${String(i + 1).padStart(2, "0")}.png`);
  fs.writeFileSync(previewPath, buffer);
  previewManifest.push({ slide: i + 1, path: previewPath, bytes: buffer.length });
}

fs.writeFileSync(
  path.join(SCRATCH, "preview_manifest.json"),
  `${JSON.stringify({ pptx: pptxPath, slide_count: presentation.slides.items.length, previews: previewManifest }, null, 2)}\n`,
);

console.log(JSON.stringify({ pptx: pptxPath, slide_count: presentation.slides.items.length, preview_dir: PREVIEWS }, null, 2));
