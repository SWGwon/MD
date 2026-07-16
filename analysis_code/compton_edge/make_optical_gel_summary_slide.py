#!/usr/bin/env python3
"""Create a one-slide PNG/PDF/PPTX summary for the optical-gel comparison."""

from __future__ import annotations

import datetime as dt
import os
import textwrap
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

from PIL import Image, ImageDraw, ImageFilter, ImageFont


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "optical_gel_comparison"
OUT_STEM = RESULTS / "optical_gel_Co60_Cs137_summary_slide"

WIDTH, HEIGHT = 1920, 1080

COLORS = {
    "bg": "#F4F7FB",
    "ink": "#18212F",
    "muted": "#566273",
    "line": "#DCE3EC",
    "blue": "#1769E0",
    "blue_dark": "#0D3F91",
    "blue_light": "#EAF2FF",
    "red": "#D82C2C",
    "red_light": "#FDECEC",
    "green": "#17845A",
    "green_light": "#E8F6F0",
    "amber": "#9A6200",
    "amber_light": "#FFF5DE",
    "white": "#FFFFFF",
}

FONT_REGULAR = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
FONT_BOLD = "/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc"


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(FONT_BOLD if bold else FONT_REGULAR, size=size)


def rounded_card(base: Image.Image, box: tuple[int, int, int, int], radius: int = 22,
                 fill: str = COLORS["white"], outline: str = COLORS["line"],
                 shadow: bool = True) -> None:
    x0, y0, x1, y1 = box
    if shadow:
        layer = Image.new("RGBA", base.size, (0, 0, 0, 0))
        ld = ImageDraw.Draw(layer)
        ld.rounded_rectangle((x0 + 5, y0 + 8, x1 + 5, y1 + 8), radius=radius,
                             fill=(24, 33, 47, 30))
        layer = layer.filter(ImageFilter.GaussianBlur(10))
        base.alpha_composite(layer)
    draw = ImageDraw.Draw(base)
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=2)


def draw_text(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str, size: int,
              color: str = COLORS["ink"], bold: bool = False, anchor: str | None = None) -> None:
    draw.text(xy, text, font=font(size, bold), fill=color, anchor=anchor)


def fit_plot(source: Path, target_size: tuple[int, int]) -> Image.Image:
    image = Image.open(source).convert("RGB")
    # Remove a small outer white margin without cutting axes or legends.
    crop = (8, 8, image.width - 8, image.height - 8)
    image = image.crop(crop)
    target_w, target_h = target_size
    scale = min(target_w / image.width, target_h / image.height)
    resized = image.resize((int(image.width * scale), int(image.height * scale)),
                           Image.Resampling.LANCZOS)
    canvas = Image.new("RGB", target_size, "white")
    canvas.paste(resized, ((target_w - resized.width) // 2, (target_h - resized.height) // 2))
    return canvas


def draw_metric_bar(draw: ImageDraw.ImageDraw, x: int, y: int, width: int,
                    source: str, before: float, after: float, gain: float) -> None:
    draw_text(draw, (x, y), source, 29, bold=True)
    draw_text(draw, (x + width, y + 4), f"×{gain:.3f}", 30, COLORS["blue_dark"], True, "ra")
    y += 47
    bar_h = 24
    max_value = 2.25
    before_w = int(width * (1.0 / max_value))
    after_w = int(width * (gain / max_value))
    draw.rounded_rectangle((x, y, x + width, y + bar_h), radius=bar_h // 2,
                           fill="#EDF1F6")
    draw.rounded_rectangle((x, y, x + before_w, y + bar_h), radius=bar_h // 2,
                           fill=COLORS["blue"])
    draw.rounded_rectangle((x, y + 34, x + width, y + 34 + bar_h), radius=bar_h // 2,
                           fill="#EDF1F6")
    draw.rounded_rectangle((x, y + 34, x + after_w, y + 34 + bar_h), radius=bar_h // 2,
                           fill=COLORS["red"])
    draw_text(draw, (x + width, y + 2), f"gel 없음  {before:.1f} NPE", 20, COLORS["muted"], anchor="ra")
    draw_text(draw, (x + width, y + 36), f"gel 있음  {after:.1f} NPE", 20, COLORS["muted"], anchor="ra")


def create_slide() -> Image.Image:
    base = Image.new("RGBA", (WIDTH, HEIGHT), COLORS["bg"])
    draw = ImageDraw.Draw(base)

    # Header
    draw_text(draw, (70, 54), "Optical gel 적용 후 검출기 응답 ≈ 2.1×", 55, bold=True)
    draw_text(draw, (72, 124), "Co-60 · Cs-137  |  background 차감 CH0+CH1 합산  |  charge-derived nominal NPE",
              23, COLORS["muted"])

    # Takeaway banner
    banner = (70, 170, 1850, 270)
    draw.rounded_rectangle(banner, radius=24, fill=COLORS["blue_dark"])
    draw_text(draw, (104, 219), "두 선원 모두 gel 적용 후 약 2배 이동",
              36, COLORS["white"], True, "lm")
    draw.rounded_rectangle((1160, 187, 1485, 253), radius=18, fill="#2457A2")
    draw_text(draw, (1185, 220), "Co-60", 25, COLORS["white"], True, "lm")
    draw_text(draw, (1460, 220), "×2.138", 31, COLORS["white"], True, "rm")
    draw.rounded_rectangle((1505, 187, 1820, 253), radius=18, fill="#2457A2")
    draw_text(draw, (1530, 220), "Cs-137", 25, COLORS["white"], True, "lm")
    draw_text(draw, (1795, 220), "×2.117", 31, COLORS["white"], True, "rm")

    # Plot cards
    plot_y0, plot_y1 = 300, 740
    co_box = (70, plot_y0, 695, plot_y1)
    cs_box = (720, plot_y0, 1345, plot_y1)
    rounded_card(base, co_box)
    rounded_card(base, cs_box)
    draw = ImageDraw.Draw(base)
    draw_text(draw, (95, 326), "Co-60 total spectrum", 28, bold=True)
    draw_text(draw, (670, 326), "gel gain ×2.138", 24, COLORS["blue_dark"], True, "ra")
    draw_text(draw, (745, 326), "Cs-137 total spectrum", 28, bold=True)
    draw_text(draw, (1320, 326), "gel gain ×2.117", 24, COLORS["blue_dark"], True, "ra")

    co_plot = fit_plot(RESULTS / "Co60_gel_comparison_total_overlay.png", (575, 360))
    cs_plot = fit_plot(RESULTS / "Cs137_gel_comparison_total_overlay.png", (575, 360))
    base.paste(co_plot, (95, 360))
    base.paste(cs_plot, (745, 360))

    # KPI card
    kpi_box = (1370, plot_y0, 1850, plot_y1)
    rounded_card(base, kpi_box)
    draw = ImageDraw.Draw(base)
    draw_text(draw, (1400, 330), "수치 비교", 31, bold=True)
    draw_text(draw, (1815, 334), "selected central fit", 17, COLORS["muted"], anchor="ra")
    draw_metric_bar(draw, 1400, 382, 410, "Co-60", 234.2, 500.8, 2.138)
    draw.line((1400, 505, 1810, 505), fill=COLORS["line"], width=2)
    draw_metric_bar(draw, 1400, 530, 410, "Cs-137", 100.6, 213.0, 2.117)

    draw.rounded_rectangle((1400, 657, 1810, 716), radius=16, fill=COLORS["green_light"])
    draw_text(draw, (1420, 686), "P90 이동", 21, COLORS["green"], True, "lm")
    draw_text(draw, (1788, 686), "Co ×2.10  ·  Cs ×2.06", 22, COLORS["green"], True, "rm")

    # Bottom interpretation and caveats
    interp_box = (70, 775, 1160, 1003)
    caveat_box = (1190, 775, 1850, 1003)
    rounded_card(base, interp_box, shadow=False)
    rounded_card(base, caveat_box, fill=COLORS["amber_light"], outline="#F1D9A5", shadow=False)
    draw = ImageDraw.Draw(base)
    draw_text(draw, (98, 808), "해석", 29, bold=True)
    draw.rounded_rectangle((98, 855, 1130, 913), radius=16, fill=COLORS["green_light"])
    draw_text(draw, (118, 884), "결론: 두 선원에서 detector-level 전하 응답이 일관되게 약 2배 증가",
              25, COLORS["green"], True, "lm")
    draw_text(draw, (102, 946), "• 채널 비대칭도 재현: CH1 shift > CH0  →  coupling·geometry 영향이 함께 반영",
              22, COLORS["ink"])
    draw_text(draw, (1218, 808), "해석 시 주의", 29, COLORS["amber"], True)
    caveats = [
        "gel 자체의 절대 photon 투과율이 아님",
        "run별 SPE gain 및 재조립 systematic 미포함",
        "Cs-137 gel fit: χ²/ndf = 3.44, local solution 존재",
        "Co-60: 두 edge를 경험적 effective feature로 요약",
    ]
    cy = 856
    for item in caveats:
        draw_text(draw, (1220, cy), f"• {item}", 20, COLORS["ink"])
        cy += 35

    # Footer
    draw.line((70, 1031, 1850, 1031), fill=COLORS["line"], width=2)
    draw_text(draw, (72, 1052), "가정: 동일 HV·gain·trigger·source geometry  |  nominal gain = 1.0×10⁷",
              18, COLORS["muted"], anchor="lm")
    draw_text(draw, (1848, 1052), "2026-07-13", 18, COLORS["muted"], anchor="rm")
    return base.convert("RGB")


def write_minimal_pptx(image_path: Path, pptx_path: Path) -> None:
    """Package a single full-slide PNG into a minimal standards-compliant PPTX."""
    now = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
    files: dict[str, str | bytes] = {
        "[Content_Types].xml": """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Default Extension="png" ContentType="image/png"/>
  <Override PartName="/ppt/presentation.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.presentation.main+xml"/>
  <Override PartName="/ppt/slides/slide1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slide+xml"/>
  <Override PartName="/ppt/slideLayouts/slideLayout1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slideLayout+xml"/>
  <Override PartName="/ppt/slideMasters/slideMaster1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slideMaster+xml"/>
  <Override PartName="/ppt/theme/theme1.xml" ContentType="application/vnd.openxmlformats-officedocument.theme+xml"/>
  <Override PartName="/ppt/presProps.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.presProps+xml"/>
  <Override PartName="/ppt/viewProps.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.viewProps+xml"/>
  <Override PartName="/ppt/tableStyles.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.tableStyles+xml"/>
  <Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
  <Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
</Types>""",
        "_rels/.rels": """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="ppt/presentation.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>""",
        "docProps/core.xml": f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <dc:title>Optical gel Co-60 Cs-137 summary</dc:title><dc:creator>Codex</dc:creator>
  <dcterms:created xsi:type="dcterms:W3CDTF">{escape(now)}</dcterms:created>
  <dcterms:modified xsi:type="dcterms:W3CDTF">{escape(now)}</dcterms:modified>
</cp:coreProperties>""",
        "docProps/app.xml": """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes"><Application>Codex</Application><PresentationFormat>Widescreen</PresentationFormat><Slides>1</Slides><Notes>0</Notes><HiddenSlides>0</HiddenSlides></Properties>""",
        "ppt/presentation.xml": """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:presentation xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main">
  <p:sldMasterIdLst><p:sldMasterId id="2147483648" r:id="rId1"/></p:sldMasterIdLst>
  <p:sldIdLst><p:sldId id="256" r:id="rId2"/></p:sldIdLst>
  <p:sldSz cx="12192000" cy="6858000" type="screen16x9"/><p:notesSz cx="6858000" cy="9144000"/>
</p:presentation>""",
        "ppt/_rels/presentation.xml.rels": """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideMaster" Target="slideMasters/slideMaster1.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slide" Target="slides/slide1.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/presProps" Target="presProps.xml"/>
  <Relationship Id="rId4" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/viewProps" Target="viewProps.xml"/>
  <Relationship Id="rId5" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/tableStyles" Target="tableStyles.xml"/>
</Relationships>""",
        "ppt/slides/slide1.xml": """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:sld xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main"><p:cSld><p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/><a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr><p:pic><p:nvPicPr><p:cNvPr id="2" name="Optical gel summary"/><p:cNvPicPr><a:picLocks noChangeAspect="1"/></p:cNvPicPr><p:nvPr/></p:nvPicPr><p:blipFill><a:blip r:embed="rId2"/><a:stretch><a:fillRect/></a:stretch></p:blipFill><p:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="12192000" cy="6858000"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></p:spPr></p:pic></p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sld>""",
        "ppt/slides/_rels/slide1.xml.rels": """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout" Target="../slideLayouts/slideLayout1.xml"/><Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="../media/image1.png"/></Relationships>""",
        "ppt/slideLayouts/slideLayout1.xml": """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:sldLayout xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main" type="blank"><p:cSld name="Blank"><p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/><a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr></p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sldLayout>""",
        "ppt/slideLayouts/_rels/slideLayout1.xml.rels": """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideMaster" Target="../slideMasters/slideMaster1.xml"/></Relationships>""",
        "ppt/slideMasters/slideMaster1.xml": """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:sldMaster xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main"><p:cSld><p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/><a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr></p:spTree></p:cSld><p:clrMap accent1="accent1" accent2="accent2" accent3="accent3" accent4="accent4" accent5="accent5" accent6="accent6" bg1="lt1" bg2="lt2" folHlink="folHlink" hlink="hlink" tx1="dk1" tx2="dk2"/><p:sldLayoutIdLst><p:sldLayoutId id="1" r:id="rId1"/></p:sldLayoutIdLst><p:txStyles><p:titleStyle/><p:bodyStyle/><p:otherStyle/></p:txStyles></p:sldMaster>""",
        "ppt/slideMasters/_rels/slideMaster1.xml.rels": """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout" Target="../slideLayouts/slideLayout1.xml"/><Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/theme" Target="../theme/theme1.xml"/></Relationships>""",
        "ppt/theme/theme1.xml": """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<a:theme xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" name="Simple"><a:themeElements><a:clrScheme name="Simple"><a:dk1><a:srgbClr val="000000"/></a:dk1><a:lt1><a:srgbClr val="FFFFFF"/></a:lt1><a:dk2><a:srgbClr val="18212F"/></a:dk2><a:lt2><a:srgbClr val="F4F7FB"/></a:lt2><a:accent1><a:srgbClr val="1769E0"/></a:accent1><a:accent2><a:srgbClr val="D82C2C"/></a:accent2><a:accent3><a:srgbClr val="17845A"/></a:accent3><a:accent4><a:srgbClr val="9A6200"/></a:accent4><a:accent5><a:srgbClr val="566273"/></a:accent5><a:accent6><a:srgbClr val="DCE3EC"/></a:accent6><a:hlink><a:srgbClr val="0000FF"/></a:hlink><a:folHlink><a:srgbClr val="800080"/></a:folHlink></a:clrScheme><a:fontScheme name="Simple"><a:majorFont><a:latin typeface="Arial"/><a:ea typeface="Noto Sans CJK KR"/><a:cs typeface="Arial"/></a:majorFont><a:minorFont><a:latin typeface="Arial"/><a:ea typeface="Noto Sans CJK KR"/><a:cs typeface="Arial"/></a:minorFont></a:fontScheme><a:fmtScheme name="Simple"><a:fillStyleLst><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:fillStyleLst><a:lnStyleLst><a:ln w="9525"><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:ln></a:lnStyleLst><a:effectStyleLst><a:effectStyle><a:effectLst/></a:effectStyle></a:effectStyleLst><a:bgFillStyleLst><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:bgFillStyleLst></a:fmtScheme></a:themeElements></a:theme>""",
        "ppt/presProps.xml": """<?xml version="1.0" encoding="UTF-8" standalone="yes"?><p:presentationPr xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main"/>""",
        "ppt/viewProps.xml": """<?xml version="1.0" encoding="UTF-8" standalone="yes"?><p:viewPr xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main" lastView="sldView"><p:normalViewPr/><p:slideViewPr><p:cSldViewPr/></p:slideViewPr><p:notesTextViewPr><p:cViewPr varScale="1"/></p:notesTextViewPr><p:gridSpacing cx="78028800" cy="78028800"/></p:viewPr>""",
        "ppt/tableStyles.xml": """<?xml version="1.0" encoding="UTF-8" standalone="yes"?><a:tblStyleLst xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" def="{5C22544A-7EE6-4342-B048-85BDC9FD1C3A}"/>""",
        "ppt/media/image1.png": image_path.read_bytes(),
    }
    pptx_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(pptx_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, content in files.items():
            archive.writestr(name, content)


def main() -> None:
    for required in [
        RESULTS / "Co60_gel_comparison_total_overlay.png",
        RESULTS / "Cs137_gel_comparison_total_overlay.png",
    ]:
        if not required.exists():
            raise SystemExit(f"Required plot not found: {required}")

    slide = create_slide()
    png_path = OUT_STEM.with_suffix(".png")
    pdf_path = OUT_STEM.with_suffix(".pdf")
    pptx_path = OUT_STEM.with_suffix(".pptx")
    slide.save(png_path, format="PNG", optimize=True)
    slide.save(pdf_path, format="PDF", resolution=144.0)
    write_minimal_pptx(png_path, pptx_path)
    print(png_path)
    print(pdf_path)
    print(pptx_path)


if __name__ == "__main__":
    main()
