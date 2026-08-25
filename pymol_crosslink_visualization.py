"""
PyMOL crosslink distance visualization script
================================================
Usage in PyMOL:
    run pymol_crosslink_visualization.py
    draw_crosslinks_for_pdb('PDB_ID')
    draw_all_structures()

Input CSV columns:
    PDB, chain_A, auth_seqid_A, chain_B, auth_seqid_B, min_distance
"""

from pymol import cmd
import os
import csv

# Paths can be supplied through environment variables. By default, the script
# reads structures and crosslinks.csv from the current working directory.
WORK_DIR = os.environ.get("XLMS_STRUCTURE_DIR", os.getcwd())
CSV_PATH = os.environ.get(
    "XLMS_CROSSLINK_CSV",
    os.path.join(WORK_DIR, "crosslinks.csv"),
)

# ── 颜色配置 ──────────────────────────────────────────────
PROTEIN_COLOR   = "skyblue"
CROSSLINK_COLOR = "orange"

# ── 线条配置 ──────────────────────────────────────────────
DASH_GAP    = 0.0
DASH_WIDTH  = 50.0
DASH_RADIUS = 0.80

# ── 小球配置 ──────────────────────────────────────────────
SPHERE_SCALE = 1
# ─────────────────────────────────────────────────────────


def _read_csv(csv_path):
    rows = []
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def _load_structure(pdb_code):
    object_name = pdb_code.lower()
    for ext in [".cif", ".pdb"]:
        fpath = os.path.join(WORK_DIR, pdb_code.lower() + ext)
        if os.path.exists(fpath):
            cmd.load(fpath, object_name)
            print(f"已加载: {fpath}")
            return object_name
    print(f"[错误] 找不到结构文件: {pdb_code.lower()}.cif 或 .pdb")
    return None


def draw_crosslinks_for_pdb(pdb_code, object_name=None, csv_path=CSV_PATH,
                              cutoff=40.0, show_labels=False):
    if object_name is None:
        object_name = pdb_code.lower()

    if object_name not in cmd.get_object_list():
        result = _load_structure(pdb_code)
        if result is None:
            return

    all_rows = _read_csv(csv_path)
    sub = [r for r in all_rows if r["PDB"] == pdb_code]

    if len(sub) == 0:
        print(f"[警告] CSV中没有找到 PDB={pdb_code} 的数据")
        return

    print(f"正在为 {pdb_code} 画 {len(sub)} 条crosslink...")

    sphere_sels = []
    n_drawn, n_skip, n_within, n_over = 0, 0, 0, 0

    for i, row in enumerate(sub):
        chain_a = row["chain_A"]
        resi_a  = int(float(row["auth_seqid_A"]))
        chain_b = row["chain_B"]
        resi_b  = int(float(row["auth_seqid_B"]))
        dist    = float(row["min_distance"])

        sel_a = f"{object_name} and chain {chain_a} and resi {resi_a} and name CA"
        sel_b = f"{object_name} and chain {chain_b} and resi {resi_b} and name CA"

        if cmd.count_atoms(sel_a) == 0 or cmd.count_atoms(sel_b) == 0:
            n_skip += 1
            continue

        dist_name = f"xl_{pdb_code}_{i}"
        cmd.distance(dist_name, sel_a, sel_b)
        cmd.color(CROSSLINK_COLOR, dist_name)
        if not show_labels:
            cmd.hide("labels", dist_name)

        sphere_sels.append(sel_a)
        sphere_sels.append(sel_b)

        n_drawn += 1
        if dist <= cutoff:
            n_within += 1
        else:
            n_over += 1

    print(f"完成：共画 {n_drawn} 条linker + {len(sphere_sels)} 个节点小球（跳过 {n_skip} 条）")
    print(f"  <= {cutoff}Å: {n_within} 条，> {cutoff}Å: {n_over} 条")

    # ── 1. 先设置蛋白外观 ──
    cmd.set_color("purplegrey", [0.75, 0.70, 0.85])
    cmd.set("dash_gap",    DASH_GAP)
    cmd.set("dash_width",  DASH_WIDTH)
    cmd.set("dash_radius", DASH_RADIUS)
    cmd.bg_color("white")
    cmd.show("cartoon", object_name)
    cmd.hide("lines",   object_name)
    cmd.color(PROTEIN_COLOR, object_name)          # 蛋白全部上紫灰色
    cmd.set("cartoon_transparency", 0.1, object_name)

    # ── 2. 再画小球（在蛋白颜色之后，防止被覆盖）──
    if sphere_sels:
        sphere_sel_str = " or ".join([f"({s})" for s in sphere_sels])
        sphere_name = f"xl_nodes_{pdb_code}"
        cmd.select(sphere_name, sphere_sel_str)
        cmd.show("spheres", sphere_name)
        cmd.set("sphere_scale", SPHERE_SCALE, sphere_name)
        cmd.color("orange", sphere_name)           # 强制橙色，覆盖蛋白颜色
        cmd.deselect()

    cmd.zoom(object_name)


def draw_all_structures(csv_path=CSV_PATH, cutoff=40.0):
    all_rows = _read_csv(csv_path)
    pdb_list = list(dict.fromkeys(r["PDB"] for r in all_rows))

    for pdb_code in pdb_list:
        cmd.reinitialize()
        draw_crosslinks_for_pdb(pdb_code, csv_path=csv_path, cutoff=cutoff)
        png_path = os.path.join(WORK_DIR, f"{pdb_code}_crosslinks.png")
        cmd.ray(2400, 2400)
        cmd.png(png_path, dpi=300)
        print(f"已保存: {png_path}\n")


print("Cross-link visualization functions loaded.")
print("  draw_crosslinks_for_pdb('PDB_ID')")
print("  draw_all_structures()")
