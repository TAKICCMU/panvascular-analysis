# 作图代码采用路径与最后修正

这些是实际历史代码，不是新的一键复现系统。本次只收集，不执行，也不根据Methods改写代码。源文件的原SHA256与原位置在06_Maps/ORIGINAL_PATH_TO_PACKAGE.tsv。

| 范围 | 原始构建与后续采用顺序 |
|---|---|
| Figure1 | current_builders/Figure1整图排版/build_figure1.py及figure1b_layout.py → update_clinical_percent_panels.py → assemble_current_approved_layout.py → final_patches/Figure1a_alignment/restore_native_alignment.py。最后一步修复人体与血管错位；旧build不能单独代表最终图。 |
| Figure2 | current_builders/Figure2整图排版/prepare_figure2.py为早期整体与a紧凑排版；最终g采用refine_g_single_row.py中的方案B，输入alternatives/source_data/option_B_final_Table10_rows.tsv。prepare脚本中的最早敏感性g不是最终采用g。2c实/虚线图例在最终导出中保留。 |
| Figure3 | visual_library/code/redraw_main.py等生成a–c；current_builders/Figure3整图排版/build_figure3.py生成d–f；updates/20261115_C_navy_gold为采用的后续颜色/标记修改；final_patches/Publication_delivery/prepare_figure_delivery.py补a符号图例等。 |
| Figure4 | current_builders/Figure4/code/build_merged_figure.py合并旧4/5保留b–d数值；final_patches/Figure4a/redraw_figure4a.py和install_plain_spatial_labels.py为最后简化workflow。 |
| ED1–10 | inherited_dependencies与visual_library/code提供原图与marker绘制；visual_library/code/refine_native_and_ed.py重排普通dot并保留原统计；final_patches/Publication_delivery/prepare_figure_delivery.py处理现行编号、内部图例与导出。原工作目录的部分panel来自已有SVG嵌入，并非全部从一个脚本一次生成。 |

当前最终SVG/PDF是04_Current_Figures中的版本。旧native_dependencies目录中的日期、旧Figure5或ED号只保留依赖追溯，不意味着另有当前版本。
最新配色/标签见palettes_and_labels及最终SVG。历史semantic_colors_current.json不必然覆盖后续BioRender局部改色。
未收集Python/R完整安装包、字体安装文件、BioRender在线工程内部状态或所有大型空间模型权重。已有requirements/environment文件随原代码保留；不声称已在新机器验证一键运行。
