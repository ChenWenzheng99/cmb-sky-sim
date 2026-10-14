# GitHub installation and validation scope

从仓库根目录执行 `python3 scripts/fetch_assets.py --cmb --theory-sources`，
再进入本目录运行 `./setup_cmb_env.sh`。该操作获取校验过的独立运行环境、
7 个依赖 wheel 和理论来源记录；大尺度结构及 CMB 天图仍由集群仓库提供。
原数值源代码和预计算输入谱保持原样。

下文的十套验收属于生成阶段的历史记录。文章采用的 50 套原生 NSIDE2048
统计摘要另见 [`../validation_50`](../validation_50/)。

# Sehgal2010 主场 + mmDL κ 的原生 CMB 生成扩展

本包从 Sehgal2010 参数的 CAMB 标量主场抽样，使用仓库中每个 ID 对应的
原生 κ 图，以 lenspyx 生成 NSIDE2048 的 T/Q/U。输出使用
`lensedcmb_new_XXXXX` 文件名前缀；原有 `lensedcmb_XXXXX.fits` 和其他源图保留。
生产脚本、理论、验证代码及离线依赖均在本包内。实际 κ FITS 为外部输入，
不包含在代码发布包内，也不读取 `A_NEW_WORK_MAIN` 的工具代码。

## 离线安装

完整包包含 `runtime_base.tar.gz`、`wheels/` 和
`evidence/offline_wheels_manifest.json`，面向 Linux x86_64/Python3.12。
先把整个包解压到最终使用的路径，再安装：

```bash
cd /your/final/path/cmb_sehgal2010_release
./setup_cmb_env.sh
./run_cmb.sh --help
```

安装器检查基础归档 SHA256 和七个 wheel 的大小、SHA256，解压到本包的
`runtime/`，在该最终前缀执行一次 `conda-unpack`，以
`pip --no-index --no-deps` 安装固定版本，并执行 `pip check` 和导入检查。
原 Ali_sky 运行环境不会被修改。已解包的 `runtime/` 不可搬迁；要换路径，
在目标路径重新解压完整包并安装。重跑安装器会复核并重装同一组离线 wheel。
也可传入基础归档的绝对路径：

```bash
./setup_cmb_env.sh /path/to/runtime_base.tar.gz
```

基础归档 SHA256：
`fa263b88f8f101e71d0f7de31cc68593202fa5d52ee7d5dc55b279d7dd94c577`。
七个覆盖包为 NumPy1.26.4、CAMB1.6.1、ducc0.38.0、contourpy1.3.3、
SymPy1.14.0、mpmath1.3.0 和 lenspyx2.0.0，详见 `requirements-extra.txt`。
基础环境还提供 healpy1.19.0、SciPy1.14.1、Astropy 和所需共享库。
本包的 lenspyx wheel 来自原始2.0.0数值/数据源文件；打包后端有调整，
原始43个库文件的核对记录见 `evidence/lenspyx_vendor_manifest.json`。

## 理论与模型假设

[Sehgal et al.2010, Section1](https://arxiv.org/abs/0908.0540) 给出
Ω_b=.044、Ω_m=.264、h=.71、n_s=.96、σ₈=.8 和四舍五入的 Ω_Λ=.736。
本包使用 H0=71、ω_b=.0221804、ω_c=.110902，并由线性 σ₈(z=0)=.8
确定 A_s。完成的理论计算得到 A_s=2.1500364844398247×10⁻⁹、
σ₈=.8000000000000018。CAMB flatness 包括辐射，实际 Ω_de=.73591700649。

未恢复原作者的完整 CAMB ini。以下是本计算明示的补全设置：τ=.08、
无质量中微子、N_eff=3.046、YHe=.24、TCMB=2.7255K、标量 pivot=.05/Mpc、
r=0、零 running、Ω_K=0、w=−1、wa=0，以及 CAMB 默认 TanhReionization/
Recfast 和 original Smith2003 Halofit。论文气体模型的近似光深为.08，
仅作为 τ 选择的动机，不等于找到原始主场的确切 τ。
作者 SED 代码中的2.726K也不能认证 CAMB 背景温度。
来源快照、未恢复设置和数值配置见 `theory/source_parameter_audit.md`、
`theory/sehgal2010_config.json` 和 `theory/source_evidence/source_manifest.json`。

`theory/outputs/sehgal2010_theory.npz` 的主谱为 raw Cℓ，单位 μK_CMB²，
零起点索引0..6500，列顺序 **TT、EE、BB、TE**。r=0 的原初 BB 严格为零。
主场截止ℓ=m=5999，κ/梯度偏折截止L=m=6143，κ `map2alm` 使用 iter3，
无权重和像素权重；κ→φ 为2/[L(L+1)]，φ₀=0而L=1保留，
梯度偏折为√[L(L+1)]φ。lenspyx 使用原生 HEALPix2048、epsilon=10⁻⁶。
输出为 float64、RING、COSMO 约定的 T/Q/U，单位 μK_CMB。
生成阶段不加 beam/mask、不作 pixel-window 修正或降采样。

CAMB 的 φ 谱是理论参考，实际重映射输入是既有 mmDL κ，二者不被混用。
本包还保存主场5999/φ6143的截止匹配理论、全角积分采样控制和零透镜控制。
完整 CAMB 参考、截止匹配 CAMB 参考、使用实测平均 κ 功率的条件参考各自
标明；后者不包含非高斯 κ 的全部高阶统计。详细数组约定见 `theory/README.md`。

## 坐标、种子和文件

生产脚本保留 κ 的方向，不执行坐标旋转。有有效 FITS `COORDSYS` 时，
该头信息优先决定输出坐标；头信息缺失时默认采用仓库所有者声明的 G。
这项声明不同于上游 mmDL ICRS 来源说明，本生成器不据此认证曾发生过未记录的
旋转。`--owner-coordinate C` 或 `E` 只改变缺失头信息时的声明，不旋转图。
各文件的 JSON 会保留该来源判断和限制。

默认主种子为20261011。`SeedSequence` 根据固定域、主种子和 ID 产生各 ID
的 uint32 种子，调度顺序和线程数不改变种子。这不能恢复旧生成器未记录的
随机抽样，也不保证 uint32 永不碰撞。要更换抽样，显式设置 `--master-seed`。

对 ID0，输入为 `DATA_ROOT/00000/kappa_00000.fits`，输出为：

- `00000/lensedcmb_new_00000.fits`
- `00000/lensedcmb_new_00000.json`
- `00000/lensedcmb_new_00000.input_spectra.npz`

三个文件记录 κ/理论/代码的哈希、种子、运行参数和全图有限值检查。
相同完整契约可校验后复用；部分输出或契约不符会拒绝。只有显式
`--execute --overwrite` 才会替换上述扩展输出，并先存入
`lensedcmb_new_backups/`；原有 CMB 文件保持原位。

## 先查看计划，再生成

以下集群数据路径是输入示例，可替换为具有同一五位 ID 目录结构的仓库：

```bash
CMB_DATA_ROOT=/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_ILC/mmDL_maps_2048
./run_cmb.sh --mmdl-root "$CMB_DATA_ROOT" \
  --theory "$PWD/theory/outputs/sehgal2010_theory.npz" \
  --ids 0,1 --master-seed 20261011 \
  --manifest "$PWD/evidence/example_plan.json"
```

默认只读检查和写出计划。正式生成需在已分配的计算节点使用 `--execute`，
并按分配的 CPU 数指定线程：

```bash
./run_cmb.sh --mmdl-root "$CMB_DATA_ROOT" \
  --theory "$PWD/theory/outputs/sehgal2010_theory.npz" \
  --ids 0 --master-seed 20261011 --threads 16 --execute \
  --manifest "$PWD/evidence/example_production.json"
```

也支持 `--start-id 0 --count 10` 和 `--ids-file /path/ids.txt`。
全局计划/结果 `--manifest` 必须写在输入数据树外。

## Slurm 工作流

生产按 ID 独立运行，不需要 MPI。下面是可在所选安装路径运行的十套作业模板；
每个数组任务使用一个进程和16CPU，排除 aliws001、003、004：

```bash
export CMB_RELEASE_ROOT=$(pwd -P)
export CMB_DATA_ROOT=/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_ILC/mmDL_maps_2048
mkdir -p logs evidence
cat > cmb_job.slurm <<'SLURM'
#!/usr/bin/env bash
#SBATCH --partition=ali
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --array=0-9%10
#SBATCH --exclude=aliws001,aliws003,aliws004
set -euo pipefail
: "${CMB_RELEASE_ROOT:?}" "${CMB_DATA_ROOT:?}"
exec "${CMB_RELEASE_ROOT}/run_cmb.sh" \
  --mmdl-root "${CMB_DATA_ROOT}" \
  --theory "${CMB_RELEASE_ROOT}/theory/outputs/sehgal2010_theory.npz" \
  --ids "${SLURM_ARRAY_TASK_ID}" --threads "${SLURM_CPUS_PER_TASK}" \
  --master-seed 20261011 --owner-coordinate G --execute \
  --manifest "${CMB_RELEASE_ROOT}/evidence/production_${SLURM_ARRAY_TASK_ID}.json"
SLURM
sbatch --export=ALL,CMB_RELEASE_ROOT,CMB_DATA_ROOT \
  --output="$CMB_RELEASE_ROOT/logs/cmb_%A_%a.log" cmb_job.slurm
```

批量其余 ID 时可用 `sbatch --array=10-499%48` 覆盖数组范围和并发上限，
其余参数同上。`slurm/` 保存实际集群作业脚本，固定为已经验证的安装和数据
绝对路径；迁移安装目录时，应采用上面的模板或同步修改这些路径与日志路径。
已提供的成功理论无需重算；另算理论时应选择新的输出目录以保留发布参考。

## 核验与试生成

环境安装后可运行轻量约定/拒绝路径测试：

```bash
source ./_cmb_env.sh
"$CMB_RUNTIME/bin/python" -m unittest discover -s theory -p 'test_*.py'
"$CMB_RUNTIME/bin/python" -m unittest discover -s production -p 'test_*.py'
"$CMB_RUNTIME/bin/python" -m unittest discover -s validation -p 'test_*.py'
"$CMB_RUNTIME/bin/python" validation/validate_inputs.py --root . \
  --report evidence/local_theory_input_audit.json
```

这些小数组测试不是降采样的原生图验证。实际原生十套 ID0..9 试生成与理论计算
已有作业记录。理论独立复核见
`theory/source_evidence/independent_theory_review.json`；试生成的原生谱、实测
κ 两点参考、图和分箱汇总保存在 `validation/evidence/pilot_review/`。
原生测量与汇总的完整命令见 `validation/design.md`。TT/EE/TE 的预先声明检查
和 BB 的模型比较分别报告，完成计算不等于对原作者原始配置或整个数据集的
科学认证；没有为匹配 BB 调整谱幅度。

十套平均相对自洽 CAMB 的拟合幅度偏差为 TT +0.0471%、EE −0.0525%、
TE +0.2088%、BB −0.6662%。BB 的亚百分比差异仍有统计显著性；以实际
十套 κ 平均二点谱作条件透镜参考时，BB 偏差为 −0.0395%。这些拟合仅用于
验证，没有用于修正生成天图。详见 [定量结果](validation/RESULTS.md) 和
[十套平均与样本标准差图](validation/figures/native_cmb_pilot_10.pdf)。
全部500套的文件契约及原图保存检查见 `evidence/campaign_status_final.json`；
该完整性检查与仅覆盖十套、ℓ<1000 的角功率谱检查范围不同。
