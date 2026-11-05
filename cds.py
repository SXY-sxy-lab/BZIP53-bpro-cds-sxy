from Bio import Entrez, SeqIO

# 1. 设置 NCBI 访问邮箱（NCBI 官方规则要求）
Entrez.email = "19065218124@163.com"

# 2. 从 nuccore 数据库在线获取 NM_116107.2 的 GenBank 格式完整记录
print("正在从 NCBI 数据库抓取 NM_116107.2 记录...")
handle = Entrez.efetch(db="nuccore", id="NM_116107.2", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()

# 3. 遍历特征表（Feature Table），定位并提取 CDS
cds_feature = None
for feature in record.features:
    if feature.type == "CDS":
        cds_feature = feature
        break

if cds_feature:
    # 提取核酸序列
    cds_sequence = cds_feature.extract(record.seq)
    gene_name = cds_feature.qualifiers.get("gene", ["Unknown"])[0]
    protein_id = cds_feature.qualifiers.get("protein_id", ["Unknown"])[0]

    print("\n" + "=" * 60)
    print("【CDS 序列提取成功】")
    print(f"目标基因: {gene_name}")
    print(f"蛋白质编号: {protein_id}")
    print(f"注释位置 (Location): {cds_feature.location}")
    print(f"序列全长: {len(cds_sequence)} bp")
    print("=" * 60 + "\n")

    # 打印 FASTA 格式
    print(f">NM_116107.2_CDS gene={gene_name} location=567..1007 length={len(cds_sequence)}bp")
    print(str(cds_sequence))
else:
    print("未在记录中找到 CDS 注释。")