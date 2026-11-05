import json
import urllib.parse
import urllib.request
from Bio import Entrez, SeqIO

Entrez.email = "19065218124@163.com"


# ==============================================================
# 步骤 1 & 2: 通过 STRING-db API 检索互作 B 蛋白
# 同时提取其标准 AGI 位点编号 (例如 3702.AT4G02640.1 -> AT4G02640)
# ==============================================================
def get_interacting_partners(gene_name="BZIP53", species_taxid=3702, top_n=5):
    string_api_url = "https://string-db.org/api/json/interaction_partners"
    params = {
        "identifiers": gene_name,
        "species": species_taxid,
        "limit": top_n,
        "caller_identity": "bio_research_lab"
    }

    url = f"{string_api_url}?{urllib.parse.urlencode(params)}"
    print(f"[步骤 1 & 2] 正在从 STRING 检索诱饵蛋白 {gene_name} 的互作网络...")

    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = json.loads(resp.read().decode('utf-8'))

    partners = []
    for item in data:
        # 判断谁是伴侣蛋白，同时提取标准 locus_tag
        if item['preferredName_A'].lower() == gene_name.lower():
            partner = item['preferredName_B']
            string_id = item['stringId_B']
        else:
            partner = item['preferredName_A']
            string_id = item['stringId_A']

        # 从 "3702.AT4G02640.1" 中提取纯位点号 "AT4G02640"
        locus = string_id.split('.')[1] if '.' in string_id else string_id
        score = item['score']
        partners.append((partner, locus, score))

    # 去重
    seen = set()
    unique_partners = []
    for p, loc, s in sorted(partners, key=lambda x: x[2], reverse=True):
        if p not in seen and p.lower() != gene_name.lower():
            seen.add(p)
            unique_partners.append((p, loc, s))

    return unique_partners


# ==============================================================
# 步骤 3 & 4: 依据 AGI 位点编号从 NCBI 精确获取转录本与 CDS
# ==============================================================
def fetch_cds_and_protein(gene_name, locus_tag, organism="Arabidopsis thaliana"):
    print(f"\n[步骤 3 & 4] 正在 NCBI 检索互作 B 蛋白 ({gene_name} / {locus_tag}) 的 mRNA 转录本...")

    # 使用标准 AGI locus_tag 检索，避免命名缩写差异
    search_queries = [
        f'"{locus_tag}"[Locus Tag] AND "{organism}"[Organism] AND biomol_mrna[Properties]',
        f'"{locus_tag}"[All Fields] AND "{organism}"[Organism] AND biomol_mrna[Properties]',
        f'"{gene_name}"[All Fields] AND "{organism}"[Organism] AND biomol_mrna[Properties]'
    ]

    gi_id = None
    for q in search_queries:
        search_handle = Entrez.esearch(db="nuccore", term=q, retmax=3, sort="relevance")
        search_results = Entrez.read(search_handle)
        search_handle.close()
        if search_results["IdList"]:
            gi_id = search_results["IdList"][0]
            break

    if not gi_id:
        print(f"未能自动检索到 {gene_name} ({locus_tag}) 的转录本。")
        return None

    # 获取 GenBank 记录
    fetch_handle = Entrez.efetch(db="nuccore", id=gi_id, rettype="gb", retmode="text")
    record = SeqIO.read(fetch_handle, "genbank")
    fetch_handle.close()

    # 提取 CDS
    cds_feat = next((f for f in record.features if f.type == "CDS"), None)
    if not cds_feat:
        print(f"在记录 {record.id} 中未找到 CDS 注释。")
        return None

    cds_seq = cds_feat.extract(record.seq)
    protein_seq = cds_seq.translate(to_stop=True)
    protein_id = cds_feat.qualifiers.get("protein_id", ["Unknown"])[0]

    return {
        "accession": record.id,
        "gene_name": gene_name,
        "locus_tag": locus_tag,
        "protein_id": protein_id,
        "location": str(cds_feat.location),
        "cds_seq": str(cds_seq),
        "protein_seq": str(protein_seq)
    }


if __name__ == "__main__":
    bait_gene = "BZIP53"
    partners = get_interacting_partners(bait_gene, top_n=5)

    print("\n" + "=" * 60)
    print(f"与诱饵蛋白 A ({bait_gene}) 存在强互作的候选 B 蛋白：")
    for idx, (p, loc, score) in enumerate(partners, 1):
        print(f"  {idx}. 基因: {p:<8} 位点号: {loc:<12} (STRING得分: {score})")
    print("=" * 60)

    if partners:
        chosen_b_gene, chosen_locus, _ = partners[0]
        print(f"\n>>> 选定构建猎物载体的目标 B 蛋白: {chosen_b_gene} ({chosen_locus})")

        res = fetch_cds_and_protein(chosen_b_gene, chosen_locus)

        if res:
            print("\n" + "=" * 65)
            print(f"【成功提取 B 蛋白 ({chosen_b_gene}) 载体构建所需序列】")
            print(f"• 转录本 Accession: {res['accession']}")
            print(f"• 基因位点编号:     {res['locus_tag']}")
            print(f"• 蛋白质编号:       {res['protein_id']}")
            print(f"• CDS 注释区间:     {res['location']}")
            print(f"• CDS 碱基长度:     {len(res['cds_seq'])} bp")
            print(f"• 翻译蛋白长度:     {len(res['protein_seq'])} aa")
            print("=" * 65)

            print(f"\n>CDS_{chosen_b_gene}_for_Prey_Vector length={len(res['cds_seq'])}bp")
            print(res['cds_seq'])

            print(f"\n>Protein_{chosen_b_gene} length={len(res['protein_seq'])}aa")
            print(res['protein_seq'])