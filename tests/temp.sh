# 查找 Viruses 的 TaxID 并筛选提取下属的所有 Kingdom
taxonkit list --ids 10239 --data-dir /data/mengxf/Database/NCBI/taxonomy \
  | taxonkit lineage --data-dir /data/mengxf/Database/NCBI/taxonomy \
  | taxonkit reformat -f "{K}" -F --data-dir /data/mengxf/Database/NCBI/taxonomy \
  | cut -f 3 \
  | sort -u \
  | grep -v "^$"
