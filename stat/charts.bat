echo "Draw charts for all .xls files"

echo "common"
forfiles /m *.common.xls -c "cmd /c perl e:/wq/Scripts/alignDB/stat/common_chart_factory.pl -i @path"

echo "gc"
forfiles /m *.gc.xls -c "cmd /c perl e:/wq/Scripts/alignDB/stat/gc_chart_factory.pl -i @path"

echo "gene"
forfiles /m *.gene.xls -c "cmd /c perl e:/wq/Scripts/alignDB/stat/gene_chart_factory.pl -i @path"

