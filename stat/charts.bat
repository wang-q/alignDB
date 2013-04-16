echo "Draw charts for all .xlsx files"

echo "common"
forfiles /m *.common.xlsx -c "cmd /c perl d:/wq/Scripts/alignDB/stat/common_chart_factory.pl -i @path"

echo "gc"
forfiles /m *.gc.xlsx -c "cmd /c perl d:/wq/Scripts/alignDB/stat/gc_chart_factory.pl -i @path"

echo "gene"
forfiles /m *.gene.xlsx -c "cmd /c perl d:/wq/Scripts/alignDB/stat/gene_chart_factory.pl -i @path"

echo "multi"
forfiles /m *.multi.xlsx -c "cmd /c perl d:/wq/Scripts/alignDB/stat/multi_chart_factory.pl -i @path"

echo "dnds"
forfiles /m *.dnds.xlsx -c "cmd /c perl d:/wq/Scripts/alignDB/stat/dnds_chart_factory.pl -i @path"

echo "ofg"
forfiles /m *.ofg.xlsx -c "cmd /c perl d:/wq/Scripts/alignDB/stat/ofg_chart_factory.pl -i @path"
