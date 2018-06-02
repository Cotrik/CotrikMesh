folders=`ls -l | grep "^d" | awk '{print$9}'`
> index.info
sheetPdfList=""
dominantSheetPdfList=""
n=1
for folder in $folders;
do
    sheetPdfList+="$folder/SheetsConnectivities.pdf "
    dominantSheetPdfList+="$folder/DominantSheetsConnectivities.pdf "
    cd $folder
    ExtractSheetDecompositions $folder.vtk
    Rscript ../correlation_sheet.r SheetsConnectivities.mat SheetsConnectivities.pdf
    Rscript ../correlation_sheet.r DominantSheetsConnectivities.mat DominantSheetsConnectivities.pdf
    echo "[/Page $n /Title ($folder) /OUT pdfmark" >> ../index.info
    cd ..
    let n=n+1
done

pdfunite $sheetPdfList CombinedSheetsConnectivities.pdf
pdfunite $dominantSheetPdfList CombinedDominantSheetsConnectivities.pdf
gs -sDEVICE=pdfwrite -q -dBATCH -dNOPAUSE -sOutputFile=CombinedSheetsConnectivitiesWithBookMark.pdf -dPDFSETTINGS=/prepress index.info -f CombinedSheetsConnectivities.pdf
gs -sDEVICE=pdfwrite -q -dBATCH -dNOPAUSE -sOutputFile=CombinedDominantSheetsConnectivitiesWithBookMark.pdf -dPDFSETTINGS=/prepress index.info -f CombinedDominantSheetsConnectivities.pdf
