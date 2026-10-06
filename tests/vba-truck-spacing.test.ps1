$ErrorActionPreference = 'Stop'
$root = Split-Path $PSScriptRoot -Parent
$excel = $null
$book = $null
try {
    $excel = New-Object -ComObject Excel.Application
    $excel.Visible = $false
    $excel.DisplayAlerts = $false
    $book = $excel.Workbooks.Add()
    $module = $book.VBProject.VBComponents.Add(1)
    $module.Name = 'TruckSpacingTests'
    $source = Get-Content -LiteralPath (Join-Path $root 'LL Analysis VBA Code.txt') -Raw
    # Suppress only test-workbook message boxes; critical messages become test failures.
    $source = $source -replace '\bMsgBox\b', 'TestMsgBox'
    $harness = Get-Content -LiteralPath (Join-Path $PSScriptRoot 'vba-truck-spacing-harness.bas') -Raw
    $module.CodeModule.AddFromString($source + "`r`n" + $harness)
    $excel.Run("'" + $book.Name + "'!TestTruckSpacing")
    Write-Output 'PASS: Excel VBA spacing enumeration, fixed-truck envelope parity, governing FEM reconstruction, and input/output integration.'
} finally {
    if ($book) { $book.Close($false) }
    if ($excel) {
        $excel.StatusBar = $false
        $excel.Quit()
    }
    if ($module) { [void][Runtime.InteropServices.Marshal]::FinalReleaseComObject($module) }
    if ($book) { [void][Runtime.InteropServices.Marshal]::FinalReleaseComObject($book) }
    if ($excel) { [void][Runtime.InteropServices.Marshal]::FinalReleaseComObject($excel) }
}
