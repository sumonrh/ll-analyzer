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
    $result = $excel.Run("'" + $book.Name + "'!RunTruckSpacingTests")
    if ($result -ne 'PASS') { throw $result }
    $result = $excel.Run("'" + $book.Name + "'!TestCustomStorage")
    if ($result -ne 'PASS') { throw $result }
    $eventsModule = $book.VBProject.VBComponents.Item($book.CodeName)
    $events = Get-Content -LiteralPath (Join-Path $root 'LL Analysis Workbook Events.txt') -Raw
    $eventsModule.CodeModule.AddFromString(($events -replace '\bMsgBox\b', 'TestMsgBox'))
    $result = $excel.Run("'" + $book.Name + "'!TestOptionalWorkbookEvents")
    if ($result -ne 'PASS') { throw $result }
    Write-Output 'PASS: All three BCL subdivisions and endpoints, spacing envelopes/FEM reconstruction, dropdown validation, editable presets with extra-axle warnings, Custom restoration and 1-20 axle analysis.'
} finally {
    if ($book) { $book.Close($false) }
    if ($excel) {
        $excel.StatusBar = $false
        $excel.Quit()
    }
    if ($module) { [void][Runtime.InteropServices.Marshal]::FinalReleaseComObject($module) }
    if ($eventsModule) { [void][Runtime.InteropServices.Marshal]::FinalReleaseComObject($eventsModule) }
    if ($book) { [void][Runtime.InteropServices.Marshal]::FinalReleaseComObject($book) }
    if ($excel) { [void][Runtime.InteropServices.Marshal]::FinalReleaseComObject($excel) }
}
