Private Function TestExtraAxleWarnings(Optional increment As Boolean = False, _
                                      Optional reset As Boolean = False) As Long
    Static count As Long
    If reset Then count = 0
    If increment Then count = count + 1
    TestExtraAxleWarnings = count
End Function

Public Sub TestMsgBox(ByVal prompt As String, Optional ByVal buttons As VbMsgBoxStyle = vbOKOnly)
    If InStr(prompt, "will be ignored and erased") > 0 Then
        TestExtraAxleWarnings True
        TestAssert Application.WorksheetFunction.CountA( _
            ThisWorkbook.Worksheets("LL Input").Range("B15:C29")) > 0, "Warn before erasing extra axles"
    End If
    If (buttons And &H70) = vbCritical Then Err.Raise vbObjectError + 100, "VBA test", prompt
End Sub

Public Function RunTruckSpacingTests() As String
    On Error GoTo TestFailed
    TestTruckSpacing
    RunTruckSpacingTests = "PASS"
    Exit Function
TestFailed:
    Application.EnableEvents = True
    RunTruckSpacingTests = "FAIL: " & Err.Source & ": " & Err.Description
End Function

Private Sub TestAssert(condition As Boolean, message As String)
    If Not condition Then Err.Raise vbObjectError + 101, "VBA test", message
End Sub

Private Sub TestSelectModel(ws As Worksheet, index As Long)
    Dim selector As DropDown
    Set selector = ws.DropDowns(TRUCK_MODEL_DROPDOWN)
    selector.ListIndex = index
    Application.Run selector.OnAction
End Sub

Private Sub TestClose(ByVal actual As Double, ByVal expected As Double, message As String)
    TestAssert Abs(actual - expected) <= 0.0000001 * (1# + Abs(expected)), _
        message & ": " & actual & " <> " & expected
End Sub

Private Sub TestPair(values() As Double, ByRef count As Long, maximum As Double, minimum As Double)
    values(count) = maximum: values(count + 1) = minimum: count = count + 2
End Sub

Private Function TestCaseVector(label As String, variableTruck As Boolean, gap As Double, _
    automatic As Boolean, multiplier As Double, factor As Double, udl As Double) As Double()
    Dim spans(2) As Double, elements(5) As Double, nodes(6) As Double, supportX(3) As Double
    Dim supports(3) As Long, constrained(13) As Boolean, K() As Double
    Dim axles(4) As Double, spacings(4) As Double, steps(72) As Double
    Dim maxV() As Double, minV() As Double, maxM() As Double, minM() As Double
    Dim maxD() As Double, minD() As Double, maxR() As Double, minR() As Double
    Dim values() As Double, count As Long, i As Long, s As Long, row As Long, r As Long
    Dim loads As Variant, gaps As Variant, rebuilt As Double
    spans(0) = 12#: spans(1) = 18#: spans(2) = 12#
    loads = Array(50#, 140#, 140#, 175#, 120#)
    gaps = Array(3.6, 1.2, gap, 6.6, 0#)
    For i = 0 To 4
        axles(i) = loads(i): spacings(i) = gaps(i)
    Next i
    For i = 0 To 5
        elements(i) = spans(i \ 2) / 2#
        nodes(i + 1) = nodes(i) + elements(i)
    Next i
    For s = 0 To 3
        supports(s) = s * 4: constrained(supports(s)) = True
        supportX(s) = nodes(s * 2)
    Next s
    For i = 0 To 72: steps(i) = -29.4 + i * 1.4: Next i
    ReDim K(13, 13)
    BuildStiffnessMatrix K, 6, elements, 200000000000#, 0.005
    FactorizeSystem K, constrained, 14
    m_IsBCLTruck = variableTruck
    RunLoadCaseEnvelope label, factor, udl, axles, spacings, 5, _
        3, spans, 42#, 6, elements, 14, 7, 12, 200000000000#, 0.005, _
        K, supports, 4, supportX, steps, 73, _
        maxV, minV, maxM, minM, maxD, minD, maxR, minR, automatic, multiplier, True
    ReDim values(2 * (12 + 7 + 7 + 73 * 4 + 4) - 1)
    For i = 0 To 11: TestPair values, count, maxV(i), minV(i): Next i
    For i = 0 To 6: TestPair values, count, maxM(i), minM(i): Next i
    For i = 0 To 6: TestPair values, count, maxD(i), minD(i): Next i
    For i = 0 To 72
        For s = 0 To 3: TestPair values, count, maxR(i, s), minR(i, s): Next s
    Next i
    For s = 0 To 3
        TestPair values, count, m_OptimizedSupport(0, s), m_OptimizedSupport(1, s)
    Next s
    If label = "Truck" Then
        For r = 0 To UBound(m_TruckGovValue, 2)
            For row = 0 To 1
                If variableTruck Then
                    TestAssert m_TruckGovSpacing(row, r) >= 6.6 And m_TruckGovSpacing(row, r) <= 18#, _
                        "Governing gap outside inclusive range"
                End If
                rebuilt = ReconstructedTruckResponse(r, row = 0, axles, spacings, 5, elements, 6, _
                    14, 200000000000#, 0.005, K, supports, 4, automatic, multiplier, factor)
                TestClose rebuilt, m_TruckGovValue(row, r), "Governing response reconstruction"
            Next row
        Next r
    End If
    TestCaseVector = values
End Function

Private Sub TestSpacingEnvelope(label As String, automatic As Boolean, multiplier As Double, _
    factor As Double, udl As Double)
    Dim actual() As Double, fixed() As Double, expected() As Double
    Dim index As Long, i As Long, gap As Double
    actual = TestCaseVector(label, True, 6.6, automatic, multiplier, factor, udl)
    ReDim expected(UBound(actual))
    Dim lastIndex As Long
    lastIndex = Int(11.4 / m_BCLSpacingIncrement) + 1
    For index = 0 To lastIndex
        gap = 6.6 + m_BCLSpacingIncrement * index
        If index = lastIndex Then gap = 18#
        fixed = TestCaseVector(label, False, gap, automatic, multiplier, factor, udl)
        For i = 0 To UBound(actual) Step 2
            If fixed(i) > expected(i) Then expected(i) = fixed(i)
            If fixed(i + 1) < expected(i + 1) Then expected(i + 1) = fixed(i + 1)
        Next i
    Next index
    For i = 0 To UBound(actual)
        TestClose actual(i), expected(i), label & " envelope ordinate " & i
    Next i
End Sub

Sub TestTruckSpacing()
    Dim i As Long, ws As Worksheet, results As Worksheet, heading As Range
    Dim model As Boolean, errorNumber As Long
    Dim increments As Variant, counts As Variant, subdivisionIndex As Long
    increments = Array(0.5, 1#, 2#): counts = Array(24, 13, 7)
    For subdivisionIndex = 0 To 2
    m_BCLSpacingIncrement = increments(subdivisionIndex)
    m_IsBCLTruck = True
    TestAssert TruckSpacingCount() = counts(subdivisionIndex), "Expected subdivision configuration count"
    For i = 0 To counts(subdivisionIndex) - 2
        TestClose TruckSpacingAt(i), 6.6 + i * increments(subdivisionIndex), "Spacing " & i
    Next i
    TestClose TruckSpacingAt(counts(subdivisionIndex) - 1), 18#, "Inclusive final spacing"
    m_IsBCLTruck = False
    TestAssert TruckSpacingCount() = 1, "CL must have one geometry"
    TestSpacingEnvelope "Truck", True, 0.75, 1#, 0#
    TestSpacingEnvelope "Truck", False, 1#, 1.25, 0#
    TestSpacingEnvelope "Lane", False, 1#, 0.8, 9#
    TestSpacingEnvelope "Lane", False, 1#, 0.8, 0#
    Next subdivisionIndex

    CreateInputSheet
    Set ws = ThisWorkbook.Worksheets("LL Input")
    TestAssert ws.Range("B8").Value = "CL-625", "New sheets default to CL"
    ws.Range("B8").ClearContents
    TestAssert Not ReadTruckModel(ws), "Blank truck model means CL"
    ws.Protect
    EnsureInputLayout ws
    TestAssert ws.Range("B8").Value = "CL-625", "Existing sheets gain CL default"
    TestAssert ws.Range("B8").Validation.Formula1 = "CL-625,BCL-625,Custom", "Selector order"
    TestAssert ws.DropDowns.Count = 2, "Setup must not duplicate either dropdown"
    TestBCLSubdivisionInput ws
    TestAssert ws.DropDowns(TRUCK_MODEL_DROPDOWN).ListIndex = 1, "Dropdown must match CL cell value"
    TestAssert Not ws.ProtectContents And Not ws.ProtectDrawingObjects, "Input worksheet must remain unprotected"
    TestAssert Not ws.Range("B15:C29").Locked, "Preset rows must allow formatting and editing"
    TestAssert Not ws.Range("B8").Locked And Not ws.Range("B41").Locked, "Other inputs remain editable"
    Application.EnableEvents = False
    ws.Range("B8").Value = "invalid truck"
    On Error Resume Next
    model = ReadTruckModel(ws)
    errorNumber = Err.Number
    Err.Clear
    On Error GoTo 0
    TestAssert errorNumber = vbObjectError + 26, "Invalid truck model must raise an explicit error"
    ws.Range("B8").Value = CVErr(xlErrValue)
    On Error Resume Next
    model = ReadTruckModel(ws)
    errorNumber = Err.Number
    Err.Clear
    On Error GoTo 0
    TestAssert errorNumber = vbObjectError + 26, "Worksheet-error truck model must be rejected"
    Application.EnableEvents = True
    TestSelectModel ws, 2
    TestClose ws.Range("B11").Value, 140#, "Live BCL preset load"
    TestClose ws.Range("B14").Value, 120#, "Live BCL final axle load"
    TestClose ws.Range("C12").Value, 6.6, "Displayed variable gap minimum"
    TestAssert Application.WorksheetFunction.CountA(ws.Range("B15:C29")) = 0, "Inactive rows must be cleared"
    ws.Range("B32").Value = 2
    ws.Range("B33").Value = 2#
    ws.Range("B7").Value = 3
    ws.Range("B15").Value = 99#
    ws.Range("C29").Value = 7#
    EnsureInputLayout ws
    TestClose ws.Range("B15").Value, 99#, "Setup must not silently erase extra axles"
    TestExtraAxleWarnings False, True
    MainAnalysis
    TestAssert TestExtraAxleWarnings() = 1, "Analyze must warn once for BCL extra loads/spacings"
    TestAssert Application.WorksheetFunction.CountA(ws.Range("B15:C29")) = 0, "Analyze must erase BCL extra axles"
    Set results = ThisWorkbook.Worksheets("LL Results")
    TestAssert InStr(CStr(results.Range("A1").Value), "BCL-625") > 0, "Results identify truck"
    TestAssert InStr(CStr(results.Range("A1").Value), "13 configurations") > 0, "Default subdivision must reach results"
    Set heading = results.Cells.Find("Governing axle 3-4 spacing (m)", LookIn:=xlValues, LookAt:=xlWhole)
    TestAssert Not heading Is Nothing, "Governing spacing column missing"
    For i = 1 To 4
        TestAssert heading.Offset(i, 0).Value >= 6.6 And heading.Offset(i, 0).Value <= 18#, _
            "Diagnostic governing gap outside inclusive range"
    Next i
    ws.DropDowns(BCL_INCREMENT_DROPDOWN).ListIndex = 3
    Application.Run ws.DropDowns(BCL_INCREMENT_DROPDOWN).OnAction
    MainAnalysis
    TestAssert InStr(CStr(results.Range("A1").Value), "7 configurations") > 0, "2 m subdivision must reach analysis/results"
    TestClose m_BCLSpacingIncrement, 2#, "Analysis reads selected subdivision"
    ws.DropDowns(BCL_INCREMENT_DROPDOWN).ListIndex = 2
    Application.Run ws.DropDowns(BCL_INCREMENT_DROPDOWN).OnAction
    TestSelectModel ws, 1
    TestClose ws.Range("B11").Value, 125#, "Live CL preset load"
    TestClose ws.Range("B14").Value, 150#, "Live CL final axle load"
    ws.Range("C15").Value = 4#
    TestExtraAxleWarnings False, True
    MainAnalysis
    TestAssert TestExtraAxleWarnings() = 1, "Spacing-only extra input must warn for CL"
    TestAssert Application.WorksheetFunction.CountA(ws.Range("B15:C29")) = 0, "Analyze must erase CL extra axles"
    TestAssert Not m_IsBCLTruck, "Switching back must restore fixed CL behavior"
    TestAssert InStr(CStr(results.Range("A1").Value), "BCL-625") = 0, "Stale BCL result label"
    TestCustomInputs ws
    Application.StatusBar = False
End Sub

Private Sub TestBCLSubdivisionInput(ws As Worksheet)
    Dim selector As DropDown, i As Long, increment As Double, errorNumber As Long, invalid As Variant
    TestClose ws.Range("B37").Value, 1#, "Default BCL subdivision is 1 m"
    Set selector = ws.DropDowns(BCL_INCREMENT_DROPDOWN)
    TestAssert selector.ListIndex = 2 And selector.ListCount = 3, "Subdivision dropdown default and options"
    TestAssert ws.Range("B37").Validation.Type = xlValidateCustom, "Manual subdivision entry requires validation"
    TestAssert Not ws.Range("B37").Validation.Value, "Manual typing must be rejected"
    TestAssert ws.Range("B37").Validation.AlertStyle = xlValidAlertStop, "Invalid entry must stop"
    For i = 1 To 3
        selector.ListIndex = i
        Application.Run selector.OnAction
        Select Case i
            Case 1: increment = 0.5
            Case 2: increment = 1#
            Case 3: increment = 2#
        End Select
        TestClose ReadBCLSubdivision(ws), increment, "Subdivision callback numeric value"
        EnsureInputLayout ws
        TestAssert selector.ListIndex = i, "Setup must preserve selected subdivision"
    Next i
    For Each invalid In Array(0.001, 0#, -1#, 3#, "bad", True, CVErr(xlErrValue))
        ws.Range("B37").Value = invalid
        On Error Resume Next
        increment = ReadBCLSubdivision(ws)
        errorNumber = Err.Number
        Err.Clear
        On Error GoTo 0
        TestAssert errorNumber = vbObjectError + 31, "Pasted or programmatic invalid subdivision must fail"
    Next invalid
    selector.ListIndex = 2
    Application.Run selector.OnAction
End Sub

Private Sub TestCustomInputs(ws As Worksheet)
    Dim i As Long, axles() As Double, spacings() As Double, count As Long, errorNumber As Long
    Dim cache As Worksheet, expected() As Variant
    Application.EnableEvents = False
    TestSelectModel ws, 3
    TestAssert Not Application.EnableEvents, "Dropdown action must preserve disabled events"
    Application.EnableEvents = True
    TestAssert Not ws.Range("B10:C29").Locked, "Custom must unlock all 20 rows"
    TestAssert ws.Range("B29").Interior.Color = RGB(230, 240, 255), "Last custom load row must be active immediately"
    TestAssert ws.Range("C29").Interior.Color = RGB(230, 240, 255), "Last custom spacing row must be active immediately"
    For i = 0 To 19
        ws.Cells(10 + i, 2).Value = 50# + i
        ws.Cells(10 + i, 3).Value = 0.2 + i * 0.1
    Next i
    expected = ws.Range("B10:C29").Value
    ReadCustomAxles ws, axles, spacings, count
    TestAssert count = 20, "Custom must read all 20 axles"
    TestSelectModel ws, 2
    TestAssert Application.WorksheetFunction.CountA(ws.Range("B15:C29")) = 0, "Preset clears all extra custom axles"
    TestSelectModel ws, 1
    TestSelectModel ws, 3
    For i = 0 To 19
        TestClose ws.Cells(10 + i, 2).Value, expected(i + 1, 1), "Restored custom load"
        TestClose ws.Cells(10 + i, 3).Value, expected(i + 1, 2), "Restored custom spacing"
    Next i
    Set cache = ThisWorkbook.Worksheets(CUSTOM_AXLE_SHEET)
    TestAssert cache.Visible = xlSheetVeryHidden, "Custom cache must be hidden"
    MainAnalysis
    TestAssert Not m_IsBCLTruck, "Custom must not use BCL spacing sweep"
    TestAssert TestExtraAxleWarnings() = 1, "Custom must not issue preset extra-axle warnings"
    TestClose ws.Range("B29").Value, expected(20, 1), "Custom Analyze must retain axle 20"
    TestClose ws.Range("C29").Value, expected(20, 2), "Custom Analyze must retain final gap"
    ws.Range("B10:C29").ClearContents
    ws.Range("B10").Value = 100#
    ws.Range("C10").Value = "ignored last gap"
    ReadCustomAxles ws, axles, spacings, count
    TestAssert count = 1, "Single custom axle must be supported"
    MainAnalysis
    ws.Range("B10:C29").ClearContents
    ws.Range("B10:B11").Value = Application.Transpose(Array(50#, 80#))
    ws.Range("C10").Value = 0#
    ReadCustomAxles ws, axles, spacings, count
    TestAssert count = 2 And spacings(0) = 0#, "Coincident custom axles must be supported"
    For i = 0 To 3
        ws.Range("B10:C29").ClearContents
        ws.Range("B10:B11").Value = 50#
        ws.Range("C10").Value = 1#
        Select Case i
            Case 0: ws.Range("B11").ClearContents: ws.Range("B12").Value = 50#
            Case 1: ws.Range("B11").Value = "bad load"
            Case 2: ws.Range("C10").Value = -1#
            Case 3: ws.Range("C10").ClearContents
        End Select
        On Error Resume Next
        ReadCustomAxles ws, axles, spacings, count
        errorNumber = Err.Number
        Err.Clear
        On Error GoTo 0
        TestAssert errorNumber = vbObjectError + 29, "Invalid custom input must fail explicitly"
    Next i
    Application.EnableEvents = False
    ws.Range("B8").Value = "CL-625"
    Application.EnableEvents = True
    EnsureInputLayout ws
    TestClose ws.Range("B11").Value, 125#, "Analyze fallback must synchronize preset"
    TestAssert Application.EnableEvents, "Events must remain enabled"
    TestSelectModel ws, 3
    For i = 0 To 19
        ws.Cells(10 + i, 2).Value = 50# + i
        ws.Cells(10 + i, 3).Value = 10# + i
    Next i
    MainAnalysis
End Sub

Public Function TestCustomStorage() As String
    Dim ws As Worksheet, cache As Worksheet, i As Long
    On Error GoTo TestFailed
    Set ws = ThisWorkbook.Worksheets("LL Input")
    Set cache = ThisWorkbook.Worksheets(CUSTOM_AXLE_SHEET)
    TestAssert ws.Range("B8").Value = "Custom", "Custom model must remain selected"
    TestAssert Not ws.ProtectContents, "Custom worksheet must remain unprotected"
    TestAssert Not ws.Range("B10:C29").Locked, "Custom rows must remain editable"
    TestAssert cache.Visible = xlSheetVeryHidden, "Custom cache must remain hidden"
    TestSelectModel ws, 2
    TestClose ws.Range("B11").Value, 140#, "Workbook events must fill BCL preset"
    For i = 0 To 19
        TestClose cache.Cells(10 + i, 2).Value, 50# + i, "Worksheet-backed custom load"
        TestClose cache.Cells(10 + i, 3).Value, 10# + i, "Worksheet-backed custom gap"
    Next i
    TestSelectModel ws, 3
    For i = 0 To 19
        TestClose ws.Cells(10 + i, 2).Value, 50# + i, "Saved custom load"
        TestClose ws.Cells(10 + i, 3).Value, 10# + i, "Saved custom gap"
    Next i
    TestCustomStorage = "PASS"
    Exit Function
TestFailed:
    TestCustomStorage = "FAIL: " & Err.Source & ": " & Err.Description
End Function

Public Function TestOptionalWorkbookEvents() As String
    Dim ws As Worksheet
    On Error GoTo TestFailed
    Set ws = ThisWorkbook.Worksheets("LL Input")
    ws.Range("B8").Value = "BCL-625"
    TestClose ws.Range("B11").Value, 140#, "Optional cell edit updates preset immediately"
    TestAssert ws.DropDowns(TRUCK_MODEL_DROPDOWN).ListIndex = 2, "Cell edits must synchronize dropdown"
    ws.Range("B8").Value = "Custom"
    TestAssert Not ws.Range("B10:C29").Locked, "Optional cell edit unlocks Custom immediately"
    TestAssert ws.DropDowns(TRUCK_MODEL_DROPDOWN).ListIndex = 3, "Custom cell edit must synchronize dropdown"
    TestSelectModel ws, 1
    TestClose ws.Range("B11").Value, 125#, "Form Control works with optional events installed"
    TestAssert Not ws.ProtectContents, "Preset switching must not reprotect the worksheet"
    ws.Range("B10:C29").Font.Name = "Arial"
    ws.Range("B10:C29").Font.Size = 12
    ws.Range("B10:C29").Borders.LineStyle = xlContinuous
    TestClose ws.Range("B29").Font.Size, 12#, "Preset table fonts must be editable"
    TestAssert Application.EnableEvents, "Optional handlers must preserve enabled events"
    TestOptionalWorkbookEvents = "PASS"
    Exit Function
TestFailed:
    TestOptionalWorkbookEvents = "FAIL: " & Err.Source & ": " & Err.Description
End Function
