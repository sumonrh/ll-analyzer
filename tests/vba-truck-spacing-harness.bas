Private Sub TestMsgBox(ByVal prompt As String, Optional ByVal buttons As VbMsgBoxStyle = vbOKOnly)
    If (buttons And vbCritical) <> 0 Then Err.Raise vbObjectError + 100, "VBA test", prompt
End Sub

Private Sub TestAssert(condition As Boolean, message As String)
    If Not condition Then Err.Raise vbObjectError + 101, "VBA test", message
End Sub

Private Sub TestClose(actual As Double, expected As Double, message As String)
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
    For index = 0 To 23
        gap = 6.6 + 0.5 * index
        If index = 23 Then gap = 18#
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
    m_IsBCLTruck = True
    TestAssert TruckSpacingCount() = 24, "Expected 24 spacing configurations"
    For i = 0 To 22
        TestClose TruckSpacingAt(i), 6.6 + i * 0.5, "Spacing " & i
    Next i
    TestClose TruckSpacingAt(23), 18#, "Inclusive final spacing"
    m_IsBCLTruck = False
    TestAssert TruckSpacingCount() = 1, "CL must have one geometry"
    TestSpacingEnvelope "Truck", True, 0.75, 1#, 0#
    TestSpacingEnvelope "Truck", False, 1#, 1.25, 0#
    TestSpacingEnvelope "Lane", False, 1#, 0.8, 9#
    TestSpacingEnvelope "Lane", False, 1#, 0.8, 0#

    CreateInputSheet
    Set ws = ThisWorkbook.Worksheets("LL Input")
    TestAssert ws.Range("B8").Value = "CL-625", "New sheets default to CL"
    ws.Range("B8").ClearContents
    TestAssert Not ReadTruckModel(ws), "Blank truck model means CL"
    EnsureInputLayout ws
    TestAssert ws.Range("B8").Value = "CL-625", "Existing sheets gain CL default"
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
    ws.Range("B8").Value = "BCL-625"
    ws.Range("B10:C29").Value = "ignored BCL axle table"
    ws.Range("B32").Value = 2
    ws.Range("B33").Value = 2#
    ws.Range("B7").Value = 3
    MainAnalysis
    TestAssert ws.Range("B10").Value = "ignored BCL axle table", "BCL must preserve CL table"
    Set results = ThisWorkbook.Worksheets("LL Results")
    TestAssert InStr(CStr(results.Range("A1").Value), "BCL-625") > 0, "Results identify truck"
    Set heading = results.Cells.Find("Governing axle 3-4 spacing (m)", LookIn:=xlValues, LookAt:=xlWhole)
    TestAssert Not heading Is Nothing, "Governing spacing column missing"
    For i = 1 To 4
        TestAssert heading.Offset(i, 0).Value >= 6.6 And heading.Offset(i, 0).Value <= 18#, _
            "Diagnostic governing gap outside inclusive range"
    Next i
    ws.Range("B8").Value = "CL-625"
    ws.Range("B10:C29").ClearContents
    ws.Range("B10:B14").Value = Application.Transpose(Array(50, 125, 125, 175, 150))
    ws.Range("C10:C14").Value = Application.Transpose(Array(3.6, 1.2, 6.6, 6.6, 0))
    MainAnalysis
    TestAssert Not m_IsBCLTruck, "Switching back must restore fixed CL behavior"
    TestAssert InStr(CStr(results.Range("A1").Value), "BCL-625") = 0, "Stale BCL result label"
    Application.StatusBar = False
End Sub
