import { useState, useEffect, useRef, useMemo } from 'react';
import { Play, RotateCcw, Plus, Trash2, Settings, AlertCircle, Download } from 'lucide-react';
import AnalysisWorker from './analysis.worker?worker&inline';
import {
    DEFAULT_SPANS, DEFAULT_AXLES, DEFAULT_CONFIG, MAX_AXLES,
    computeAutoDlaInfo, computeEffectiveIncrement,
} from './beam-engine';
import type {
    Span, Axle, AnalysisConfig, AnalysisResults, AnalysisResponse,
    EnvelopePoint, ReactionEnvelope, LoadCase, CaseResults, AnalysisProgress,
} from './beam-engine';

type SheetJs = {
    utils: {
        book_new(): object;
        json_to_sheet(rows: object[]): object;
        book_append_sheet(workbook: object, sheet: object, name: string): void;
    };
    writeFile(workbook: object, name: string): void;
};
declare global {
    interface Window { XLSX?: SheetJs }
}

// --- COMPONENTS ---

// Helper to calculate nice ticks for chart axes
const calculateTicks = (min: number, max: number, targetCount: number) => {
    if (min === max) return [min];
    const span = max - min;
    const step = Math.pow(10, Math.floor(Math.log10(span / targetCount)));
    const err = targetCount / (span / step);

    let finalStep = step;
    if (err <= .15) finalStep *= 10;
    else if (err <= .35) finalStep *= 5;
    else if (err <= .75) finalStep *= 2;

    const start = Math.ceil(min / finalStep) * finalStep;
    const end = Math.floor(max / finalStep) * finalStep;

    const ticks = [];
    const decimals = Math.max(0, -Math.floor(Math.log10(finalStep)));

    for (let val = start; val <= end + (finalStep / 2); val += finalStep) {
        const cleanVal = parseFloat(val.toFixed(decimals));
        if (cleanVal >= min && cleanVal <= max) ticks.push(cleanVal);
    }
    return ticks;
};

const EnvelopeChart = ({
    data,
    dataKeyMax,
    dataKeyMin,
    title,
    unit,
    color,
    flipY = false,
    xAxisTitle = 'Span position (m)'
}: {
    data: EnvelopePoint[],
    dataKeyMax: 'max',
    dataKeyMin: 'min',
    title: string,
    unit: string,
    color: string,
    flipY?: boolean,
    xAxisTitle?: string
}) => {
    const containerRef = useRef<HTMLDivElement>(null);
    const [width, setWidth] = useState(600);
    const height = 300;
    const padding = { top: 40, right: 30, bottom: 50, left: 70 };
    const [hoverData, setHoverData] = useState<{
        x: number;
        svgX: number;
        max: number;
        min: number;
    } | null>(null);

    useEffect(() => {
        const container = containerRef.current;
        if (!container) return;
        const observer = new ResizeObserver(entries => {
            const measured = entries[0].contentRect.width;
            if (measured > 0) setWidth(measured);
        });
        observer.observe(container);
        return () => observer.disconnect();
    }, []);

    if (!data || data.length === 0) return <div className="h-[300px] flex items-center justify-center text-gray-400" > No Data </div>;

    const xVals = data.map(d => d.x);
    const maxVals = data.map(d => d[dataKeyMax]);
    const minVals = data.map(d => d[dataKeyMin]);
    const allY = [...maxVals, ...minVals];

    const xMin = Math.min(...xVals);
    const xMax = Math.max(...xVals);
    let yMin = Math.min(...allY);
    let yMax = Math.max(...allY);

    // Find global extrema
    let globalMaxVal = -Infinity;
    let globalMaxX = xMin;
    let globalMinVal = Infinity;
    let globalMinX = xMin;

    for (let i = 0; i < data.length; i++) {
        const pt = data[i];
        const curMax = pt[dataKeyMax];
        const curMin = pt[dataKeyMin];
        if (curMax > globalMaxVal) {
            globalMaxVal = curMax;
            globalMaxX = pt.x;
        }
        if (curMin < globalMinVal) {
            globalMinVal = curMin;
            globalMinX = pt.x;
        }
    }

    const yRange = yMax - yMin;
    if (yRange === 0) { yMax += 1; yMin -= 1; }
    else { yMax += yRange * 0.1; yMin -= yRange * 0.1; }

    const xTicks = calculateTicks(xMin, xMax, 8);
    const yTicks = calculateTicks(yMin, yMax, 6);

    const xScale = (val: number) => padding.left + ((val - xMin) / (xMax - xMin)) * (width - padding.left - padding.right);
    const yScale = (val: number) =>
        flipY
            ? padding.top + ((val - yMin) / (yMax - yMin)) * (height - padding.top - padding.bottom)
            : height - padding.bottom - ((val - yMin) / (yMax - yMin)) * (height - padding.top - padding.bottom);

    const createPath = (vals: number[]) => {
        return vals.map((y, i) => `${i === 0 ? 'M' : 'L'} ${xScale(data[i].x)} ${yScale(y)}`).join(' ');
    };

    const pathMax = createPath(maxVals);
    const pathMin = createPath(minVals);

    const pathFill = `${pathMax} L ${xScale(data[data.length - 1].x)} ${yScale(minVals[minVals.length - 1])} ` +
        minVals.slice().reverse().map((y, i) => `L ${xScale(data[data.length - 1 - i].x)} ${yScale(y)}`).join(' ') + " Z";

    const zeroY = yScale(0);

    const unitSymbol = unit.includes('(') ? unit.split('(')[1].replace(')', '') : unit;

    const formatVal = (val: number) => {
        if (!Number.isFinite(val)) return '0.00';
        if (Math.abs(val) < 1e-6) return '0.00';
        if (Math.abs(val) < 0.005) return val.toFixed(4);
        return val.toFixed(2);
    };

    const formatValWithDetail = (val: number) => {
        if (!Number.isFinite(val)) return '0.00 ' + unitSymbol;
        if (unitSymbol.toLowerCase() === 'm') {
            const mm = (val * 1000).toFixed(2);
            return `${val.toFixed(4)} m (${mm} mm)`;
        }
        return `${formatVal(val)} ${unitSymbol}`;
    };

    const handlePointerMove = (e: React.PointerEvent<SVGSVGElement>) => {
        const svgRect = e.currentTarget.getBoundingClientRect();
        if (svgRect.width <= 0) return;
        const clientX = e.clientX - svgRect.left;
        const plotWidth = width - padding.left - padding.right;
        if (plotWidth <= 0) return;

        const rawRatio = (clientX - padding.left) / plotWidth;
        const clampedRatio = Math.max(0, Math.min(1, rawRatio));
        const targetX = xMin + clampedRatio * (xMax - xMin);

        // Find interpolated values or nearest segment
        let bestMax = 0;
        let bestMin = 0;
        let found = false;

        for (let i = 0; i < data.length - 1; i++) {
            const x1 = data[i].x;
            const x2 = data[i + 1].x;
            if (x1 === x2) continue; // vertical step boundary in stepped shear
            const minSegX = Math.min(x1, x2);
            const maxSegX = Math.max(x1, x2);
            if (targetX >= minSegX && targetX <= maxSegX) {
                const t = (targetX - x1) / (x2 - x1);
                bestMax = data[i][dataKeyMax] + t * (data[i + 1][dataKeyMax] - data[i][dataKeyMax]);
                bestMin = data[i][dataKeyMin] + t * (data[i + 1][dataKeyMin] - data[i][dataKeyMin]);
                found = true;
                break;
            }
        }

        if (!found) {
            let nearestDist = Infinity;
            let nearestIdx = 0;
            for (let i = 0; i < data.length; i++) {
                const dist = Math.abs(data[i].x - targetX);
                if (dist < nearestDist) {
                    nearestDist = dist;
                    nearestIdx = i;
                }
            }
            bestMax = data[nearestIdx][dataKeyMax];
            bestMin = data[nearestIdx][dataKeyMin];
        }

        const currentSvgX = xScale(targetX);
        setHoverData({
            x: targetX,
            svgX: currentSvgX,
            max: bestMax,
            min: bestMin
        });
    };

    const handlePointerLeave = () => {
        setHoverData(null);
    };

    const tooltipWidth = 205;
    const tooltipHeight = 68;
    const tooltipX = hoverData
        ? hoverData.svgX > width - tooltipWidth - 20
            ? Math.max(padding.left + 5, hoverData.svgX - tooltipWidth - 12)
            : Math.min(width - padding.right - tooltipWidth - 5, hoverData.svgX + 12)
        : 0;
    const tooltipY = Math.max(
        padding.top + 6,
        Math.min(height - padding.bottom - tooltipHeight - 6, padding.top + 8)
    );

    return (
        <div ref={containerRef} className="w-full bg-white rounded-lg shadow-sm border border-gray-200 p-4 mb-6" >
            <div className="flex flex-wrap items-center justify-between gap-2 mb-2">
                <h3 className="text-lg font-semibold text-gray-800">{title}</h3>
                <div className="flex flex-wrap items-center gap-2 text-xs">
                    <div className="flex items-center gap-1.5 px-2.5 py-1 rounded-md bg-blue-50 border border-blue-200 text-blue-800">
                        <span className="font-medium text-gray-500">Max:</span>
                        <span className="font-bold font-mono">{formatVal(globalMaxVal)} {unitSymbol}</span>
                        <span className="text-gray-500">@ {globalMaxX.toFixed(2)}m</span>
                    </div>
                    <div className="flex items-center gap-1.5 px-2.5 py-1 rounded-md bg-red-50 border border-red-200 text-red-800">
                        <span className="font-medium text-gray-500">Min:</span>
                        <span className="font-bold font-mono">{formatVal(globalMinVal)} {unitSymbol}</span>
                        <span className="text-gray-500">@ {globalMinX.toFixed(2)}m</span>
                    </div>
                </div>
            </div>

            <svg
                width={width}
                height={height}
                className="overflow-visible select-none"
                onPointerMove={handlePointerMove}
                onPointerLeave={handlePointerLeave}
            >

                {/* X-Grid & Labels */}
                {
                    xTicks.map(tick => {
                        const xPos = xScale(tick);
                        return (
                            <g key={`x-${tick}`
                            }>
                                <line x1={xPos} y1={padding.top} x2={xPos} y2={height - padding.bottom} stroke="#e5e7eb" strokeWidth="1" />
                                <text x={xPos} y={height - padding.bottom + 15
                                } textAnchor="middle" fontSize="10" fill="#6b7280" > {tick} </text>
                            </g>
                        );
                    })}

                {/* Y-Grid & Labels */}
                {
                    yTicks.map(tick => {
                        const yPos = yScale(tick);
                        return (
                            <g key={`y-${tick}`
                            }>
                                <line x1={padding.left} y1={yPos} x2={width - padding.right} y2={yPos} stroke="#e5e7eb" strokeWidth="1" />
                                <text x={padding.left - 8} y={yPos + 3} textAnchor="end" fontSize="10" fill="#6b7280" > {tick} </text>
                            </g>
                        );
                    })}

                <line x1={padding.left} y1={padding.top} x2={padding.left} y2={height - padding.bottom} stroke="#374151" strokeWidth="1" />
                <line x1={padding.left} y1={height - padding.bottom} x2={width - padding.right} y2={height - padding.bottom} stroke="#374151" strokeWidth="1" />

                {zeroY > padding.top && zeroY < height - padding.bottom && (
                    <line x1={padding.left} y1={zeroY} x2={width - padding.right} y2={zeroY} stroke="#9ca3af" strokeWidth="1.5" strokeDasharray="4 4" />
                )}

                <text
                    x={padding.left + (width - padding.left - padding.right) / 2}
                    y={height - 10}
                    textAnchor="middle"
                    fontSize="12"
                    fontWeight="500"
                    fill="#374151"
                >
                    {xAxisTitle}
                </text>

                <text
                    x={15}
                    y={padding.top + (height - padding.top - padding.bottom) / 2}
                    textAnchor="middle"
                    fontSize="12"
                    fontWeight="500"
                    fill="#374151"
                    className="transform -rotate-90 origin-center"
                    style={{ transformBox: 'fill-box' }}
                >
                    {unit}
                </text>

                <path d={pathFill} fill={color} fillOpacity="0.1" />
                <path d={pathMax} fill="none" stroke={color} strokeWidth="2" />
                <path d={pathMin} fill="none" stroke="#ef4444" strokeWidth="2" strokeDasharray="4 2" />

                {/* Extrema markers on the curves */}
                {Number.isFinite(globalMaxVal) && (
                    <circle
                        cx={xScale(globalMaxX)}
                        cy={yScale(globalMaxVal)}
                        r="4.5"
                        fill={color}
                        stroke="#ffffff"
                        strokeWidth="1.5"
                    >
                        <title>{`Max: ${formatVal(globalMaxVal)} ${unitSymbol} at x = ${globalMaxX.toFixed(2)} m`}</title>
                    </circle>
                )}
                {Number.isFinite(globalMinVal) && (
                    <circle
                        cx={xScale(globalMinX)}
                        cy={yScale(globalMinVal)}
                        r="4.5"
                        fill="#ef4444"
                        stroke="#ffffff"
                        strokeWidth="1.5"
                    >
                        <title>{`Min: ${formatVal(globalMinVal)} ${unitSymbol} at x = ${globalMinX.toFixed(2)} m`}</title>
                    </circle>
                )}

                {/* Dynamic Cursor Overlay */}
                {hoverData && (
                    <g pointerEvents="none">
                        <line
                            x1={hoverData.svgX}
                            y1={padding.top}
                            x2={hoverData.svgX}
                            y2={height - padding.bottom}
                            stroke="#475569"
                            strokeWidth="1.2"
                            strokeDasharray="3 3"
                        />
                        <circle
                            cx={hoverData.svgX}
                            cy={yScale(hoverData.max)}
                            r="4.5"
                            fill={color}
                            stroke="#ffffff"
                            strokeWidth="2"
                        />
                        <circle
                            cx={hoverData.svgX}
                            cy={yScale(hoverData.min)}
                            r="4.5"
                            fill="#ef4444"
                            stroke="#ffffff"
                            strokeWidth="2"
                        />

                        {/* Tooltip Card */}
                        <rect
                            x={tooltipX}
                            y={tooltipY}
                            width={tooltipWidth}
                            height={tooltipHeight}
                            rx="6"
                            fill="#0f172a"
                            opacity="0.94"
                        />
                        <text x={tooltipX + 10} y={tooltipY + 18} fill="#f8fafc" fontSize="11" fontWeight="700">
                            x = {hoverData.x.toFixed(2)} m
                        </text>
                        <text x={tooltipX + 10} y={tooltipY + 36} fill="#60a5fa" fontSize="10" fontWeight="500">
                            Max: {formatValWithDetail(hoverData.max)}
                        </text>
                        <text x={tooltipX + 10} y={tooltipY + 54} fill="#f87171" fontSize="10" fontWeight="500">
                            Min: {formatValWithDetail(hoverData.min)}
                        </text>
                    </g>
                )}

                {/* Interactive Overlay Layer */}
                <rect
                    x={padding.left}
                    y={padding.top}
                    width={width - padding.left - padding.right}
                    height={height - padding.top - padding.bottom}
                    fill="transparent"
                    className="cursor-crosshair"
                />

            </svg>
            <div className="flex justify-center gap-6 mt-2 text-sm" >
                <div className="flex items-center" > <div className="w-4 h-0.5 bg-[color:var(--color)] mr-2" style={{ backgroundColor: color }}> </div> Max Envelope</div >
                <div className="flex items-center" > <div className="w-4 h-0.5 bg-red-500 mr-2 border-dashed border-t-2 border-red-500" > </div> Min Envelope</div >
            </div>
        </div>
    );
};

// Beam Schematic for Configuration View (shows beam with supports)
const BeamSchematic = ({ spans }: { spans: Span[] }) => {
    const containerRef = useRef<HTMLDivElement>(null);
    const [width, setWidth] = useState(600);
    const height = 120;
    const padding = { top: 30, right: 40, bottom: 40, left: 40 };

    useEffect(() => {
        if (containerRef.current) setWidth(containerRef.current.clientWidth);
        const handleResize = () => containerRef.current && setWidth(containerRef.current.clientWidth);
        window.addEventListener('resize', handleResize);
        return () => window.removeEventListener('resize', handleResize);
    }, []);

    const totalLen = spans.reduce((a, b) => a + b.length, 0);
    if (totalLen === 0) return null;

    const beamY = height / 2;
    const xScale = (val: number) => padding.left + (val / totalLen) * (width - padding.left - padding.right);

    // Calculate support positions
    const supportPositions: number[] = [0];
    let cumLen = 0;
    for (const span of spans) {
        cumLen += span.length;
        supportPositions.push(cumLen);
    }

    // Triangle support symbol
    const triangleSize = 12;
    const renderSupport = (x: number, idx: number) => {
        const xPos = xScale(x);
        return (
            <g key={`support-${idx}`}>
                {/* Triangle */}
                <polygon
                    points={`${xPos},${beamY + 4} ${xPos - triangleSize},${beamY + triangleSize + 8} ${xPos + triangleSize},${beamY + triangleSize + 8}`}
                    fill="#374151"
                    stroke="#1f2937"
                    strokeWidth="1"
                />
                {/* Ground line */}
                <line
                    x1={xPos - triangleSize - 4}
                    y1={beamY + triangleSize + 10}
                    x2={xPos + triangleSize + 4}
                    y2={beamY + triangleSize + 10}
                    stroke="#1f2937"
                    strokeWidth="2"
                />
                {/* Support label */}
                <text
                    x={xPos}
                    y={beamY + triangleSize + 25}
                    textAnchor="middle"
                    fontSize="10"
                    fill="#6b7280"
                >
                    {idx === 0 ? 'A' : String.fromCharCode(65 + idx)}
                </text>
            </g>
        );
    };

    return (
        <div ref={containerRef} className="w-full bg-white rounded-lg shadow-sm border border-gray-200 p-4 mb-4">
            <h3 className="text-sm font-semibold text-gray-700 mb-2">Beam Configuration</h3>
            <svg width={width} height={height} className="overflow-visible">
                {/* Beam line */}
                <line
                    x1={xScale(0)}
                    y1={beamY}
                    x2={xScale(totalLen)}
                    y2={beamY}
                    stroke="#2563eb"
                    strokeWidth="6"
                    strokeLinecap="round"
                />

                {/* Span labels and dimension lines */}
                {spans.map((span, idx) => {
                    let startX = 0;
                    for (let i = 0; i < idx; i++) startX += spans[i].length;
                    const endX = startX + span.length;
                    const midX = (startX + endX) / 2;

                    return (
                        <g key={`span-${idx}`}>
                            {/* Dimension line */}
                            <line
                                x1={xScale(startX) + 2}
                                y1={beamY - 20}
                                x2={xScale(endX) - 2}
                                y2={beamY - 20}
                                stroke="#9ca3af"
                                strokeWidth="1"
                                markerStart="url(#arrowLeft)"
                                markerEnd="url(#arrowRight)"
                            />
                            {/* Span length label */}
                            <text
                                x={xScale(midX)}
                                y={beamY - 25}
                                textAnchor="middle"
                                fontSize="11"
                                fill="#374151"
                                fontWeight="500"
                            >
                                {span.length}m
                            </text>
                        </g>
                    );
                })}

                {/* Arrow markers definition */}
                <defs>
                    <marker id="arrowLeft" markerWidth="6" markerHeight="6" refX="0" refY="3" orient="auto">
                        <path d="M6,0 L0,3 L6,6" fill="none" stroke="#9ca3af" strokeWidth="1" />
                    </marker>
                    <marker id="arrowRight" markerWidth="6" markerHeight="6" refX="6" refY="3" orient="auto">
                        <path d="M0,0 L6,3 L0,6" fill="none" stroke="#9ca3af" strokeWidth="1" />
                    </marker>
                </defs>

                {/* Supports */}
                {supportPositions.map((pos, idx) => renderSupport(pos, idx))}
            </svg>
        </div>
    );
};

// Beam Reaction Diagram for Results View (shows beam with max reaction values)
const BeamReactionDiagram = ({
    spans,
    reactions,
    supportPositions,
    dla,
    dlaAuto,
    dlaMultiplier,
}: {
    spans: Span[],
    reactions: ReactionEnvelope[],
    supportPositions: number[],
    dla?: number,
    dlaAuto?: boolean,
    dlaMultiplier?: number,
}) => {
    const containerRef = useRef<HTMLDivElement>(null);
    const [width, setWidth] = useState(600);
    const height = 200;
    const padding = { top: 30, right: 60, bottom: 60, left: 60 };

    useEffect(() => {
        if (containerRef.current) setWidth(containerRef.current.clientWidth);
        const handleResize = () => containerRef.current && setWidth(containerRef.current.clientWidth);
        window.addEventListener('resize', handleResize);
        return () => window.removeEventListener('resize', handleResize);
    }, []);

    const totalLen = spans.reduce((a, b) => a + b.length, 0);
    if (totalLen === 0 || !reactions || reactions.length === 0) return null;

    const beamY = 55;
    const xScale = (val: number) => padding.left + (val / totalLen) * (width - padding.left - padding.right);

    const triangleSize = 10;

    const renderSupportWithReaction = (x: number, idx: number, reaction: ReactionEnvelope) => {
        const xPos = xScale(x);
        // Reactions are upward-positive kN (sign fixed). Show governing max; min (uplift) in title.
        const displayR = reaction.max;

        return (
            <g key={`support-${idx}`}>
                {/* Triangle support */}
                <polygon
                    points={`${xPos},${beamY + 4} ${xPos - triangleSize},${beamY + triangleSize + 6} ${xPos + triangleSize},${beamY + triangleSize + 6}`}
                    fill="#374151"
                    stroke="#1f2937"
                    strokeWidth="1"
                />

                {/* Ground line under triangle */}
                <line
                    x1={xPos - triangleSize - 3}
                    y1={beamY + triangleSize + 8}
                    x2={xPos + triangleSize + 3}
                    y2={beamY + triangleSize + 8}
                    stroke="#1f2937"
                    strokeWidth="2"
                />

                {/* Support label (A, B, C, D) - positioned to the side */}
                <text
                    x={xPos}
                    y={beamY + triangleSize + 25}
                    textAnchor="middle"
                    fontSize="11"
                    fill="#374151"
                    fontWeight="600"
                >
                    {String.fromCharCode(65 + idx)}
                </text>

                {/* Reaction arrow (upward) - line and arrowhead */}
                <line
                    x1={xPos}
                    y1={beamY + triangleSize + 55}
                    x2={xPos}
                    y2={beamY + triangleSize + 35}
                    stroke="#dc2626"
                    strokeWidth="2.5"
                />
                {/* Arrowhead pointing UP */}
                <polygon
                    points={`${xPos - 5},${beamY + triangleSize + 35} ${xPos},${beamY + triangleSize + 27} ${xPos + 5},${beamY + triangleSize + 35}`}
                    fill="#dc2626"
                />

                {/* Reaction value - at bottom */}
                <text
                    x={xPos}
                    y={beamY + triangleSize + 72}
                    textAnchor="middle"
                    fontSize="11"
                    fill="#dc2626"
                    fontWeight="600"
                >
                    <title>{`Max ${reaction.max.toFixed(1)} kN, Min ${reaction.min.toFixed(1)} kN${reaction.govPos !== undefined ? `, gov at ${reaction.govPos.toFixed(2)}m` : ''}`}</title>
                    {displayR.toFixed(1)} kN
                </text>
            </g>
        );
    };

    return (
        <div ref={containerRef} className="w-full bg-white rounded-lg shadow-sm border border-gray-200 p-4 mb-6">
            <h3 className="text-lg font-semibold text-gray-800 mb-2">Maximum Support Reactions</h3>
            <svg width={width} height={height} className="overflow-visible">
                {/* Reaction arrow marker - pointing upward */}
                <defs>
                    <marker id="reactionArrow" markerWidth="10" markerHeight="10" refX="5" refY="10" orient="auto">
                        <path d="M0,10 L5,0 L10,10 Z" fill="#dc2626" />
                    </marker>
                </defs>

                {/* Beam line */}
                <line
                    x1={xScale(0)}
                    y1={beamY}
                    x2={xScale(totalLen)}
                    y2={beamY}
                    stroke="#2563eb"
                    strokeWidth="6"
                    strokeLinecap="round"
                />

                {/* Span labels - above beam */}
                {spans.map((span, idx) => {
                    let startX = 0;
                    for (let i = 0; i < idx; i++) startX += spans[i].length;
                    const endX = startX + span.length;
                    const midX = (startX + endX) / 2;

                    return (
                        <g key={`span-label-${idx}`}>
                            {/* Dimension line */}
                            <line
                                x1={xScale(startX) + 5}
                                y1={beamY - 18}
                                x2={xScale(endX) - 5}
                                y2={beamY - 18}
                                stroke="#9ca3af"
                                strokeWidth="1"
                            />
                            {/* End ticks */}
                            <line x1={xScale(startX) + 5} y1={beamY - 14} x2={xScale(startX) + 5} y2={beamY - 22} stroke="#9ca3af" strokeWidth="1" />
                            <line x1={xScale(endX) - 5} y1={beamY - 14} x2={xScale(endX) - 5} y2={beamY - 22} stroke="#9ca3af" strokeWidth="1" />
                            {/* Label */}
                            <text
                                x={xScale(midX)}
                                y={beamY - 25}
                                textAnchor="middle"
                                fontSize="10"
                                fill="#6b7280"
                            >
                                {span.length}m
                            </text>
                        </g>
                    );
                })}

                {/* Supports with reactions */}
                {supportPositions.map((pos, idx) =>
                    reactions[idx] && renderSupportWithReaction(pos, idx, reactions[idx])
                )}
            </svg>
            <div className="flex justify-center gap-4 mt-1 text-xs text-gray-500">
                <span>↑ Maximum reaction forces shown{dlaAuto ? ` (Auto truck DLA 40%/30%/25% × d=${((dlaMultiplier ?? 1) * 100).toFixed(0)}%; support maxima are continuous optima)` : dla !== undefined ? ` (Applied DLA = ${(dla * 100).toFixed(0)}%)` : ''}</span>
            </div>
        </div>
    );
};

export default function BeamAnalysisApp() {
    const [spans, setSpans] = useState<Span[]>(DEFAULT_SPANS);
    const [axles, setAxles] = useState<Axle[]>(DEFAULT_AXLES);
    const [config, setConfig] = useState<AnalysisConfig>(DEFAULT_CONFIG);
    const [results, setResults] = useState<AnalysisResults | null>(null);
    const [isAnalyzing, setIsAnalyzing] = useState(false);
    const [activeTab, setActiveTab] = useState<'config' | 'results'>('config');
    const [resultCase, setResultCase] = useState<LoadCase>('truck');
    const [progress, setProgress] = useState<AnalysisProgress | null>(null);
    const [analysisError, setAnalysisError] = useState<string | null>(null);
    const workerRef = useRef<Worker | null>(null);
    const displayed: CaseResults | null = results ? results.cases[resultCase] ?? results : null;

    useEffect(() => () => workerRef.current?.terminate(), []);

    const autoDlaInfo = useMemo(() => computeAutoDlaInfo(spans, axles), [spans, axles]);
    const stepInfo = useMemo(
        () => computeEffectiveIncrement(spans, axles, config.truckIncrement, config.nElemsPerSpan),
        [spans, axles, config.truckIncrement, config.nElemsPerSpan]
    );

    // Inject SheetJS script for Excel Export
    useEffect(() => {
        const script = document.createElement('script');
        script.src = "https://cdn.sheetjs.com/xlsx-0.20.1/package/dist/xlsx.full.min.js";
        script.async = true;
        script.onerror = () => {
            console.error('Could not load the Excel export library.');
            setAnalysisError('Excel export is unavailable. Check your internet connection and reload the app.');
        };
        document.body.appendChild(script);
        return () => {
            document.body.removeChild(script);
        }
    }, []);

    const addSpan = () => {
        const newId = `s${Date.now()}`;
        setSpans([...spans, { id: newId, length: 20 }]);
    };

    const removeSpan = (id: string) => {
        if (spans.length <= 1) return;
        setSpans(spans.filter(s => s.id !== id));
    };

    const updateSpan = (id: string, val: number) => {
        setSpans(spans.map(s => s.id === id ? { ...s, length: val } : s));
    };

    const runAnalysis = () => {
        workerRef.current?.terminate();
        setAnalysisError(null);
        setIsAnalyzing(true);
        setProgress({ fraction: 0, message: 'Starting analysis...' });
        const fail = (message: string) => {
            console.error('Analysis failed:', message);
            setAnalysisError(message);
            setIsAnalyzing(false);
            setProgress(null);
            workerRef.current?.terminate();
            workerRef.current = null;
        };
        try {
            const worker = new AnalysisWorker();
            workerRef.current = worker;
            worker.onmessage = (event: MessageEvent<AnalysisResponse>) => {
                const response = event.data;
                if (response.type === 'progress') {
                    setProgress(response.progress);
                    return;
                }
                if (response.type === 'error') {
                    fail(response.message);
                    return;
                }
                setResults(response.result);
                setResultCase(response.result.loadCase);
                setActiveTab('results');
                setIsAnalyzing(false);
                setProgress(null);
                worker.terminate();
                workerRef.current = null;
            };
            worker.onerror = event => {
                event.preventDefault();
                fail(event.message || 'Could not run the analysis worker.');
            };
            worker.onmessageerror = () => fail('Could not read the analysis worker results.');
            worker.postMessage({ spans, axles, config });
        } catch (error) {
            fail(error instanceof Error ? error.message : 'Could not start the analysis.');
        }
    };

    const downloadExcel = () => {
        if (!results) return;

        if (typeof window === 'undefined' || !window.XLSX) {
            alert("Excel export library is loading. Please try again in a few seconds.");
            return;
        }

        const xlsx = window.XLSX;
        const wb = xlsx.utils.book_new();

        const formatData = (data: EnvelopePoint[]) => data.map(d => ({
            "Position (m)": d.x,
            "Max": d.max,
            "Min": d.min
        }));

        const append = (name: string, rows: object[]) =>
            xlsx.utils.book_append_sheet(wb, xlsx.utils.json_to_sheet(rows), name);
        append('Analysis Settings', [{
            'Load Case': results.loadCase,
            'Base Increment (m)': results.baseIncrement,
            'Effective Increment (m)': results.incrementUsed,
            'Step Control': results.incrementReason,
            'Elements per Span': results.config.nElemsPerSpan,
            'Elastic Modulus (Pa)': results.config.E,
            'Moment of Inertia (m^4)': results.config.I,
            'DLA Multiplier d': results.config.dlaMultiplier ?? 1,
            'Lane UDL (kN/m)': results.config.laneUdl ?? 9,
            'Method': 'FEA + influence-line UDL zones + continuous truck optimisation (VBA parity)',
            'Elapsed (ms)': results.elapsedMs,
        }]);
        append('Spans', results.spans.map((span, i) => ({ 'Span': i + 1, 'Length (m)': span.length })));
        append('Axles', results.axles.map((axle, i) => ({
            'Axle': i + 1, 'Load (kN)': axle.load,
            'Spacing to Next (m)': i < results.axles.length - 1 ? axle.spacing : 0,
        })));
        for (const name of ['truck', 'lane', 'envelope'] as const) {
            const data = results.cases[name];
            if (!data) continue;
            append(`${name} Shear`, formatData(data.shear));
            append(`${name} Moment`, formatData(data.moment));
            append(`${name} Deflection`, formatData(data.deflection));
            append(`${name} Reactions`, data.reactions.map((r, i) => ({
                'Support': `Support ${i + 1}`, 'Location (m)': r.x,
                'Max Reaction (kN)': r.max, 'Min Reaction (kN)': r.min,
                'Gov Truck Pos (m)': r.govPos, 'Applied DLA': data.dlaUsed,
                'DLA Auto': data.dlaAuto, 'DLA Base': data.dlaBase, 'DLA Multiplier d': data.dlaMultiplier,
            })));
            append(`${name} Reaction Diagrams`, results.truckPositions.map((x, p) => {
                const row: Record<string, number> = { 'Truck Position (m)': x };
                data.reactionDiagrams.forEach((diagram, s) => {
                    row[`Support ${s + 1} Max (kN)`] = diagram[p].max;
                    row[`Support ${s + 1} Min (kN)`] = diagram[p].min;
                });
                return row;
            }));
        }
        try {
            xlsx.writeFile(wb, 'beam_analysis_results.xlsx');
        } catch (error) {
            console.error('Excel export failed:', error);
            setAnalysisError(`Excel export failed: ${error instanceof Error ? error.message : 'Unexpected export error.'}`);
        }
    };

    return (
        <div className="min-h-screen bg-gray-50 text-slate-800 font-sans" >
            {/* Header */}
            < header className="bg-blue-700 text-white p-4 shadow-md sticky top-0 z-10" >
                <div className="max-w-5xl mx-auto flex justify-between items-center" >
                    <h1 className="text-xl font-bold flex items-center gap-2" >
                        <span className="bg-white text-blue-700 p-1 rounded font-black text-xs" > FEM </span>
                        Beam Analysis < span className="text-blue-200 font-normal text-sm hidden sm:inline" >| CL - 625 Truck Moving Load </span>
                    </h1>
                    < button
                        onClick={runAnalysis}
                        disabled={isAnalyzing}
                        className={`flex items-center gap-2 px-4 py-2 rounded font-medium transition-colors ${isAnalyzing ? 'bg-blue-800 cursor-wait' : 'bg-white text-blue-700 hover:bg-blue-50'}`
                        }
                    >
                        {isAnalyzing ? <RotateCcw className="animate-spin w-4 h-4" /> : <Play className="w-4 h-4" />}
                        {isAnalyzing ? 'Calculating...' : 'Run Analysis'}
                    </button>
                </div>
            </header>

            < main className="max-w-5xl mx-auto p-4 md:p-6" >
                {analysisError && (
                    <div role="alert" className="bg-red-50 border border-red-200 text-red-800 p-3 rounded mb-4">
                        {analysisError}
                    </div>
                )}
                {progress && (
                    <div role="status" className="bg-blue-50 border border-blue-200 text-blue-800 p-3 rounded mb-4">
                        {progress.message}
                        <progress className="w-full mt-2" value={progress.fraction} max={1} />
                    </div>
                )}

                {/* Tabs */}
                < div className="flex gap-4 border-b border-gray-200 mb-6" >
                    <button
                        onClick={() => setActiveTab('config')}
                        className={`pb-2 px-1 font-medium text-sm transition-colors ${activeTab === 'config' ? 'text-blue-600 border-b-2 border-blue-600' : 'text-gray-500 hover:text-gray-700'}`}
                    >
                        Configuration
                    </button>
                    < button
                        onClick={() => setActiveTab('results')}
                        disabled={!results}
                        className={`pb-2 px-1 font-medium text-sm transition-colors ${activeTab === 'results' ? 'text-blue-600 border-b-2 border-blue-600' : 'text-gray-500 hover:text-gray-700 disabled:opacity-50'}`}
                    >
                        Results
                    </button>
                </div>

                {
                    activeTab === 'config' && (
                        <>
                            <BeamSchematic spans={spans} />
                            <div className="grid grid-cols-1 md:grid-cols-2 gap-6" >
                                {/* Spans Card */}
                                < div className="bg-white p-6 rounded-lg shadow-sm border border-gray-200" >
                                    <div className="flex justify-between items-center mb-4" >
                                        <h2 className="text-lg font-semibold flex items-center gap-2" >
                                            <Settings className="w-5 h-5 text-gray-500" /> Geometry
                                        </h2>
                                        < button onClick={addSpan} className="text-sm bg-blue-50 text-blue-600 px-3 py-1 rounded hover:bg-blue-100 flex items-center gap-1" >
                                            <Plus className="w-3 h-3" /> Add Span
                                        </button>
                                    </div>

                                    < div className="space-y-3" >
                                        {
                                            spans.map((span, idx) => (
                                                <div key={span.id} className="flex items-center gap-3 p-3 bg-gray-50 rounded border border-gray-100" >
                                                    <span className="text-sm font-bold text-gray-400 w-8" >#{idx + 1} </span>
                                                    < div className="flex-1" >
                                                        <label className="text-xs text-gray-500 block" > Length(m) </label>
                                                        < input
                                                            type="number"
                                                            value={span.length}
                                                            onChange={(e) => updateSpan(span.id, parseFloat(e.target.value) || 0)
                                                            }
                                                            className="w-full bg-white border border-gray-300 rounded px-2 py-1 text-sm focus:ring-2 focus:ring-blue-500 outline-none"
                                                        />
                                                    </div>
                                                    < button onClick={() => removeSpan(span.id)} className="text-gray-400 hover:text-red-500 p-2" >
                                                        <Trash2 className="w-4 h-4" />
                                                    </button>
                                                </div>
                                            ))}
                                    </div>
                                    < div className="mt-4 pt-4 border-t border-gray-100 text-sm text-gray-500 flex justify-between" >
                                        <span>Total Length: </span>
                                        < span className="font-mono font-bold text-gray-800" > {spans.reduce((a, b) => a + b.length, 0).toFixed(2)} m </span>
                                    </div>
                                </div>

                                {/* Config Card */}
                                <div className="space-y-6" >
                                    <div className="bg-white p-6 rounded-lg shadow-sm border border-gray-200" >
                                        <h2 className="text-lg font-semibold mb-4" > Analysis Settings </h2>
                                        < div className="space-y-4" >
                                            <div>
                                                <label className="block text-sm font-medium text-gray-700 mb-1" > Load Case </label>
                                                < select
                                                    value={config.loadCase}
                                                    onChange={(e) => setConfig({ ...config, loadCase: e.target.value as 'truck' | 'lane' | 'envelope' })}
                                                    className="w-full bg-white border border-gray-300 rounded px-3 py-2 text-sm focus:ring-2 focus:ring-blue-500 outline-none"
                                                >
                                                    <option value="truck" > CL-625 Truck Only (Standard) </option>
                                                    <option value="lane" >{`CL-625 Lane Load (80% Truck, no DLA + ${(config.laneUdl ?? 9)} kN/m patterned)`}</option>
                                                    <option value="envelope" > Envelope (max of Truck and Lane) </option>
                                                </select>
                                            </div>

                                            <div>
                                                <label className="block text-sm font-medium text-gray-700 mb-1" > Lane UDL (kN/m) </label>
                                                <div className="flex items-center gap-2">
                                                    <input
                                                        type="number"
                                                        step="0.5"
                                                        min="0"
                                                        value={config.laneUdl ?? 9}
                                                        onChange={(e) => {
                                                            const val = parseFloat(e.target.value);
                                                            setConfig({
                                                                ...config,
                                                                laneUdl: isNaN(val) ? 9 : Math.max(0, val),
                                                            });
                                                        }}
                                                        className="w-1/3 bg-white border border-gray-300 rounded px-2.5 py-1 text-sm focus:ring-2 focus:ring-blue-500 outline-none font-semibold"
                                                        placeholder="e.g. 9"
                                                    />
                                                    <span className="text-xs text-gray-600">
                                                        Default 9 kN/m (CL-625); use 7 or 8 for evaluation. Applies to Lane and Envelope cases.
                                                    </span>
                                                </div>
                                            </div>

                                            <div>
                                                <div className="flex items-center justify-between mb-1.5">
                                                    <label className="block text-sm font-medium text-gray-700">
                                                        Dynamic Load Allowance (DLA)
                                                    </label>
                                                    <label className="flex items-center gap-1.5 text-xs text-gray-600 cursor-pointer select-none">
                                                        <input
                                                            type="checkbox"
                                                            checked={config.dlaOverride !== null && config.dlaOverride !== undefined}
                                                            onChange={(e) => {
                                                                const checked = e.target.checked;
                                                                setConfig({
                                                                    ...config,
                                                                    dlaOverride: checked ? (config.dlaOverride ?? autoDlaInfo.dla) : null,
                                                                });
                                                            }}
                                                            className="rounded border-gray-300 text-blue-600 focus:ring-blue-500 h-3.5 w-3.5 cursor-pointer"
                                                        />
                                                        <span className="font-medium">Override DLA</span>
                                                    </label>
                                                </div>

                                                {config.dlaOverride === null || config.dlaOverride === undefined ? (
                                                    <div className="bg-slate-50 border border-slate-200 rounded p-2.5 flex items-center justify-between">
                                                        <div className="flex items-center gap-2">
                                                            <span className="inline-flex items-center px-2 py-0.5 rounded text-xs font-bold bg-blue-100 text-blue-800">
                                                                Auto: 40%/30%/25% × d={(config.dlaMultiplier ?? 1).toFixed(2)}
                                                            </span>
                                                            <span className="text-xs text-slate-600">
                                                                Span-based uniform ({autoDlaInfo.desc})
                                                            </span>
                                                        </div>
                                                        <span className="text-[10px] text-slate-400 font-medium">CSA S6 Cl. 3.8.4.5</span>
                                                    </div>
                                                ) : (
                                                    <div className="space-y-1.5">
                                                        <div className="flex items-center gap-2">
                                                            <input
                                                                type="number"
                                                                step="0.01"
                                                                min="0"
                                                                max="1.0"
                                                                value={config.dlaOverride}
                                                                onChange={(e) => {
                                                                    const val = parseFloat(e.target.value);
                                                                    setConfig({
                                                                        ...config,
                                                                        dlaOverride: isNaN(val) ? 0 : val,
                                                                    });
                                                                }}
                                                                className="w-1/3 bg-white border border-blue-400 rounded px-2.5 py-1 text-sm focus:ring-2 focus:ring-blue-500 outline-none font-semibold text-blue-900"
                                                                placeholder="e.g. 0.30"
                                                            />
                                                            <span className="text-xs text-amber-800 bg-amber-50 border border-amber-200 rounded px-2.5 py-1">
                                                                Manual override active: <strong>{((config.dlaOverride ?? 0) * 100).toFixed(1)}%</strong>
                                                            </span>
                                                        </div>
                                                    </div>
                                                )}
                                                <span className="text-[11px] text-gray-500 mt-1 block">
                                                    Automated per CSA S6 Cl. 3.8.4.5 with FEA + influence-line placement (verified vs Midas Civil).
                                                    Uniform span-based DLA: 40% (1 axle), 30% (2 axles / tandem), 25% (≥3 axles), × d multiplier.
                                                    Lane UDL (9 kN/m) uses exact positive/negative influence zones incl. partial elements, no DLA; lane truck uses 80% with no DLA.
                                                </span>
                                                <div className="mt-2">
                                                    <label className="block text-sm font-medium text-gray-700 mb-1">
                                                        Truck DLA Multiplier, d (0–1)
                                                    </label>
                                                    <div className="flex items-center gap-2">
                                                        <input
                                                            type="number"
                                                            step="0.01"
                                                            min="0"
                                                            max="1"
                                                            value={config.dlaMultiplier ?? 1}
                                                            onChange={(e) => {
                                                                const val = parseFloat(e.target.value);
                                                                setConfig({
                                                                    ...config,
                                                                    dlaMultiplier: isNaN(val) ? 1 : Math.max(0, Math.min(1, val)),
                                                                });
                                                            }}
                                                            className="w-1/3 bg-white border border-gray-300 rounded px-2.5 py-1 text-sm focus:ring-2 focus:ring-blue-500 outline-none font-semibold"
                                                            placeholder="e.g. 1.00"
                                                        />
                                                        <span className="text-xs text-gray-600">
                                                            d=0 off, d=1 full, d=0.75 = 25% less. Applies to truck DLA (auto or override).
                                                        </span>
                                                    </div>
                                                </div>
                                            </div>

                                            < div className="grid grid-cols-2 gap-4" >
                                                <div>
                                                    <label className="block text-sm font-medium text-gray-700 mb-1" > Elastic Modulus(Pa) </label>
                                                    < input
                                                        type="number"
                                                        value={config.E}
                                                        onChange={(e) => setConfig({ ...config, E: parseFloat(e.target.value) })}
                                                        className="w-full border border-gray-300 rounded px-2 py-1 text-sm"
                                                    />
                                                </div>
                                                < div >
                                                    <label className="block text-sm font-medium text-gray-700 mb-1" > Inertia(m⁴) </label>
                                                    < input
                                                        type="number"
                                                        value={config.I}
                                                        onChange={(e) => setConfig({ ...config, I: parseFloat(e.target.value) })}
                                                        className="w-full border border-gray-300 rounded px-2 py-1 text-sm"
                                                    />
                                                </div>
                                            </div>

                                            <div className="grid grid-cols-2 gap-4" >
                                                <div>
                                                    <label className="block text-sm font-medium text-gray-700 mb-1" > Elements per span </label>
                                                    < input
                                                        type="number"
                                                        min="2"
                                                        max="200"
                                                        step="1"
                                                        value={config.nElemsPerSpan}
                                                        onChange={(e) => setConfig({ ...config, nElemsPerSpan: Number(e.target.value) })}
                                                        className="w-full border border-gray-300 rounded px-2 py-1 text-sm"
                                                    />
                                                </div>
                                                <div>
                                                    <label className="block text-sm font-medium text-gray-700 mb-1" > Truck step, base (m) </label>
                                                    < input
                                                        type="number"
                                                        min="0.02"
                                                        max="2"
                                                        step="0.05"
                                                        value={config.truckIncrement}
                                                        onChange={(e) => setConfig({ ...config, truckIncrement: Number(e.target.value) })}
                                                        className="w-full border border-gray-300 rounded px-2 py-1 text-sm"
                                                    />
                                                </div>
                                            </div>

                                            < div className="bg-blue-50 p-3 rounded text-sm text-blue-800 flex gap-2 items-start" >
                                                <AlertCircle className="w-4 h-4 mt-0.5 shrink-0" />
                                                <p>Mesh: {config.nElemsPerSpan} elements per span. Truck sweep uses <strong>{stepInfo.effective.toFixed(3)}m</strong> steps{stepInfo.wasAdjusted ? <> (adjusted from {config.truckIncrement}m — {stepInfo.reason})</> : <> (base {config.truckIncrement}m)</>}. Exact support alignments are included in both directions.</p>
                                            </div>
                                        </div>
                                    </div>

                                    < div className="bg-white p-6 rounded-lg shadow-sm border border-gray-200" >
                                        <h2 className="text-lg font-semibold mb-2" > Truck Configuration </h2>
                                        <p className="text-sm text-gray-500 mb-4">CL-625 defaults; customize up to {MAX_AXLES} axles. Spacing is to the next axle.</p>
                                        <div className="space-y-2">
                                            {
                                                axles.map((axle, i) => (
                                                    <div key={axle.id} className="flex items-end gap-2 bg-gray-50 rounded p-2">
                                                        <label className="flex-1 text-xs text-gray-600">
                                                            Axle {i + 1} load (kN)
                                                            <input type="number" min="0" value={axle.load}
                                                                onChange={e => setAxles(axles.map(a => a.id === axle.id ? { ...a, load: Number(e.target.value) } : a))}
                                                                className="block w-full border border-gray-300 rounded p-1 mt-1" />
                                                        </label>
                                                        {i < axles.length - 1 && (
                                                            <label className="flex-1 text-xs text-gray-600">
                                                                Spacing to next (m)
                                                                <input type="number" min="0" step="0.1" value={axle.spacing}
                                                                    onChange={e => setAxles(axles.map(a => a.id === axle.id ? { ...a, spacing: Number(e.target.value) } : a))}
                                                                    className="block w-full border border-gray-300 rounded p-1 mt-1" />
                                                            </label>
                                                        )}
                                                        <button aria-label={`Remove axle ${i + 1}`} disabled={axles.length <= 1}
                                                            onClick={() => setAxles(axles.filter(a => a.id !== axle.id))}
                                                            className="p-2 text-gray-400 hover:text-red-500 disabled:opacity-30">
                                                            <Trash2 className="w-4 h-4" />
                                                        </button>
                                                    </div>
                                                ))}
                                        </div>
                                        <div className="flex gap-3 mt-3">
                                            <button disabled={axles.length >= MAX_AXLES}
                                                onClick={() => setAxles([...axles.map((a, i) => i === axles.length - 1 ? { ...a, spacing: 3.6 } : a),
                                                    { id: `a${Date.now()}`, load: 100, spacing: 0 }])}
                                                className="text-sm text-blue-600 disabled:opacity-30">Add Axle</button>
                                            <button onClick={() => setAxles(DEFAULT_AXLES)} className="text-sm text-gray-600">Reset CL-625</button>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </>
                    )}

                {
                    activeTab === 'results' && results && displayed && (
                        <div className="animate-in fade-in slide-in-from-bottom-4 duration-500" >
                            <div className="flex flex-wrap justify-between items-center gap-3 mb-4">
                                <label className="text-sm font-medium">
                                    Display load case
                                    <select aria-label="Display load case" value={resultCase}
                                        onChange={e => setResultCase(e.target.value as LoadCase)}
                                        className="ml-2 border border-gray-300 rounded p-2">
                                        {(['truck', 'lane', 'envelope'] as const).filter(name => results.cases[name]).map(name =>
                                            <option key={name} value={name}>{name === 'truck' ? 'Truck' : name === 'lane' ? 'Lane' : 'Combined Envelope'}</option>)}
                                    </select>
                                </label>
                                <span className="text-xs text-gray-500">
                                    {results.elapsedMs.toFixed(0)}ms | {results.stats.truckSolves} truck positions | {results.stats.udlSolves} span UDL solves
                                </span>
                            </div>
                            {results.incrementUsed !== undefined && (
                                <div className="w-full bg-blue-50 border border-blue-200 rounded-lg px-4 py-2 mb-4 text-xs text-blue-800">
                                    Sweep step used: <strong>{results.incrementUsed.toFixed(3)}m</strong>
                                    {results.baseIncrement !== undefined && Math.abs(results.incrementUsed - results.baseIncrement) > 1e-9
                                        ? <> (adjusted from base {results.baseIncrement.toFixed(3)}m: {results.incrementReason})</>
                                        : <> (base setting)</>}.
                                    Exact axle/support alignment positions included. Envelopes use continuous truck optimisation; UDL uses exact influence zones.
                                    {displayed.dlaAuto
                                        ? <> Truck DLA auto (span-based uniform): 40%/30%/25% × d={(displayed.dlaMultiplier ?? 1).toFixed(2)}.</>
                                        : <> Truck DLA effective: {((displayed.dlaUsed ?? 0) * 100).toFixed(2)}%{resultCase !== 'lane' ? ` (base ${((displayed.dlaBase ?? 0) * 100).toFixed(2)}% × d=${(displayed.dlaMultiplier ?? 1).toFixed(2)})` : ' (lane: no DLA)' }.</>}
                                </div>
                            )}
                            <BeamReactionDiagram
                                spans={results.spans}
                                reactions={displayed.reactions}
                                supportPositions={results.supportPositions}
                                dla={displayed.dlaUsed}
                                dlaAuto={displayed.dlaAuto}
                                dlaMultiplier={displayed.dlaMultiplier}
                            />
                            <EnvelopeChart
                                title="Shear Force Envelope"
                                data={displayed.shear}
                                dataKeyMax="max"
                                dataKeyMin="min"
                                unit="Shear (kN)"
                                color="#2563eb"
                            />

                            <EnvelopeChart
                                title="Bending Moment Envelope"
                                data={displayed.moment}
                                dataKeyMax="max"
                                dataKeyMin="min"
                                unit="Moment (kNm)"
                                color="#059669"
                                flipY={true}
                            />

                            <EnvelopeChart
                                title="Deflection Envelope"
                                data={displayed.deflection}
                                dataKeyMax="max"
                                dataKeyMin="min"
                                unit="Deflection (m)"
                                color="#9333ea"
                            />
                            <div className="bg-white rounded-lg border border-gray-200 p-4 mb-6 overflow-x-auto">
                                <h3 className="text-lg font-semibold mb-3">Support Reaction Summary</h3>
                                <table className="w-full text-sm text-right">
                                    <thead><tr className="border-b">
                                        <th className="text-left">Support</th><th>Location (m)</th>
                                        <th>Max (kN)</th><th>Min / uplift (kN)</th><th>Governing truck position (m)</th>
                                    </tr></thead>
                                    <tbody>{displayed.reactions.map((r, s) => (
                                        <tr key={s} className="border-b border-gray-100">
                                            <td className="text-left py-2">Support {s + 1}</td>
                                            <td>{r.x.toFixed(2)}</td><td>{r.max.toFixed(2)}</td>
                                            <td>{r.min.toFixed(2)}</td><td>{r.govPos.toFixed(3)}</td>
                                        </tr>
                                    ))}</tbody>
                                </table>
                            </div>

                            <div className="bg-white p-6 rounded-lg shadow-sm border border-gray-200 mt-6" >
                                <h3 className="font-semibold mb-4" > Export Data </h3>
                                < div className="text-sm text-gray-600 mb-4" >
                                    Download full-precision Excel results: settings, spans, axles, shear, moment, deflection, support summaries and reaction diagrams. Envelope mode includes Truck, Lane and Combined Envelope sheets.
                                </div>
                                < button
                                    onClick={downloadExcel}
                                    className="bg-green-700 text-white px-4 py-2 rounded text-sm hover:bg-green-800 flex items-center gap-2"
                                >
                                    <Download className="w-4 h-4" />
                                    Download Excel(XLSX)
                                </button>
                            </div>
                        </div>
                    )
                }
            </main>
        </div>
    );
}