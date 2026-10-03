import { analyzeBeam } from './beam-engine';
import type { AnalysisRequest, AnalysisResponse } from './beam-engine';

const send = (response: AnalysisResponse) => self.postMessage(response);

self.onmessage = (event: MessageEvent<AnalysisRequest>) => {
    try {
        const result = analyzeBeam(event.data, progress => send({ type: 'progress', progress }));
        send({ type: 'result', result });
    } catch (error) {
        console.error('Live load analysis failed:', error);
        send({ type: 'error', message: error instanceof Error ? error.message : 'Unexpected analysis error.' });
    }
};
