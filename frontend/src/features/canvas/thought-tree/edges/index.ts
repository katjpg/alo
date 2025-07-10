import DecompositionEdge from './decomposition-edge';
import GenerationEdge from './generation-edge';
import EditingEdge from './editing-edge';

export { DecompositionEdge, GenerationEdge, EditingEdge };

export const edgeTypes = {
  'decomposition-edge': DecompositionEdge,
  'generation-edge': GenerationEdge,
  'editing-edge': EditingEdge,
};

// Export types
export type { DecompositionEdge as DecompositionEdgeType } from './decomposition-edge';
export type { GenerationEdge as GenerationEdgeType } from './generation-edge';
export type { EditingEdge as EditingEdgeType } from './editing-edge';