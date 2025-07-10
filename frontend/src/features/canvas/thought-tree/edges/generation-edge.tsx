import {
  BaseEdge,
  getBezierPath,
  type EdgeProps,
  type Edge,
} from '@xyflow/react';

type GenerationEdgeData = {};

export type GenerationEdge = Edge<GenerationEdgeData>;

export default function GenerationEdge({
  sourceX,
  sourceY,
  targetX,
  targetY,
  sourcePosition,
  targetPosition,
  style = {},
  markerEnd,
}: EdgeProps<GenerationEdge>) {
  const [edgePath] = getBezierPath({
    sourceX,
    sourceY,
    sourcePosition,
    targetX,
    targetY,
    targetPosition,
  });

  return (
    <BaseEdge 
      path={edgePath} 
      markerEnd={markerEnd} 
      style={{
        ...style,
        stroke: '#8b5cf6', // Purple color for generation
        strokeWidth: 2,
      }} 
    />
  );
}