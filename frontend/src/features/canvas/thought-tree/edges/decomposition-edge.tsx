import {
  BaseEdge,
  getBezierPath,
  type EdgeProps,
  type Edge,
} from '@xyflow/react';

type DecompositionEdgeData = {};

export type DecompositionEdge = Edge<DecompositionEdgeData>;

export default function DecompositionEdge({
  sourceX,
  sourceY,
  targetX,
  targetY,
  sourcePosition,
  targetPosition,
  style = {},
  markerEnd,
}: EdgeProps<DecompositionEdge>) {
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
        stroke: '#3b82f6', // Blue color for decomposition
        strokeWidth: 2,
        strokeDasharray: '5,5', // Dashed line
      }} 
    />
  );
}