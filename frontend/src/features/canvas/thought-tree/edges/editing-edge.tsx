import {
  BaseEdge,
  getBezierPath,
  type EdgeProps,
  type Edge,
} from '@xyflow/react';

type EditingEdgeData = {};

export type EditingEdge = Edge<EditingEdgeData>;

export default function EditingEdge({
  sourceX,
  sourceY,
  targetX,
  targetY,
  sourcePosition,
  targetPosition,
  style = {},
  markerEnd,
}: EdgeProps<EditingEdge>) {
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
        stroke: '#10b981', // Green color for editing
        strokeWidth: 3,
        strokeLinecap: 'round',
      }} 
    />
  );
}