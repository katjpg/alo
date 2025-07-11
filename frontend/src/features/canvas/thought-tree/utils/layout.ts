import dagre from '@dagrejs/dagre';
import type { Node, Edge } from '@xyflow/react';
import type { ThoughtTreeNode, ThoughtTreeEdge } from '../data/mock-data';

const dagreGraph = new dagre.graphlib.Graph().setDefaultEdgeLabel(() => ({}));

// Node dimensions
const nodeWidth = 256; // Based on w-64 (16rem = 256px)
const nodeHeight = 320; // Approximate height of molecule nodes

export const getLayoutedElements = (
  nodes: ThoughtTreeNode[],
  edges: ThoughtTreeEdge[],
  direction = 'LR'
) => {
  const isHorizontal = direction === 'LR';
  
  // Configure dagre graph
  dagreGraph.setGraph({ 
    rankdir: direction,
    nodesep: 100, // Space between nodes in same rank
    ranksep: 200, // Space between ranks
    marginx: 50,
    marginy: 50,
  });

  // Add nodes to dagre
  nodes.forEach((node) => {
    dagreGraph.setNode(node.id, { width: nodeWidth, height: nodeHeight });
  });

  // Add edges to dagre
  edges.forEach((edge) => {
    dagreGraph.setEdge(edge.source, edge.target);
  });

  // Run the layout algorithm
  dagre.layout(dagreGraph);

  // Map the results back to React Flow nodes
  const layoutedNodes = nodes.map((node) => {
    const nodeWithPosition = dagreGraph.node(node.id);
    
    // Center the node position (dagre gives top-left)
    const x = nodeWithPosition.x - nodeWidth / 2;
    const y = nodeWithPosition.y - nodeHeight / 2;

    return {
      ...node,
      position: { x, y },
      // Update handle positions for horizontal flow
      targetPosition: isHorizontal ? 'left' : 'top',
      sourcePosition: isHorizontal ? 'right' : 'bottom',
    } as ThoughtTreeNode;
  });

  return { nodes: layoutedNodes, edges };
};

// Helper function to relayout when nodes/edges change
export const applyDagreLayout = (
  nodes: ThoughtTreeNode[],
  edges: ThoughtTreeEdge[],
  direction = 'LR'
): { nodes: ThoughtTreeNode[]; edges: ThoughtTreeEdge[] } => {
  // Clear the graph for fresh layout
  dagreGraph.nodes().forEach((n) => dagreGraph.removeNode(n));
  dagreGraph.edges().forEach((e) => dagreGraph.removeEdge(e));
  
  return getLayoutedElements(nodes, edges, direction);
};