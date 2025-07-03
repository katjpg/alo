import type { Node, Edge } from '@xyflow/react';

export interface MoleculeNodeData extends Record<string, unknown> {
  label: string;
  type: 'initial' | 'generated' | 'edited';
}

export interface ClusterNodeData extends Record<string, unknown> {
  label: string;
  clusterId: string;
  moleculeCount: number;
}

export type ThoughtTreeNode = 
  | Node<MoleculeNodeData, 'molecule-node'>
  | Node<ClusterNodeData, 'cluster-node'>;

export type ThoughtTreeEdge = Edge<any, 'decomposition-edge' | 'generation-edge' | 'editing-edge'>;

// Mock nodes representing the ligand optimization workflow
export const initialNodes: ThoughtTreeNode[] = [
  // Initial molecule
  {
    id: 'initial-1',
    type: 'molecule-node',
    position: { x: 800, y: 100 },
    data: {
      label: 'Starting Ligand',
      type: 'initial'
    }
  },

  // Cluster nodes (G-1 through G-4) - extremely spaced out
  {
    id: 'cluster-g1',
    type: 'cluster-node',
    position: { x: 100, y: 500 },
    data: {
      label: 'G-1',
      clusterId: 'G-1',
      moleculeCount: 52
    }
  },
  {
    id: 'cluster-g2',
    type: 'cluster-node',
    position: { x: 600, y: 500 },
    data: {
      label: 'G-2',
      clusterId: 'G-2',
      moleculeCount: 48
    }
  },
  {
    id: 'cluster-g3',
    type: 'cluster-node',
    position: { x: 1100, y: 500 },
    data: {
      label: 'G-3',
      clusterId: 'G-3',
      moleculeCount: 55
    }
  },
  {
    id: 'cluster-g4',
    type: 'cluster-node',
    position: { x: 1600, y: 500 },
    data: {
      label: 'G-4',
      clusterId: 'G-4',
      moleculeCount: 45
    }
  },

  // Generated molecules from G-1 cluster - extremely spaced out
  {
    id: 'gen-g1-1',
    type: 'molecule-node',
    position: { x: 0, y: 900 },
    data: {
      label: 'G-1-1',
      type: 'generated'
    }
  },
  {
    id: 'gen-g1-2',
    type: 'molecule-node',
    position: { x: 250, y: 900 },
    data: {
      label: 'G-1-2',
      type: 'generated'
    }
  },

  // Generated molecules from G-2 cluster - extremely spaced out
  {
    id: 'gen-g2-1',
    type: 'molecule-node',
    position: { x: 500, y: 900 },
    data: {
      label: 'G-2-1',
      type: 'generated'
    }
  },
  {
    id: 'gen-g2-2',
    type: 'molecule-node',
    position: { x: 750, y: 900 },
    data: {
      label: 'G-2-2',
      type: 'generated'
    }
  },

  // Edited molecules - extremely spaced out
  {
    id: 'edit-g1-x',
    type: 'molecule-node',
    position: { x: 125, y: 1300 },
    data: {
      label: 'G-1-X',
      type: 'edited'
    }
  },
  {
    id: 'edit-g2-x',
    type: 'molecule-node',
    position: { x: 625, y: 1300 },
    data: {
      label: 'G-2-X',
      type: 'edited'
    }
  }
];

// Mock edges representing the workflow connections
export const initialEdges: ThoughtTreeEdge[] = [
  // Initial to clusters (decomposition)
  {
    id: 'initial-to-g1',
    source: 'initial-1',
    target: 'cluster-g1',
    type: 'decomposition-edge',
    animated: true
  },
  {
    id: 'initial-to-g2',
    source: 'initial-1',
    target: 'cluster-g2',
    type: 'decomposition-edge',
    animated: true
  },
  {
    id: 'initial-to-g3',
    source: 'initial-1',
    target: 'cluster-g3',
    type: 'decomposition-edge',
    animated: true
  },
  {
    id: 'initial-to-g4',
    source: 'initial-1',
    target: 'cluster-g4',
    type: 'decomposition-edge',
    animated: true
  },

  // Clusters to generated molecules (generation)
  {
    id: 'g1-to-gen1',
    source: 'cluster-g1',
    target: 'gen-g1-1',
    type: 'generation-edge'
  },
  {
    id: 'g1-to-gen2',
    source: 'cluster-g1',
    target: 'gen-g1-2',
    type: 'generation-edge'
  },
  {
    id: 'g2-to-gen1',
    source: 'cluster-g2',
    target: 'gen-g2-1',
    type: 'generation-edge'
  },
  {
    id: 'g2-to-gen2',
    source: 'cluster-g2',
    target: 'gen-g2-2',
    type: 'generation-edge'
  },

  // Generated to edited molecules (editing)
  {
    id: 'gen-g1-1-to-edit',
    source: 'gen-g1-1',
    target: 'edit-g1-x',
    type: 'editing-edge',
    style: { stroke: '#10b981', strokeWidth: 2 }
  },
  {
    id: 'gen-g2-1-to-edit',
    source: 'gen-g2-1',
    target: 'edit-g2-x',
    type: 'editing-edge',
    style: { stroke: '#10b981', strokeWidth: 2 }
  }
];

// Node types configuration
export const nodeTypes = {
  'molecule-node': 'molecule-node',
  'cluster-node': 'cluster-node'
};

// Edge types configuration
export const edgeTypes = {
  'decomposition-edge': 'decomposition-edge',
  'generation-edge': 'generation-edge',
  'editing-edge': 'editing-edge'
};