import MoleculeNode from './molecule-node';
import ClusterNode from './cluster-node';

export { MoleculeNode, ClusterNode };

export const nodeTypes = {
  'molecule-node': MoleculeNode,
  'cluster-node': ClusterNode,
} as const;

// Export types
export type { MoleculeNodeData, ClusterNodeData } from '../data/mock-data';