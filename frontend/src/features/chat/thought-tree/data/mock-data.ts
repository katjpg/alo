import type { Node, Edge } from '@xyflow/react';

// Extended interfaces to include molecular properties
export interface MoleculeNodeData extends Record<string, unknown> {
  label: string;
  type: 'initial' | 'generated' | 'edited' | 'optimized';
  smiles: string;
  properties: {
    bbbp: number;
    mutagenicity: number;
    hia: number;
    plogp: number;
    qed: number;
    drd2: number;
    sas: number;
  };
}

export interface ClusterNodeData extends Record<string, unknown> {
  label: string;
  clusterId: string;
  moleculeCount: number;
  avgProperties?: {
    qed: number;
    drd2: number;
    plogp: number;
  };
}

export type ThoughtTreeNode = 
  | Node<MoleculeNodeData, 'molecule-node'>
  | Node<ClusterNodeData, 'cluster-node'>;

export type ThoughtTreeEdge = Edge<any, 'decomposition-edge' | 'generation-edge' | 'editing-edge'>;

// Your actual molecular data
const molecularData = [
  {smiles: "c1ccc2c(c1)ncc(C(=O)O)c2Cl", bbbp: 0.42, mutagenicity: 0.65, hia: 0.71, plogp: 1.85, qed: 0.48, drd2: 0.18, sas: 2.45, node_id: "starting", node_type: "root"},
  {smiles: "c1ccc2c(c1)ncc(C(=O)N)c2Cl", bbbp: 0.68, mutagenicity: 0.52, hia: 0.78, plogp: 2.12, qed: 0.61, drd2: 0.22, sas: 2.38, node_id: "G-1", node_type: "cluster"},
  {smiles: "c1ccc2c(c1)ncc(C(=O)NC)c2Cl", bbbp: 0.71, mutagenicity: 0.48, hia: 0.82, plogp: 2.35, qed: 0.65, drd2: 0.25, sas: 2.41, node_id: "G-1-1", node_type: "generated"},
  {smiles: "c1ccc2c(c1)ncc(C(=O)NCC)c2Cl", bbbp: 0.74, mutagenicity: 0.45, hia: 0.85, plogp: 2.58, qed: 0.67, drd2: 0.28, sas: 2.44, node_id: "G-1-2", node_type: "generated"},
  {smiles: "c1ccc2c(c1)ncc(C(=O)NCCO)c2Cl", bbbp: 0.72, mutagenicity: 0.38, hia: 0.88, plogp: 2.21, qed: 0.71, drd2: 0.31, sas: 2.52, node_id: "G-1-X", node_type: "optimized"},
  {smiles: "c1ccc2c(c1)ncc(C(=O)O)c2F", bbbp: 0.58, mutagenicity: 0.41, hia: 0.79, plogp: 1.92, qed: 0.72, drd2: 0.35, sas: 2.31, node_id: "G-2", node_type: "cluster"},
  {smiles: "c1ccc2c(c1)ncc(C(=O)O)c2OC", bbbp: 0.61, mutagenicity: 0.38, hia: 0.83, plogp: 2.08, qed: 0.74, drd2: 0.38, sas: 2.48, node_id: "G-2-1", node_type: "generated"},
  {smiles: "c1ccc2c(c1)ncc(C(=O)N)c2OC", bbbp: 0.65, mutagenicity: 0.35, hia: 0.86, plogp: 2.15, qed: 0.76, drd2: 0.42, sas: 2.51, node_id: "G-2-2", node_type: "generated"},
  {smiles: "c1ccc2c(c1)ncc(C(=O)NCCO)c2OC", bbbp: 0.73, mutagenicity: 0.28, hia: 0.91, plogp: 2.28, qed: 0.81, drd2: 0.48, sas: 2.68, node_id: "G-2-X", node_type: "optimized"},
  {smiles: "c1ccc2c(c1)nc(C)c(C(=O)O)c2Cl", bbbp: 0.55, mutagenicity: 0.42, hia: 0.81, plogp: 2.31, qed: 0.69, drd2: 0.31, sas: 2.72, node_id: "G-3", node_type: "cluster"},
  {smiles: "c1ccc2c(c1)nc(CC)c(C(=O)O)c2Cl", bbbp: 0.59, mutagenicity: 0.39, hia: 0.84, plogp: 2.54, qed: 0.71, drd2: 0.34, sas: 2.78, node_id: "G-3-1", node_type: "generated"},
  {smiles: "c1ccc2c(c1)nc(C(C)C)c(C(=O)O)c2Cl", bbbp: 0.62, mutagenicity: 0.36, hia: 0.87, plogp: 2.78, qed: 0.73, drd2: 0.37, sas: 2.85, node_id: "G-3-2", node_type: "generated"},
  {smiles: "Cc1cc2ccccc2nc1C(=O)O", bbbp: 0.48, mutagenicity: 0.71, hia: 0.75, plogp: 1.68, qed: 0.52, drd2: 0.24, sas: 2.21, node_id: "G-4", node_type: "cluster"},
  {smiles: "Cc1cc2ccccc2nc1C(=O)N", bbbp: 0.52, mutagenicity: 0.68, hia: 0.78, plogp: 1.75, qed: 0.55, drd2: 0.27, sas: 2.24, node_id: "G-4-1", node_type: "generated"},
  {smiles: "Cc1cc2ccccc2nc1C(=O)NC", bbbp: 0.56, mutagenicity: 0.64, hia: 0.82, plogp: 1.98, qed: 0.58, drd2: 0.31, sas: 2.28, node_id: "G-4-2", node_type: "generated"},
];

// Helper function to calculate positions
function calculateNodePosition(nodeId: string, index: number): { x: number; y: number } {
  // Starting molecule at the top center
  if (nodeId === 'starting') {
    return { x: 800, y: 100 };
  }
  
  // Clusters in a row
  if (nodeId.match(/^G-\d+$/)) {
    const clusterNum = parseInt(nodeId.split('-')[1]);
    return { x: 200 + (clusterNum - 1) * 400, y: 400 };
  }
  
  // Generated molecules
  if (nodeId.match(/^G-\d+-\d+$/)) {
    const [, cluster, num] = nodeId.split('-');
    const clusterNum = parseInt(cluster);
    const molNum = parseInt(num);
    return { 
      x: 200 + (clusterNum - 1) * 400 + (molNum - 1) * 150 - 75, 
      y: 700 
    };
  }
  
  // Optimized molecules
  if (nodeId.match(/^G-\d+-X$/)) {
    const clusterNum = parseInt(nodeId.split('-')[1]);
    return { x: 200 + (clusterNum - 1) * 400, y: 1000 };
  }
  
  return { x: 0, y: 0 };
}

// Convert molecular data to nodes
export const initialNodes: ThoughtTreeNode[] = molecularData.map((mol, index) => {
  const position = calculateNodePosition(mol.node_id, index);
  
  if (mol.node_type === 'cluster') {
    // For clusters, calculate average properties from their generated molecules
    const clusterMolecules = molecularData.filter(m => 
      m.node_id.startsWith(mol.node_id + '-') && m.node_type === 'generated'
    );
    
    const avgQed = clusterMolecules.reduce((sum, m) => sum + m.qed, 0) / clusterMolecules.length || mol.qed;
    const avgDrd2 = clusterMolecules.reduce((sum, m) => sum + m.drd2, 0) / clusterMolecules.length || mol.drd2;
    const avgPlogp = clusterMolecules.reduce((sum, m) => sum + m.plogp, 0) / clusterMolecules.length || mol.plogp;
    
    return {
      id: mol.node_id,
      type: 'cluster-node',
      position,
      data: {
        label: mol.node_id,
        clusterId: mol.node_id,
        moleculeCount: clusterMolecules.length,
        avgProperties: {
          qed: parseFloat(avgQed.toFixed(2)),
          drd2: parseFloat(avgDrd2.toFixed(2)),
          plogp: parseFloat(avgPlogp.toFixed(2))
        }
      }
    };
  }
  
  // Molecule nodes
  return {
    id: mol.node_id,
    type: 'molecule-node',
    position,
    data: {
      label: mol.node_id === 'starting' ? 'Starting Ligand' : mol.node_id,
      type: mol.node_type === 'root' ? 'initial' : 
            mol.node_type === 'optimized' ? 'optimized' : 
            mol.node_type as 'generated',
      smiles: mol.smiles,
      properties: {
        bbbp: mol.bbbp,
        mutagenicity: mol.mutagenicity,
        hia: mol.hia,
        plogp: mol.plogp,
        qed: mol.qed,
        drd2: mol.drd2,
        sas: mol.sas
      }
    }
  };
});

// Generate edges based on node relationships
export const initialEdges: ThoughtTreeEdge[] = [];

// Add decomposition edges (starting → clusters)
const clusters = molecularData.filter(m => m.node_type === 'cluster');
clusters.forEach(cluster => {
  initialEdges.push({
    id: `starting-to-${cluster.node_id}`,
    source: 'starting',
    target: cluster.node_id,
    type: 'decomposition-edge',
    animated: true
  });
});

// Add generation edges (clusters → generated molecules)
molecularData.forEach(mol => {
  if (mol.node_type === 'generated') {
    const clusterId = mol.node_id.split('-').slice(0, 2).join('-');
    initialEdges.push({
      id: `${clusterId}-to-${mol.node_id}`,
      source: clusterId,
      target: mol.node_id,
      type: 'generation-edge'
    });
  }
});

// Add editing edges (generated → optimized)
molecularData.forEach(mol => {
  if (mol.node_type === 'optimized') {
    const clusterId = mol.node_id.replace('-X', '');
    // Find the best performing generated molecule to connect to optimized
    const generatedMols = molecularData.filter(m => 
      m.node_id.startsWith(clusterId + '-') && 
      m.node_type === 'generated'
    );
    
    if (generatedMols.length > 0) {
      // Connect from the first generated molecule (you could choose based on properties)
      const sourceMol = generatedMols[0];
      initialEdges.push({
        id: `${sourceMol.node_id}-to-${mol.node_id}`,
        source: sourceMol.node_id,
        target: mol.node_id,
        type: 'editing-edge',
        style: { stroke: '#10b981', strokeWidth: 2 }
      });
    }
  }
});

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

// Export the molecular data for use in other components
export { molecularData };