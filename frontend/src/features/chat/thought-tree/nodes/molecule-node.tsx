import type { NodeProps } from '@xyflow/react';
import { Handle, Position } from '@xyflow/react';
import { Card, CardContent } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { cn } from '@/lib/utils';
import { useSelectedNode } from '@/hooks/use-selected-node';
import type { MoleculeNodeData } from '../data/mock-data';

export default function MoleculeNode({ data, selected, id }: NodeProps) {
  const nodeData = data as MoleculeNodeData;
  const { selectedNode } = useSelectedNode();
  const isSelected = selectedNode?.id === id;
  const getNodeStyle = () => {
    switch (nodeData.type) {
      case 'initial':
        return 'border-blue-500 bg-blue-50 dark:bg-blue-950';
      case 'generated':
        return 'border-purple-500 bg-purple-50 dark:bg-purple-950';
      case 'edited':
        return 'border-green-500 bg-green-50 dark:bg-green-950';
      case 'optimized':
        return 'border-emerald-500 bg-emerald-50 dark:bg-emerald-950';
      default:
        return 'border-gray-300';
    }
  };

  const getTypeColor = () => {
    switch (nodeData.type) {
      case 'initial':
        return 'bg-blue-500';
      case 'generated':
        return 'bg-purple-500';
      case 'edited':
        return 'bg-green-500';
      case 'optimized':
        return 'bg-emerald-500';
      default:
        return 'bg-gray-500';
    }
  };

  return (
    <>
      {/* Input handle */}
      <Handle
        type="target"
        position={Position.Top}
        className="w-3 h-3 !bg-gray-400 border-2 border-white"
      />
      
      <Card
        className={cn(
          "w-64 shadow-md transition-all duration-200 p-4",
          getNodeStyle(),
          selected && "ring-2 ring-blue-400 ring-offset-2",
          isSelected && "ring-2 ring-primary ring-offset-2 shadow-lg"
        )}
      >
        <CardContent className="p-0 flex flex-col">
          {/* Molecule structure placeholder */}
          <div className="w-full h-48 bg-gray-200 rounded-md flex items-center justify-center mb-2 overflow-hidden">
            <span className="text-gray-500 text-sm text-center px-4 break-all">
              {nodeData.smiles ? nodeData.smiles.slice(0, 60) + '...' : 'Molecule'}
            </span>
          </div>
          
          {/* Label */}
          <div className="text-center mb-2">
            <h3 className="text-sm font-semibold">{nodeData.label}</h3>
          </div>
          
          {/* Key Properties */}
          <div className="flex-1 space-y-1">
            <div className="flex justify-between text-xs">
              <span className="text-gray-600">QED:</span>
              <span className={cn(
                "font-medium",
                nodeData.properties?.qed >= 0.7 ? "text-green-600" : 
                nodeData.properties?.qed >= 0.5 ? "text-yellow-600" : "text-red-600"
              )}>
                {nodeData.properties?.qed.toFixed(2)}
              </span>
            </div>
            <div className="flex justify-between text-xs">
              <span className="text-gray-600">DRD2:</span>
              <span className={cn(
                "font-medium",
                nodeData.properties?.drd2 >= 0.4 ? "text-green-600" : 
                nodeData.properties?.drd2 >= 0.25 ? "text-yellow-600" : "text-red-600"
              )}>
                {nodeData.properties?.drd2.toFixed(2)}
              </span>
            </div>
            <div className="flex justify-between text-xs">
              <span className="text-gray-600">pLogP:</span>
              <span className="font-medium">{nodeData.properties?.plogp.toFixed(2)}</span>
            </div>
          </div>
          
          {/* Type Badge */}
          <Badge 
            variant="secondary" 
            className={cn("text-white text-xs mt-2 self-center", getTypeColor())}
          >
            {nodeData.type}
          </Badge>
        </CardContent>
      </Card>

      {/* Output handle */}
      <Handle
        type="source"
        position={Position.Bottom}
        className="w-3 h-3 !bg-gray-400 border-2 border-white"
      />
    </>
  );
}