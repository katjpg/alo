import type { NodeProps } from '@xyflow/react';
import { Handle, Position } from '@xyflow/react';
import { Card, CardContent } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { cn } from '@/lib/utils';
import { useSelectedNode } from '@/hooks/use-selected-node';
import { IconFlask } from '@tabler/icons-react';
import type { MoleculeNodeData } from '../data/mock-data';

export default function MoleculeNode({ data, selected, id }: NodeProps) {
  const nodeData = data as MoleculeNodeData;
  const { selectedNode } = useSelectedNode();
  const isSelected = selectedNode?.id === id;

  const getTypeBadgeColor = () => {
    switch (nodeData.type) {
      case 'initial':
        return 'bg-blue-100 text-blue-700 border-blue-200';
      case 'generated':
        return 'bg-purple-100 text-purple-700 border-purple-200';
      case 'edited':
        return 'bg-green-100 text-green-700 border-green-200';
      case 'optimized':
        return 'bg-emerald-100 text-emerald-700 border-emerald-200';
      default:
        return 'bg-gray-100 text-gray-700 border-gray-200';
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
          "w-64 bg-white border-gray-200 shadow-sm hover:shadow-md transition-all duration-200",
          selected && "ring-2 ring-gray-400 ring-offset-2",
          isSelected && "ring-2 ring-primary ring-offset-2 shadow-lg"
        )}
      >
        <CardContent className="p-0">
          {/* Header section */}
          <div className="flex items-center gap-2 px-3 py-2 border-b border-gray-100">
            <IconFlask className="h-5 w-5 text-gray-500 flex-shrink-0" />
            <span className="flex-1 truncate text-sm font-medium text-gray-900">
              {nodeData.label}
            </span>
            <Badge 
              variant="outline"
              className={cn("text-xs capitalize", getTypeBadgeColor())}
            >
              {nodeData.type}
            </Badge>
          </div>

          {/* Molecule structure */}
          <div className="px-3 py-3 border-b border-gray-100">
            <div className="w-full h-40 bg-gray-50 rounded-md flex items-center justify-center overflow-hidden">
              <span className="text-gray-400 text-xs text-center px-4 font-mono break-all">
                {nodeData.smiles ? nodeData.smiles.slice(0, 50) + '...' : 'Molecule Structure'}
              </span>
            </div>
          </div>
          
          {/* Properties section */}
          <div className="px-3 py-3 space-y-2">
            <div className="flex justify-between items-center text-xs">
              <span className="text-gray-600">QED:</span>
              <span className={cn(
                "font-medium",
                nodeData.properties?.qed >= 0.7 ? "text-green-600" : 
                nodeData.properties?.qed >= 0.5 ? "text-yellow-600" : "text-red-600"
              )}>
                {nodeData.properties?.qed.toFixed(2)}
              </span>
            </div>
            <div className="flex justify-between items-center text-xs">
              <span className="text-gray-600">DRD2:</span>
              <span className={cn(
                "font-medium",
                nodeData.properties?.drd2 >= 0.4 ? "text-green-600" : 
                nodeData.properties?.drd2 >= 0.25 ? "text-yellow-600" : "text-red-600"
              )}>
                {nodeData.properties?.drd2.toFixed(2)}
              </span>
            </div>
            <div className="flex justify-between items-center text-xs">
              <span className="text-gray-600">pLogP:</span>
              <span className="font-medium text-gray-900">
                {nodeData.properties?.plogp.toFixed(2)}
              </span>
            </div>
          </div>
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