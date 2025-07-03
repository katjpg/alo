import type { NodeProps } from '@xyflow/react';
import { Handle, Position } from '@xyflow/react';
import { Card, CardContent } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { cn } from '@/lib/utils';
import type { MoleculeNodeData } from '../data/mock-data';

export default function MoleculeNode({ data, selected }: NodeProps) {
  const nodeData = data as MoleculeNodeData;
  const getNodeStyle = () => {
    switch (nodeData.type) {
      case 'initial':
        return 'border-blue-500 bg-blue-50 dark:bg-blue-950';
      case 'generated':
        return 'border-purple-500 bg-purple-50 dark:bg-purple-950';
      case 'edited':
        return 'border-green-500 bg-green-50 dark:bg-green-950';
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
          "w-36 h-40 shadow-md transition-all duration-200 p-3",
          getNodeStyle(),
          selected && "ring-2 ring-blue-400 ring-offset-2"
        )}
      >
        <CardContent className="p-0 flex flex-col items-center justify-between h-full">
          {/* Molecule image placeholder - takes majority of space */}
          <div className="w-full flex-1 bg-gray-300 rounded-md flex items-center justify-center mb-2">
            <span className="text-gray-500 text-xs">Molecule</span>
          </div>
          
          {/* Label - centered */}
          <div className="text-center mb-1">
            <h3 className="text-sm font-semibold">{nodeData.label}</h3>
          </div>
          
          {/* Badge - centered */}
          <Badge 
            variant="secondary" 
            className={cn("text-white text-xs", getTypeColor())}
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