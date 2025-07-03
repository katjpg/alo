import type { NodeProps } from '@xyflow/react';
import { Handle, Position } from '@xyflow/react';
import { Card, CardContent } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { cn } from '@/lib/utils';
import { IconAtom } from '@tabler/icons-react';
import type { ClusterNodeData } from '../data/mock-data';

export default function ClusterNode({ data, selected }: NodeProps) {
  const nodeData = data as ClusterNodeData;
  return (
    <>
      {/* Input handle */}
      <Handle
        type="target"
        position={Position.Top}
        className="w-3 h-3 !bg-orange-400 border-2 border-white"
      />
      
      <Card
        className={cn(
          "w-36 h-40 shadow-md transition-all duration-200 p-3",
          "border-orange-500 bg-orange-50 dark:bg-orange-950",
          selected && "ring-2 ring-orange-400 ring-offset-2"
        )}
      >
        <CardContent className="p-0 flex flex-col items-center justify-between h-full">
          {/* Cluster visualization placeholder - takes majority of space */}
          <div className="w-full flex-1 bg-orange-100 dark:bg-orange-900 rounded-md flex items-center justify-center border-2 border-dashed border-orange-300 dark:border-orange-700 mb-2">
            <div className="text-center">
              <IconAtom className="h-8 w-8 text-orange-600 dark:text-orange-400 mx-auto mb-1" />
              <div className="text-orange-600 dark:text-orange-400 text-xs font-medium">
                {nodeData.moleculeCount}
              </div>
            </div>
          </div>
          
          {/* Label - centered */}
          <div className="text-center mb-1">
            <h3 className="text-sm font-semibold">{nodeData.label}</h3>
          </div>
          
          {/* Badge - centered */}
          <Badge variant="secondary" className="bg-orange-500 text-white text-xs">
            Cluster
          </Badge>
        </CardContent>
      </Card>

      {/* Output handle */}
      <Handle
        type="source"
        position={Position.Bottom}
        className="w-3 h-3 !bg-orange-400 border-2 border-white"
      />
    </>
  );
}