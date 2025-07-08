import type { NodeProps } from '@xyflow/react';
import { Handle, Position } from '@xyflow/react';
import { Card, CardContent } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { cn } from '@/lib/utils';
import { IconAtom } from '@tabler/icons-react';
import { useSelectedNode } from '@/hooks/use-selected-node';
import type { ClusterNodeData } from '../data/mock-data';

export default function ClusterNode({ data, selected, id }: NodeProps) {
  const nodeData = data as ClusterNodeData;
  const { selectedNode } = useSelectedNode();
  const isSelected = selectedNode?.id === id;
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
          "w-48 h-52 shadow-md transition-all duration-200 p-3",
          "border-orange-500 bg-orange-50 dark:bg-orange-950",
          selected && "ring-2 ring-orange-400 ring-offset-2",
          isSelected && "ring-2 ring-primary ring-offset-2 shadow-lg"
        )}
      >
        <CardContent className="p-0 flex flex-col h-full">
          {/* Cluster visualization */}
          <div className="w-full h-20 bg-orange-100 dark:bg-orange-900 rounded-md flex items-center justify-center border-2 border-dashed border-orange-300 dark:border-orange-700 mb-2">
            <div className="text-center">
              <IconAtom className="h-6 w-6 text-orange-600 dark:text-orange-400 mx-auto" />
              <div className="text-orange-600 dark:text-orange-400 text-xs font-medium">
                {nodeData.moleculeCount} molecules
              </div>
            </div>
          </div>
          
          {/* Label */}
          <div className="text-center mb-2">
            <h3 className="text-sm font-semibold">{nodeData.label}</h3>
          </div>
          
          {/* Average Properties */}
          <div className="flex-1 space-y-1">
            <div className="text-xs font-medium text-gray-600 mb-1">Avg Properties:</div>
            <div className="flex justify-between text-xs">
              <span className="text-gray-600">QED:</span>
              <span className={cn(
                "font-medium",
                nodeData.avgProperties?.qed >= 0.7 ? "text-green-600" : 
                nodeData.avgProperties?.qed >= 0.5 ? "text-yellow-600" : "text-red-600"
              )}>
                {nodeData.avgProperties?.qed.toFixed(2) || 'N/A'}
              </span>
            </div>
            <div className="flex justify-between text-xs">
              <span className="text-gray-600">DRD2:</span>
              <span className={cn(
                "font-medium",
                nodeData.avgProperties?.drd2 >= 0.4 ? "text-green-600" : 
                nodeData.avgProperties?.drd2 >= 0.25 ? "text-yellow-600" : "text-red-600"
              )}>
                {nodeData.avgProperties?.drd2.toFixed(2) || 'N/A'}
              </span>
            </div>
            <div className="flex justify-between text-xs">
              <span className="text-gray-600">pLogP:</span>
              <span className="font-medium">{nodeData.avgProperties?.plogp?.toFixed(2) || 'N/A'}</span>
            </div>
          </div>
          
          {/* Badge */}
          <Badge variant="secondary" className="bg-orange-500 text-white text-xs mt-2 self-center">
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