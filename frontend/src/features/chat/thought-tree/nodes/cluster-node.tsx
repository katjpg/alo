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
            <IconAtom className="h-5 w-5 text-gray-500 flex-shrink-0" />
            <span className="flex-1 truncate text-sm font-medium text-gray-900">
              {nodeData.label}
            </span>
            <Badge 
              variant="outline"
              className="text-xs capitalize bg-orange-100 text-orange-700 border-orange-200"
            >
              Cluster
            </Badge>
          </div>

          {/* Cluster visualization */}
          <div className="px-3 py-3 border-b border-gray-100">
            <div className="w-full h-40 bg-gray-50 rounded-md flex items-center justify-center border border-gray-200">
              <div className="text-center">
                <IconAtom className="h-12 w-12 text-gray-400 mx-auto mb-2" />
                <div className="text-gray-600 text-sm font-medium">
                  {nodeData.moleculeCount} molecules
                </div>
              </div>
            </div>
          </div>
          
          {/* Average Properties section */}
          <div className="px-3 py-3 space-y-2">
            <div className="text-xs font-medium text-gray-500 uppercase tracking-wider">Average Properties</div>
            <div className="flex justify-between items-center text-xs">
              <span className="text-gray-600">QED:</span>
              <span className={cn(
                "font-medium",
                nodeData.avgProperties?.qed >= 0.7 ? "text-green-600" : 
                nodeData.avgProperties?.qed >= 0.5 ? "text-yellow-600" : "text-red-600"
              )}>
                {nodeData.avgProperties?.qed?.toFixed(2) || 'N/A'}
              </span>
            </div>
            <div className="flex justify-between items-center text-xs">
              <span className="text-gray-600">DRD2:</span>
              <span className={cn(
                "font-medium",
                nodeData.avgProperties?.drd2 >= 0.4 ? "text-green-600" : 
                nodeData.avgProperties?.drd2 >= 0.25 ? "text-yellow-600" : "text-red-600"
              )}>
                {nodeData.avgProperties?.drd2?.toFixed(2) || 'N/A'}
              </span>
            </div>
            <div className="flex justify-between items-center text-xs">
              <span className="text-gray-600">pLogP:</span>
              <span className="font-medium text-gray-900">
                {nodeData.avgProperties?.plogp?.toFixed(2) || 'N/A'}
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