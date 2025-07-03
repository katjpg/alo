"use client";

import { useCallback } from "react";
import {
  Background,
  Controls,
  MiniMap,
  ReactFlow,
  addEdge,
  useNodesState,
  useEdgesState,
  type OnConnect,
  BackgroundVariant,
} from "@xyflow/react";

import "@xyflow/react/dist/style.css";

import { nodeTypes } from "./nodes";
import { edgeTypes } from "./edges";
import { initialNodes, initialEdges, type ThoughtTreeNode, type ThoughtTreeEdge } from "./data/mock-data";

export default function TreeHierarchy() {
  const [nodes, , onNodesChange] = useNodesState<ThoughtTreeNode>(initialNodes);
  const [edges, setEdges, onEdgesChange] = useEdgesState<ThoughtTreeEdge>(initialEdges);
  
  const onConnect: OnConnect = useCallback(
    (connection) => setEdges((edges) => addEdge(connection, edges)),
    [setEdges]
  );

  return (
    <div className="w-full h-full">
      <ReactFlow<ThoughtTreeNode, ThoughtTreeEdge>
        nodes={nodes}
        nodeTypes={nodeTypes}
        onNodesChange={onNodesChange}
        edges={edges}
        edgeTypes={edgeTypes}
        onEdgesChange={onEdgesChange}
        onConnect={onConnect}
        fitView
        proOptions={{ hideAttribution: true }}
      >
        <Background 
          variant={BackgroundVariant.Dots} 
          gap={20} 
          size={1}
        />
        <MiniMap 
          className="!bg-white border border-gray-200"
          nodeStrokeWidth={3}
          nodeColor={(node) => {
            switch (node.type) {
              case 'molecule-node':
                return '#8b5cf6';
              case 'cluster-node':
                return '#f97316';
              default:
                return '#6b7280';
            }
          }}
        />
        <Controls 
          className="!bg-white border border-gray-200"
          showInteractive={false}
        />
      </ReactFlow>
    </div>
  );
}