"use client";

import { useCallback, useEffect } from "react";
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
  type NodeMouseHandler,
  type ReactFlowProps,
} from "@xyflow/react";

import "@xyflow/react/dist/style.css";
import "@/styles/canvas.css";

import { nodeTypes } from "./nodes";
import { edgeTypes } from "./edges";
import { 
  initialNodes, 
  initialEdges, 
  type ThoughtTreeNode, 
  type ThoughtTreeEdge 
} from "./data/mock-data";
import type { ToolType } from "@/features/canvas/components/thought-tree-canvas";
import { useSelectedNode } from "@/hooks/use-selected-node";
import { useRightPanelState } from "@/hooks/use-right-panel";

interface TreeHierarchyProps {
  activeTool: ToolType;
  onToolChange: (tool: ToolType) => void;
}

export default function TreeHierarchy({ activeTool, onToolChange }: TreeHierarchyProps) {
  const [nodes, setNodes, onNodesChange] = useNodesState<ThoughtTreeNode>(initialNodes);
  const [edges, setEdges, onEdgesChange] = useEdgesState<ThoughtTreeEdge>(initialEdges);
  const { setSelectedNode, clearSelection } = useSelectedNode();
  const { setIsOpen: setRightPanelOpen } = useRightPanelState();
  
  const onConnect: OnConnect = useCallback(
    (connection) => setEdges((edges) => addEdge(connection, edges)),
    [setEdges]
  );

  // Handle canvas click based on active tool
  const handlePaneClick: ReactFlowProps["onPaneClick"] = useCallback((event) => {
    if (activeTool === "select") {
      // Clear selection when clicking on canvas
      clearSelection();
      setRightPanelOpen(false);
    } else if (activeTool === "text") {
      // Get click position
      const reactFlowBounds = (event.target as HTMLElement).getBoundingClientRect();
      const position = {
        x: event.clientX - reactFlowBounds.left,
        y: event.clientY - reactFlowBounds.top,
      };
      
      // Add a text node
      const newNode: ThoughtTreeNode = {
        id: `text-${Date.now()}`,
        type: "molecule-node", // Temporarily using molecule-node, you can create a text-node type
        position,
        data: {
          label: "New Text",
          type: "edited",
          smiles: "",
          properties: {
            bbbp: 0.5,
            mutagenicity: 0.5,
            hia: 0.5,
            plogp: 2.0,
            qed: 0.5,
            drd2: 0.3,
            sas: 3.0
          },
        },
      };
      
      setNodes((nds) => [...nds, newNode]);
      // Switch back to select mode
      onToolChange('select');
    } else if (activeTool === "node") {
      // Get click position
      const reactFlowBounds = (event.target as HTMLElement).getBoundingClientRect();
      const position = {
        x: event.clientX - reactFlowBounds.left,
        y: event.clientY - reactFlowBounds.top,
      };
      
      // Add a molecule node
      const newNode: ThoughtTreeNode = {
        id: `node-${Date.now()}`,
        type: "molecule-node",
        position,
        data: {
          label: "New Molecule",
          type: "generated",
          smiles: "C1=CC=CC=C1",
          properties: {
            bbbp: 0.65,
            mutagenicity: 0.45,
            hia: 0.75,
            plogp: 2.13,
            qed: 0.68,
            drd2: 0.35,
            sas: 2.5
          },
        },
      };
      
      setNodes((nds) => [...nds, newNode]);
      // Switch back to select mode
      onToolChange('select');
    }
  }, [activeTool, setNodes, clearSelection, setRightPanelOpen, onToolChange]);

  // Handle node click based on active tool
  const handleNodeClick: NodeMouseHandler = useCallback((event, node) => {
    if (activeTool === "select") {
      // Store the selected node and open right panel
      setSelectedNode(node as ThoughtTreeNode);
      setRightPanelOpen(true);
    } else if (activeTool === "hand") {
      // Prevent selection when using hand tool
      event.stopPropagation();
    }
  }, [activeTool, setSelectedNode, setRightPanelOpen]);

  // Configure interaction based on active tool
  const panOnDrag = activeTool === "hand" ? true : false;
  const selectionOnDrag = activeTool === "select" ? true : false;
  const nodesDraggable = activeTool === "select" ? true : false;
  const nodesConnectable = activeTool === "select" ? true : false;
  const elementsSelectable = activeTool === "select" ? true : false;

  // Add keyboard shortcuts
  useEffect(() => {
    const handleKeyPress = (event: KeyboardEvent) => {
      // Let parent component handle tool switching
      // This is just for tool-specific actions
      if (event.key === "Delete" || event.key === "Backspace") {
        if (activeTool === "select") {
          // Delete selected nodes
          setNodes((nds) => nds.filter((node) => !node.selected));
          setEdges((eds) => eds.filter((edge) => !edge.selected));
        }
      }
    };

    window.addEventListener("keydown", handleKeyPress);
    return () => window.removeEventListener("keydown", handleKeyPress);
  }, [activeTool, setNodes, setEdges]);

  return (
    <div className="w-full h-full" data-tool={activeTool}>
      <ReactFlow<ThoughtTreeNode, ThoughtTreeEdge>
        nodes={nodes}
        nodeTypes={nodeTypes}
        onNodesChange={onNodesChange}
        edges={edges}
        edgeTypes={edgeTypes}
        onEdgesChange={onEdgesChange}
        onConnect={onConnect}
        onPaneClick={handlePaneClick}
        onNodeClick={handleNodeClick}
        panOnDrag={panOnDrag}
        selectionOnDrag={selectionOnDrag}
        nodesDraggable={nodesDraggable}
        nodesConnectable={nodesConnectable}
        elementsSelectable={elementsSelectable}
        fitView
        minZoom={0.1}
        maxZoom={4}
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
          orientation="horizontal"
        />
      </ReactFlow>
    </div>
  );
}