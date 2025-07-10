"use client";

import { useSelectedNode } from "@/hooks/use-selected-node";
import { useRightPanelState } from "@/hooks/use-right-panel";
import { Button } from "@/components/ui/button";
import { ScrollArea } from "@/components/ui/scroll-area";
import { Separator } from "@/components/ui/separator";
import { Badge } from "@/components/ui/badge";
import { cn } from "@/lib/utils";
import { IconX, IconCopy, IconDownload } from "@tabler/icons-react";
import type { MoleculeNodeData, ClusterNodeData } from "@/features/canvas/thought-tree/data/mock-data";

export function NodeDetailsPanel() {
  const { selectedNode } = useSelectedNode();
  const { isOpen, setIsOpen } = useRightPanelState();

  if (!isOpen || !selectedNode) return null;

  const isMoleculeNode = selectedNode.type === "molecule-node";
  const isClusterNode = selectedNode.type === "cluster-node";

  return (
    <div className="flex h-full flex-col">
      {/* Header */}
      <div className="flex items-center justify-between border-b border-border/50 pl-4 pr-3 py-3">
        <h2 className="text-sm font-semibold">
          {isMoleculeNode ? "Molecule Details" : isClusterNode ? "Cluster Details" : "Details"}
        </h2>
        <Button
          variant="ghost"
          size="sm"
          onClick={() => {
            setIsOpen(false)
            const { clearSelection } = useSelectedNode.getState()
            clearSelection()
          }}
          className="h-6 w-6 p-0 ml-1"
        >
          <IconX size={16} strokeWidth={2} />
        </Button>
      </div>

      {/* Content */}
      <ScrollArea className="flex-1">
        <div className="p-4 space-y-6">
          {isMoleculeNode && <MoleculeDetails data={selectedNode.data as MoleculeNodeData} />}
          {isClusterNode && <ClusterDetails data={selectedNode.data as ClusterNodeData} />}
        </div>
      </ScrollArea>

      {/* Actions */}
      <div className="p-4 border-t border-border/50 flex gap-2">
        <Button variant="outline" size="sm" className="flex-1">
          <IconCopy className="h-4 w-4 mr-1" />
          Copy
        </Button>
        <Button variant="outline" size="sm" className="flex-1">
          <IconDownload className="h-4 w-4 mr-1" />
          Export
        </Button>
      </div>
    </div>
  );
}

function MoleculeDetails({ data }: { data: MoleculeNodeData }) {
  const getPropertyColor = (property: string, value: number) => {
    switch (property) {
      case "qed":
        return value >= 0.7 ? "text-green-600" : value >= 0.5 ? "text-yellow-600" : "text-red-600";
      case "drd2":
        return value >= 0.4 ? "text-green-600" : value >= 0.25 ? "text-yellow-600" : "text-red-600";
      case "bbbp":
      case "hia":
        return value >= 0.7 ? "text-green-600" : value >= 0.5 ? "text-yellow-600" : "text-red-600";
      case "mutagenicity":
        return value <= 0.3 ? "text-green-600" : value <= 0.7 ? "text-yellow-600" : "text-red-600";
      case "sas":
        return value <= 3 ? "text-green-600" : value <= 4 ? "text-yellow-600" : "text-red-600";
      default:
        return "text-gray-600";
    }
  };

  const getPropertyLabel = (property: string) => {
    return property.toUpperCase();
  };

  return (
    <>
      {/* Basic Info */}
      <div>
        <h3 className="font-medium mb-2">Basic Information</h3>
        <div className="space-y-2">
          <div>
            <span className="text-sm text-muted-foreground">Name:</span>
            <p className="font-medium">{data.label}</p>
          </div>
          <div>
            <span className="text-sm text-muted-foreground">Type:</span>
            <Badge variant="secondary" className="ml-2">
              {data.type}
            </Badge>
          </div>
        </div>
      </div>

      <Separator />

      {/* SMILES */}
      <div>
        <h3 className="font-medium mb-2">Structure (SMILES)</h3>
        <div className="bg-muted p-3 rounded-md font-mono text-sm break-all">
          {data.smiles}
        </div>
      </div>

      <Separator />

      {/* Properties */}
      <div>
        <h3 className="font-medium mb-3">Molecular Properties</h3>
        <div className="space-y-3">
          {Object.entries(data.properties).map(([key, value]) => (
            <div key={key} className="space-y-1">
              <div className="flex justify-between items-start">
                <span className="text-sm text-muted-foreground">
                  {getPropertyLabel(key)}
                </span>
                <span className={cn("font-medium text-sm", getPropertyColor(key, value))}>
                  {typeof value === "number" ? value.toFixed(3) : value}
                </span>
              </div>
              {/* Property bar visualization */}
              <div className="w-full h-2 bg-gray-200 rounded-full overflow-hidden">
                <div
                  className={cn("h-full transition-all", 
                    getPropertyColor(key, value).replace("text-", "bg-")
                  )}
                  style={{ 
                    width: `${Math.min(100, 
                      key === "plogp" ? (value / 5) * 100 : 
                      key === "sas" ? ((6 - value) / 6) * 100 : 
                      value * 100
                    )}%` 
                  }}
                />
              </div>
            </div>
          ))}
        </div>
      </div>
    </>
  );
}

function ClusterDetails({ data }: { data: ClusterNodeData }) {
  return (
    <>
      {/* Basic Info */}
      <div>
        <h3 className="font-medium mb-2">Cluster Information</h3>
        <div className="space-y-2">
          <div>
            <span className="text-sm text-muted-foreground">Name:</span>
            <p className="font-medium">{data.label}</p>
          </div>
          <div>
            <span className="text-sm text-muted-foreground">Molecules:</span>
            <p className="font-medium">{data.moleculeCount} molecules</p>
          </div>
        </div>
      </div>

      <Separator />

      {/* Average Properties */}
      <div>
        <h3 className="font-medium mb-3">Average Properties</h3>
        <div className="space-y-3">
          {data.avgProperties && Object.entries(data.avgProperties).map(([key, value]) => (
            <div key={key} className="flex justify-between">
              <span className="text-sm text-muted-foreground capitalize">
                {key.toUpperCase()}:
              </span>
              <span className="font-medium">
                {value.toFixed(3)}
              </span>
            </div>
          ))}
        </div>
      </div>

      <Separator />

      {/* Generated Molecules */}
      <div>
        <h3 className="font-medium mb-2">Generated Molecules</h3>
        <p className="text-sm text-muted-foreground">
          This cluster contains {data.moleculeCount} generated molecules. Click on individual molecules to view their details.
        </p>
      </div>
    </>
  );
}