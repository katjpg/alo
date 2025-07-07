"use client";

import { Button } from "@/components/ui/button";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import { 
  IconPointer, 
  IconHandMove, 
  IconLetterT,
  IconHexagon
} from "@tabler/icons-react";
import { cn } from "@/lib/utils";
import type { ToolType } from "./thought-tree-canvas";

interface CanvasToolbarProps {
  activeTool: ToolType;
  onToolChange: (tool: ToolType) => void;
}

const tools = [
  {
    id: "select" as ToolType,
    icon: IconPointer,
    label: "Select",
    shortcut: "V"
  },
  {
    id: "hand" as ToolType,
    icon: IconHandMove,
    label: "Hand (Pan)",
    shortcut: "H"
  },
  {
    id: "text" as ToolType,
    icon: IconLetterT,
    label: "Text",
    shortcut: "T"
  },
  {
    id: "node" as ToolType,
    icon: IconHexagon,
    label: "Add Node",
    shortcut: "N"
  }
];

export function CanvasToolbar({ activeTool, onToolChange }: CanvasToolbarProps) {
  return (
    <TooltipProvider delayDuration={0}>
      <div className="absolute bottom-4 left-1/2 -translate-x-1/2 z-10">
        <div className="flex items-center gap-1 rounded-full border bg-background/95 p-1 shadow-lg backdrop-blur-sm">
          {tools.map((tool) => {
            const Icon = tool.icon;
            return (
              <Tooltip key={tool.id}>
                <TooltipTrigger asChild>
                  <Button
                    variant={activeTool === tool.id ? "default" : "ghost"}
                    size="icon"
                    className={cn(
                      "h-9 w-9 rounded-full",
                      activeTool === tool.id && "shadow-sm"
                    )}
                    onClick={() => onToolChange(tool.id)}
                  >
                    <Icon size={18} />
                  </Button>
                </TooltipTrigger>
                <TooltipContent>
                  <p className="text-xs">
                    {tool.label} <kbd className="ml-1 text-xs opacity-60">{tool.shortcut}</kbd>
                  </p>
                </TooltipContent>
              </Tooltip>
            );
          })}
        </div>
      </div>
    </TooltipProvider>
  );
}