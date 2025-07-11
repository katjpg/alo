"use client";

import { useState, useCallback, useEffect } from "react";
import { Button } from "@/components/ui/button";
import { IconMessage, IconX } from "@tabler/icons-react";
import { ChatInput } from "@/features/canvas/components/chat-input";
import TreeHierarchy from "@/features/canvas/thought-tree/tree-hierarchy";
import { CanvasToolbar } from "./canvas-toolbar";
import { NodeDetailsPanel } from "./node-details-panel";
import { useRightPanelState } from "@/hooks/use-right-panel";
import { cn } from "@/lib/utils";

export type ToolType = "select" | "hand" | "text" | "node";

export function ThoughtTreeCanvas() {
  const [showChat, setShowChat] = useState(false);
  const [activeTool, setActiveTool] = useState<ToolType>("select");
  const [isMobile, setIsMobile] = useState(false);
  const { isOpen: showRightPanel } = useRightPanelState();

  const toggleChat = useCallback(() => {
    setShowChat(!showChat);
  }, [showChat]);

  // Check for mobile viewport
  useEffect(() => {
    const checkMobile = () => {
      setIsMobile(window.innerWidth <= 640);
    };
    
    checkMobile();
    window.addEventListener('resize', checkMobile);
    return () => window.removeEventListener('resize', checkMobile);
  }, []);

  // Keyboard shortcuts for tool switching
  useEffect(() => {
    const handleKeyPress = (event: KeyboardEvent) => {
      // Prevent shortcuts when typing in input fields
      if (event.target instanceof HTMLInputElement || 
          event.target instanceof HTMLTextAreaElement) {
        return;
      }

      switch (event.key.toLowerCase()) {
        case 'v':
          setActiveTool('select');
          break;
        case 'h':
          setActiveTool('hand');
          break;
        case 't':
          setActiveTool('text');
          break;
        case 'n':
          setActiveTool('node');
          break;
      }
    };

    window.addEventListener('keydown', handleKeyPress);
    return () => window.removeEventListener('keydown', handleKeyPress);
  }, []);

  return (
    <div className="relative h-screen w-full overflow-hidden bg-background">
      {/* Main Canvas - Full Width */}
      <div className="relative w-full h-full">
        {/* Thought Tree Canvas */}
        <TreeHierarchy activeTool={activeTool} onToolChange={setActiveTool} />
        
        {/* Canvas Toolbar */}
        <CanvasToolbar 
          activeTool={activeTool}
          onToolChange={setActiveTool}
        />
      </div>

      {/* Chat Popup - Fixed Position Left */}
      {showChat && (
        <div className={cn(
          "absolute left-4 top-[65%] -translate-y-1/2 z-20",
          "w-80 h-[500px] max-h-[80vh]",
          "bg-background rounded-lg border border-border shadow-2xl",
          "flex flex-col overflow-hidden",
          "transition-all duration-300 ease-in-out",
          isMobile && "w-[calc(100%-2rem)] max-w-sm"
        )}>
          {/* Chat Header */}
          <div className="flex items-center justify-between border-b border-border/50 px-4 py-3 bg-muted/50">
            <h2 className="text-sm font-semibold">ALO</h2>
            <Button
              variant="ghost"
              size="icon"
              onClick={toggleChat}
              className="h-6 w-6"
            >
              <IconX size={16} />
            </Button>
          </div>
          
          {/* Chat Content */}
          <div className="flex-1 overflow-auto p-4 bg-background">
          </div>
          
          {/* Chat Input */}
          <div className="border-t border-border/50 p-2 bg-muted/30">
            <ChatInput />
          </div>
        </div>
      )}
      
      {/* Chat Toggle Button - Bottom Left */}
      <Button
        variant="outline"
        size="icon"
        onClick={toggleChat}
        className={cn(
          "absolute bottom-4 left-4 rounded-full shadow-lg z-10",
          "bg-background/95 backdrop-blur-sm",
          "hover:bg-accent hover:text-accent-foreground",
          "transition-all duration-200",
          showChat && "bg-accent text-accent-foreground"
        )}
      >
        <IconMessage size={20} />
      </Button>

      {/* Right Panel - Pushes content */}
      {showRightPanel && !isMobile && (
        <div className="absolute top-0 right-0 bottom-0 w-[30%] min-w-[280px] bg-background border-l border-border shadow-lg z-10 transition-transform duration-200 transform translate-x-0">
          <NodeDetailsPanel />
        </div>
      )}
    </div>
  );
}