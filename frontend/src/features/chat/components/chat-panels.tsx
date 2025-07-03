"use client"

import * as React from "react"
import {
  ResizableHandle,
  ResizablePanel,
  ResizablePanelGroup,
} from "@/components/ui/resizable"
import { TooltipProvider } from "@/components/ui/tooltip"
import { Button } from "@/components/ui/button"
import { Tabs, TabsList, TabsTrigger, TabsContent } from "@/components/ui/tabs"
import { IconX, IconArrowBarLeft } from "@tabler/icons-react"
import { ChatInput } from "./chat-input"
import { useRightPanelState } from "@/hooks/use-right-panel"

interface ChatPanelsProps {
  defaultLayout?: number[]
  children: React.ReactNode
}

export function ChatPanels({
  defaultLayout = [26, 52, 22],
  children,
}: ChatPanelsProps) {
  const [layoutMode, setLayoutMode] = React.useState<'desktop' | 'chat-only'>('desktop')
  const [isInitialized, setIsInitialized] = React.useState(false)
  const [activeTab, setActiveTab] = React.useState('chat')
  const { isOpen: isRightPanelOpen, toggle: toggleRightPanel } = useRightPanelState()

  React.useEffect(() => {
    const updateLayout = () => {
      const width = window.innerWidth
      
      // Fixed viewport width breakpoints
      if (width <= 640) {
        setLayoutMode('chat-only')
      } else {
        setLayoutMode('desktop')
      }
      
      // Mark as initialized after first layout calculation
      if (!isInitialized) {
        setIsInitialized(true)
      }
    }
    
    updateLayout()
    window.addEventListener("resize", updateLayout)
    return () => window.removeEventListener("resize", updateLayout)
  }, [])

  // Don't render anything until layout is initialized
  if (!isInitialized) {
    return null
  }


  // Chat-only layout - mobile tabs
  if (layoutMode === 'chat-only') {
    return (
      <TooltipProvider delayDuration={0}>
        <div className="relative flex h-screen w-full flex-col overflow-hidden bg-background overscroll-behavior-none">
          <Tabs value={activeTab} onValueChange={setActiveTab} className="flex h-full flex-col">
            <TabsContent value="chat" className="flex-1 flex flex-col m-0 data-[state=inactive]:hidden">
              <div className="flex items-center justify-between border-b border-border/50 px-4 py-3">
                <h2 className="text-sm font-semibold">Chat</h2>
                <div className="h-6 w-6" />
              </div>
              
              <div className="flex-1 overflow-auto p-4">
              </div>
              
              <div className="border-t border-border/50 p-4">
                <ChatInput mode="chat" />
              </div>
            </TabsContent>

            <TabsContent value="thought-tree" className="flex-1 flex flex-col m-0 data-[state=inactive]:hidden">
              <div className="flex items-center justify-between border-b border-border/50 px-4 py-3">
                <h2 className="text-sm font-semibold">Thought Tree</h2>
                <div className="h-6 w-6" />
              </div>
              
              <div className="flex-1 overflow-hidden">
                {children}
              </div>
            </TabsContent>

            <TabsList className="absolute bottom-8 left-1/2 transform -translate-x-1/2 z-20">
              <TabsTrigger value="chat">Chat</TabsTrigger>
              <TabsTrigger value="thought-tree">Thought Tree</TabsTrigger>
            </TabsList>
          </Tabs>
        </div>
      </TooltipProvider>
    )
  }



  // Desktop layout - resizable left/middle, overlay right
  return (
    <TooltipProvider delayDuration={0}>
      <div className="relative flex h-screen overflow-hidden bg-background overscroll-behavior-none">
        <ResizablePanelGroup direction="horizontal" className="flex h-full">
          {/* Left Panel - Resizable */}
          <ResizablePanel
            defaultSize={defaultLayout[0]}
            minSize={30}
            maxSize={55}
            className="bg-background"
          >
            <div className="flex h-full flex-col">
              {/* Left Panel Header */}
              <div className="flex items-center justify-between border-b border-border/50 px-4 py-3">
                <h2 className="text-sm font-semibold">Chat</h2>
                <div className="h-6 w-6" />
              </div>
              
              {/* Left Panel Content */}
              <div className="flex-1 overflow-auto p-4">
              </div>
              
              {/* Chat Input at bottom */}
              <div className="border-t border-border/50 p-4">
                <ChatInput mode="chat" />
              </div>
            </div>
          </ResizablePanel>

          <ResizableHandle withHandle />

          {/* Middle Panel - Resizable */}
          <ResizablePanel
            defaultSize={defaultLayout[1]}
            minSize={25}
            maxSize={75}
            className="bg-background"
          >
            <div className="flex h-full flex-col">
              {/* Middle Panel Header */}
              <div className="flex items-center justify-between border-b border-border/50 pl-4 pr-3 py-3">
                <h2 className="text-sm font-semibold">Thought Tree</h2>
                
                {/* Show expand button when right panel is collapsed, invisible spacer when open */}
                {!isRightPanelOpen ? (
                  <Button
                    variant="ghost"
                    size="sm"
                    onClick={toggleRightPanel}
                    className="h-6 w-6 p-0 ml-4"
                  >
                    <IconArrowBarLeft size={16} strokeWidth={2} />
                  </Button>
                ) : (
                  <div className="h-6 w-6" />
                )}
              </div>
              
              {/* Middle Panel Content - Thought Tree */}
              <div className="flex-1 overflow-hidden">
                {children}
              </div>
            </div>
          </ResizablePanel>
        </ResizablePanelGroup>

        {/* Right Panel - Overlay */}
        {isInitialized && isRightPanelOpen && layoutMode === 'desktop' && (
          <div className="absolute top-0 right-0 bottom-0 w-[280px] bg-background border-l border-border shadow-lg z-10 transition-transform duration-200 transform translate-x-0">
            <div className="flex h-full flex-col">
              {/* Right Panel Header */}
              <div className="flex items-center justify-between border-b border-border/50 pl-4 pr-3 py-3">
                <h2 className="text-sm font-semibold">Details</h2>
                <Button
                  variant="ghost"
                  size="sm"
                  onClick={toggleRightPanel}
                  className="h-6 w-6 p-0 ml-1"
                >
                  <IconX size={16} strokeWidth={2} />
                </Button>
              </div>
              
              {/* Right Panel Content */}
              <div className="flex-1 overflow-auto p-4">
                <div className="space-y-4">
                  <div className="rounded-lg border p-4">
                    <p className="text-sm text-muted-foreground">
                      Details panel
                    </p>
                    <p className="text-xs text-muted-foreground mt-2">
                      This panel shows additional information and properties
                    </p>
                  </div>
                </div>
              </div>
            </div>
          </div>
        )}
      </div>
    </TooltipProvider>
  )
}