"use client"

import { useEffect } from 'react'
import { 
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from '@/components/ui/tooltip'
import { Button } from '@/components/ui/button'
import { cn } from '@/lib/utils'
import { IconPlus, IconAtom, IconArrowBarLeft, IconArrowBarRight, IconBookmark } from '@tabler/icons-react'
import Link from 'next/link'
import { Suspense } from 'react'
import { HistorySection } from '@/components/sidebar/history-section'
import { HistorySkeleton } from '@/components/sidebar/history-skeleton'
import { useSidebarState } from '@/hooks/use-sidebar'

export function AppSidebar() {
  const { isOpen, isHydrated, setHydrated, toggle } = useSidebarState()

  useEffect(() => {
    setHydrated(true)
  }, [setHydrated])

  if (!isHydrated) {
    return (
      <TooltipProvider>
        <div className="fixed bottom-0 left-0 top-0 z-[50] flex h-screen flex-shrink-0 flex-col border-r border-border bg-background w-[60px]">
          <div className="flex flex-row justify-between items-center border-b border-border/50 px-2 py-3">
            <div className="h-6 w-6" />
            <div className="flex items-center justify-center w-full">
              <Link href="/" className="flex items-center justify-center">
                <IconAtom className="size-5" />
              </Link>
            </div>
            <div className="h-6 w-6" />
          </div>
          <div className="flex flex-col h-full">
            <div className="px-2 py-4 space-y-3">
              <Button variant="outline" size="icon" className="w-full">
                <IconPlus className="size-4" />
              </Button>
              <Button variant="ghost" size="icon" className="w-full">
                <IconBookmark className="size-4" />
              </Button>
            </div>
            <div className="mt-auto p-2 border-t border-border/50">
              <Button variant="ghost" size="icon" className="w-full">
                <IconArrowBarRight size={16} strokeWidth={2} />
              </Button>
            </div>
          </div>
        </div>
      </TooltipProvider>
    )
  }

  return (
    <TooltipProvider>
      <div
        className={cn(
          "fixed bottom-0 left-0 top-0 z-[50] flex h-screen flex-shrink-0 flex-col border-r border-border transition-all duration-200 bg-background",
          isOpen
            ? "w-[240px]"
            : "w-[60px]"
        )}
      >
        <div className="flex flex-row justify-between items-center border-b border-border/50 px-2 py-3">
          {isOpen ? (
            <>
              <Link href="/" className="flex items-center gap-2">
                <IconAtom className="size-5" />
                <span className="font-semibold text-sm">ALO</span>
              </Link>
              <div className="h-6 w-6" />
            </>
          ) : (
            <>
              <div className="h-6 w-6" />
              <div className="flex items-center justify-center w-full">
                <Tooltip>
                  <TooltipTrigger asChild>
                    <Link href="/" className="flex items-center justify-center">
                      <IconAtom className="size-5" />
                    </Link>
                  </TooltipTrigger>
                  <TooltipContent side="right">
                    <p>ALO - AI-based Ligand Optimization</p>
                  </TooltipContent>
                </Tooltip>
              </div>
              <div className="h-6 w-6" />
            </>
          )}
        </div>

        <div className="flex flex-col h-full">
          <div className="px-2 py-4 space-y-3">
            {isOpen ? (
              <Button asChild variant="outline" className="w-full justify-start">
                <Link href="/" className="flex items-center gap-2">
                  <IconPlus className="size-4" />
                  <span>New Discovery</span>
                </Link>
              </Button>
            ) : (
              <Tooltip>
                <TooltipTrigger asChild>
                  <Button asChild variant="outline" size="icon" className="w-full">
                    <Link href="/" className="flex items-center justify-center">
                      <IconPlus className="size-4" />
                    </Link>
                  </Button>
                </TooltipTrigger>
                <TooltipContent side="right">
                  <p>New Discovery</p>
                </TooltipContent>
              </Tooltip>
            )}

            {isOpen ? (
              <Button asChild variant="ghost" className="w-full justify-start">
                <Link href="/artifacts" className="flex items-center gap-2">
                  <IconBookmark className="size-4" />
                  <span>Saved artifacts</span>
                </Link>
              </Button>
            ) : (
              <Tooltip>
                <TooltipTrigger asChild>
                  <Button asChild variant="ghost" size="icon" className="w-full">
                    <Link href="/artifacts" className="flex items-center justify-center">
                      <IconBookmark className="size-4" />
                    </Link>
                  </Button>
                </TooltipTrigger>
                <TooltipContent side="right">
                  <p>Saved artifacts</p>
                </TooltipContent>
              </Tooltip>
            )}
          </div>
          
          {isOpen && (
            <div className="flex-1 overflow-y-auto px-2">
              <Suspense fallback={<HistorySkeleton />}>
                <HistorySection />
              </Suspense>
            </div>
          )}
          
          <div className="mt-auto p-2 border-t border-border/50">
            {isOpen ? (
              <Button
                variant="ghost"
                size="sm"
                onClick={toggle}
                className="w-full justify-start"
              >
                <IconArrowBarLeft size={16} strokeWidth={2} />
                <span className="ml-2">Close</span>
              </Button>
            ) : (
              <Tooltip>
                <TooltipTrigger asChild>
                  <Button
                    variant="ghost"
                    size="icon"
                    onClick={toggle}
                    className="w-full"
                  >
                    <IconArrowBarRight size={16} strokeWidth={2} />
                  </Button>
                </TooltipTrigger>
                <TooltipContent side="right">
                  <p>Expand Sidebar</p>
                </TooltipContent>
              </Tooltip>
            )}
          </div>
        </div>
      </div>
    </TooltipProvider>
  )
}