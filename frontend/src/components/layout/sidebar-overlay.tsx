"use client"

import { useEffect } from 'react'
import { useSidebarState } from '@/hooks/use-sidebar'
import { cn } from '@/lib/utils'

interface SidebarOverlayProps {
  className?: string
}

export function SidebarOverlay({ className }: SidebarOverlayProps = {}) {
  const { isOpen, isHydrated, setIsOpen } = useSidebarState()

  useEffect(() => {
    const handleKeyDown = (event: KeyboardEvent) => {
      if (event.key === 'Escape' && isOpen) {
        setIsOpen(false)
      }
    }

    document.addEventListener('keydown', handleKeyDown)
    return () => document.removeEventListener('keydown', handleKeyDown)
  }, [isOpen, setIsOpen])

  if (!isHydrated || !isOpen) return null

  return (
    <div
      className={cn(
        "fixed inset-0 z-40 bg-black/10",
        "transition-opacity duration-200 ease-in-out",
        "animate-in fade-in-0",
        className
      )}
      onClick={() => setIsOpen(false)}
      aria-hidden="true"
      aria-label="Close sidebar"
    />
  )
}