'use client'

import {
  SidebarGroup,
  SidebarGroupLabel,
  SidebarMenu
} from '@/components/ui/sidebar'
import { useCallback, useEffect, useRef, useState, useTransition } from 'react'
import { HistorySkeleton } from './history-skeleton'
import { MenuItem } from './menu-item'
import { ClearAction } from './clear-action'

interface DiscoverySession {
  id: string
  title: string
  createdAt: string
  smiles?: string
  properties?: Record<string, number>
}

interface DiscoveryPageResponse {
  sessions: DiscoverySession[]
  nextOffset: number | null
}

export function HistorySection() {
  const [sessions, setSessions] = useState<DiscoverySession[]>([])
  const [nextOffset, setNextOffset] = useState<number | null>(null)
  const [isLoading, setIsLoading] = useState(true)
  const loadMoreRef = useRef<HTMLDivElement>(null)
  const [isPending, startTransition] = useTransition()

  const fetchInitialSessions = useCallback(async () => {
    setIsLoading(true)
    try {
      // API CALL here
      await new Promise(resolve => setTimeout(resolve, 1000))
      
      // Using mock data as placeholder
      const mockSessions: DiscoverySession[] = []
      
      setSessions(mockSessions)
      setNextOffset(null)
    } catch (error) {
      console.error('Failed to load initial sessions:', error)
      setNextOffset(null)
    } finally {
      setIsLoading(false)
    }
  }, [])

  useEffect(() => {
    fetchInitialSessions()
  }, [fetchInitialSessions])

  useEffect(() => {
    const handleHistoryUpdate = () => {
      startTransition(() => {
        fetchInitialSessions()
      })
    }
    window.addEventListener('discovery-history-updated', handleHistoryUpdate)
    return () => {
      window.removeEventListener('discovery-history-updated', handleHistoryUpdate)
    }
  }, [fetchInitialSessions])

  const fetchMoreSessions = useCallback(async () => {
    if (isLoading || nextOffset === null) return

    setIsLoading(true)
    try {
      // API Call 
      await new Promise(resolve => setTimeout(resolve, 500))
      
      // Mock implementation
      setSessions(prevSessions => [...prevSessions])
      setNextOffset(null)
    } catch (error) {
      console.error('Failed to load more sessions:', error)
      setNextOffset(null)
    } finally {
      setIsLoading(false)
    }
  }, [nextOffset, isLoading])

  useEffect(() => {
    const observerRefValue = loadMoreRef.current
    if (!observerRefValue || nextOffset === null || isPending) return

    const observer = new IntersectionObserver(
      ([entry]) => {
        if (entry.isIntersecting && !isLoading && !isPending) {
          fetchMoreSessions()
        }
      },
      { threshold: 0.1 }
    )

    observer.observe(observerRefValue)

    return () => {
      if (observerRefValue) {
        observer.unobserve(observerRefValue)
      }
    }
  }, [fetchMoreSessions, nextOffset, isLoading, isPending])

  const isHistoryEmpty = !isLoading && !sessions.length && nextOffset === null

  return (
    <div className="flex flex-col flex-1 h-full">
      <SidebarGroup>
        <div className="flex items-center justify-between w-full">
          <SidebarGroupLabel className="p-0">Discovery History</SidebarGroupLabel>
          <ClearAction empty={isHistoryEmpty} />
        </div>
      </SidebarGroup>
      <div className="flex-1 overflow-y-auto mb-2 relative">
        {isHistoryEmpty && !isPending ? (
          <div className="px-2 text-foreground/30 text-sm text-center py-4">
            No discovery sessions yet
          </div>
        ) : (
          <SidebarMenu>
            {sessions.map(
              (session: DiscoverySession) => session && <MenuItem key={session.id} session={session} />
            )}
          </SidebarMenu>
        )}
        <div ref={loadMoreRef} style={{ height: '1px' }} />
        {(isLoading || isPending) && (
          <div className="py-2">
            <HistorySkeleton />
          </div>
        )}
      </div>
    </div>
  )
}