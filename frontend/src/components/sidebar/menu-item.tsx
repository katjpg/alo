'use client'

import {
  SidebarMenuItem,
  SidebarMenuButton
} from '@/components/ui/sidebar'
import { IconAtom, IconDots } from '@tabler/icons-react'
import Link from 'next/link'
import { cn } from '@/lib/utils'

interface DiscoverySession {
  id: string
  title: string
  createdAt: string
  smiles?: string
  properties?: Record<string, number>
}

interface DiscoveryMenuItemProps {
  session: DiscoverySession
}

export function MenuItem({ session }: DiscoveryMenuItemProps) {
  const formatDate = (dateString: string) => {
    const date = new Date(dateString)
    const now = new Date()
    const diffInHours = Math.floor((now.getTime() - date.getTime()) / (1000 * 60 * 60))
    
    if (diffInHours < 1) return 'Just now'
    if (diffInHours < 24) return `${diffInHours}h ago`
    if (diffInHours < 48) return 'Yesterday'
    return date.toLocaleDateString()
  }

  return (
    <SidebarMenuItem>
      <SidebarMenuButton asChild className="group">
        <Link href={`/discovery/${session.id}`} className="flex items-center justify-between w-full">
          <div className="flex items-center gap-2 min-w-0 flex-1">
            <IconAtom className="size-4 flex-shrink-0" />
            <div className="min-w-0 flex-1">
              <div className="truncate text-sm font-medium">
                {session.title}
              </div>
              <div className="text-xs text-muted-foreground">
                {formatDate(session.createdAt)}
              </div>
            </div>
          </div>
          <button 
            className="opacity-0 group-hover:opacity-100 transition-opacity p-1 hover:bg-accent rounded"
            onClick={(e) => {
              e.preventDefault()
              e.stopPropagation()
              // Handle menu actions
            }}
          >
            <IconDots className="size-3" />
          </button>
        </Link>
      </SidebarMenuButton>
    </SidebarMenuItem>
  )
}