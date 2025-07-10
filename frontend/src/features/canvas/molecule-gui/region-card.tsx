"use client"

import { Button } from '@/components/ui/button'

interface RegionCardProps {
  region: 'core' | 'r1'
  value: string
  color: string
  onEdit: () => void
  onDelete: () => void
}

export function RegionCard({ region, value, color, onEdit, onDelete }: RegionCardProps) {
  return (
    <div className="bg-gray-50 border border-gray-200 rounded-md p-3 flex items-center justify-between">
      <div className="flex items-center gap-3">
        <span 
          className="w-4 h-4 rounded-full"
          style={{ backgroundColor: color }}
        />
        <div>
          <div className="text-sm font-medium capitalize">{region}</div>
          <div className="text-xs text-gray-600 font-mono">{value}</div>
        </div>
      </div>
      
      <div className="flex gap-2">
        <Button
          variant="outline"
          size="sm"
          onClick={onEdit}
          className="h-7 px-2 text-xs"
        >
          Edit
        </Button>
        <Button
          variant="outline"
          size="sm"
          onClick={onDelete}
          className="h-7 px-2 text-xs hover:text-red-600 hover:border-red-600"
        >
          Delete
        </Button>
      </div>
    </div>
  )
}