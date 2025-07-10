"use client"

import { Button } from '@/components/ui/button'
import { cn } from '@/lib/utils'

interface KetcherHeaderProps {
  selectedRegion: 'core' | 'r1' | null
  regions: {
    core: string | null
    r1: string | null
  }
  selectRegion: (region: 'core' | 'r1') => void
  isValidated: boolean
}

export function KetcherHeader({
  selectedRegion,
  regions,
  selectRegion,
  isValidated
}: KetcherHeaderProps) {
  const getRegionColor = (region: 'core' | 'r1') => {
    const colors = {
      core: '#e74c3c',
      r1: '#f39c12'
    }
    return colors[region]
  }
  
  const getRegionPreview = (region: 'core' | 'r1', regionData: string) => {
    try {
      const data = JSON.parse(regionData)
      return data.description || `${data.atomIds?.length || 0} atoms selected`
    } catch {
      return regionData
    }
  }

  return (
    <div className="p-3 border-b border-border bg-gray-50">
      <h2 className="text-base font-semibold text-gray-900 mb-2">Molecule Editor</h2>
      <div className="flex flex-wrap items-center gap-2">
        <span className="text-sm text-gray-600 w-full sm:w-auto">
          Select atoms in the editor, then click a button to save as region:
        </span>
        
        <Button
          variant="outline"
          size="sm"
          onClick={() => selectRegion('core')}
          disabled={!isValidated}
          className={cn(
            "flex items-center gap-2",
            regions.core && "bg-green-50 border-green-500 text-green-700"
          )}
          title={regions.core ? `Core region saved: ${getRegionPreview('core', regions.core)}` : 'Click to save selected atoms as Core region'}
        >
          <span 
            className="w-3 h-3 rounded-full border border-gray-300"
            style={{ backgroundColor: getRegionColor('core') }}
          />
          Core {regions.core && '✓'}
        </Button>
        
        <Button
          variant="outline"
          size="sm"
          onClick={() => selectRegion('r1')}
          disabled={!isValidated}
          className={cn(
            "flex items-center gap-2",
            regions.r1 && "bg-green-50 border-green-500 text-green-700"
          )}
          title={regions.r1 ? `R1 region saved: ${getRegionPreview('r1', regions.r1)}` : 'Click to save selected atoms as R1 region'}
        >
          <span 
            className="w-3 h-3 rounded-full border border-gray-300"
            style={{ backgroundColor: getRegionColor('r1') }}
          />
          R1 {regions.r1 && '✓'}
        </Button>
      </div>
    </div>
  )
}