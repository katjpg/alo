"use client"

import React from 'react'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { ProgressIndicator } from './progress-indicator'
import { PropertySlider } from './property-slider'
import { IconBolt, IconInfoCircle } from '@tabler/icons-react'
import { cn } from '@/lib/utils'
import type { MoleculeEditorState } from '@/hooks/use-molecule-editor'
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from '@/components/ui/tooltip'

interface ParameterPanelProps {
  smiles: string
  isValidated: boolean
  isValidating: boolean
  regions: MoleculeEditorState['regions']
  properties: MoleculeEditorState['properties']
  setSmiles: (smiles: string) => void
  validateMolecule: () => Promise<boolean>
  selectRegion: (region: 'core' | 'r1') => void
  removeRegion: (region: 'core' | 'r1') => void
  updateProperty: (property: keyof MoleculeEditorState['properties'], value: number) => void
  toggleProperty: (property: keyof MoleculeEditorState['properties']) => void
  isComplete: () => boolean
  getProgress: () => boolean[]
  debugMode?: boolean
  getSelectionInfo?: () => any
  getSelectionAsSmiles?: () => Promise<any>
}

export function ParameterPanel({
  smiles,
  isValidated,
  isValidating,
  regions,
  properties,
  setSmiles,
  validateMolecule,
  selectRegion,
  removeRegion,
  updateProperty,
  toggleProperty,
  isComplete,
  getProgress,
  debugMode,
  getSelectionInfo,
  getSelectionAsSmiles
}: ParameterPanelProps) {
  const handleStartDiscovery = () => {
    console.log('Starting discovery with:', {
      smiles,
      regions,
      properties
    })
  }

  return (
    <div className="flex flex-col h-full overflow-hidden">
      {/* Header */}
      <div className="flex-shrink-0 p-3 border-b border-border bg-gray-50">
        <h2 className="text-base font-semibold text-gray-900 truncate">Optimization Parameters</h2>
      </div>

      {/* Progress Indicator */}
      <div className="flex-shrink-0 px-3 border-b border-border">
        <ProgressIndicator steps={getProgress()} />
      </div>

      {/* Scrollable Content */}
      <div className="flex-1 overflow-y-auto">
        {/* SMILES Input Section */}
        <div className="p-3 border-b border-border">
          <div className="space-y-4">
            <div>
              <label className="block text-sm font-medium text-gray-700 mb-2">
                SMILES String
              </label>
              <Input
                type="text"
                placeholder="Enter SMILES or draw in editor..."
                value={smiles}
                onChange={(e) => setSmiles(e.target.value)}
                className="font-mono text-sm"
                onKeyPress={(e) => {
                  if (e.key === 'Enter' && !isValidating && smiles.trim()) {
                    e.preventDefault()
                    validateMolecule()
                  }
                }}
              />
            </div>
            <Button 
              onClick={validateMolecule}
              disabled={isValidating || !smiles.trim()}
              className="w-full bg-blue-600 hover:bg-blue-700 disabled:opacity-50 disabled:cursor-not-allowed"
            >
              {isValidating ? (
                <span className="flex items-center justify-center gap-2">
                  <span className="animate-spin rounded-full h-4 w-4 border-b-2 border-white"></span>
                  Validating...
                </span>
              ) : (
                'Validate Molecule'
              )}
            </Button>
            {isValidating && (
              <p className="text-xs text-gray-500 text-center">
                Setting molecule structure in editor...
              </p>
            )}
            {debugMode && isValidated && (
              <>
                <Button
                  onClick={() => {
                    if (getSelectionInfo) {
                      const selectionData = getSelectionInfo()
                      console.log('=== Selection Info Debug ===')
                      console.log(selectionData)
                    }
                  }}
                  variant="outline"
                  className="w-full border-gray-300 text-gray-700 hover:bg-gray-50"
                >
                  Log Selection Info (Debug)
                </Button>
                <Button
                  onClick={async () => {
                    if (getSelectionAsSmiles) {
                      console.log('=== Getting Fragment SMILES ===')
                      const result = await getSelectionAsSmiles()
                      console.log('Fragment SMILES Result:', result)
                    }
                  }}
                  variant="outline"
                  className="w-full border-gray-300 text-gray-700 hover:bg-gray-50"
                >
                  Get Fragment SMILES (Debug)
                </Button>
              </>
            )}
          </div>
        </div>

        {/* Properties Section */}
        <div className="p-3">
          <h3 className="text-base font-semibold text-gray-900 mb-3">Target Properties</h3>
          <div className="space-y-3">
            {Object.entries(properties).map(([key, prop]) => (
              <div key={key} className="space-y-1">
                <div className="flex items-center justify-between">
                  <label className="flex items-center gap-2">
                    <input
                      type="checkbox"
                      checked={prop.enabled}
                      onChange={() => toggleProperty(key as keyof MoleculeEditorState['properties'])}
                      className="rounded border-gray-300 text-blue-600 focus:ring-blue-500"
                    />
                    <span className="text-sm font-medium text-gray-700">{prop.label}</span>
                  </label>
                  <TooltipProvider>
                    <Tooltip>
                      <TooltipTrigger asChild>
                        <IconInfoCircle className="w-4 h-4 text-gray-400 hover:text-gray-600 cursor-help" />
                      </TooltipTrigger>
                      <TooltipContent>
                        <p className="text-sm">{prop.tooltip}</p>
                      </TooltipContent>
                    </Tooltip>
                  </TooltipProvider>
                </div>
                {prop.enabled && (
                  <div className="ml-6">
                    <PropertySlider
                      label=""
                      value={prop.value}
                      min={prop.min}
                      max={prop.max}
                      step={prop.step}
                      unit={prop.unit}
                      onChange={(value) => updateProperty(key as keyof MoleculeEditorState['properties'], value)}
                      tooltip={prop.tooltip}
                    />
                  </div>
                )}
              </div>
            ))}
          </div>
        </div>
      </div>

      {/* Bottom Actions */}
      <div className="flex-shrink-0 p-3 border-t border-border bg-gray-50">
        <Button
          onClick={handleStartDiscovery}
          disabled={!isComplete()}
          className={cn(
            "w-full flex items-center justify-center gap-2",
            isComplete() 
              ? "bg-green-600 hover:bg-green-700" 
              : "bg-gray-400 cursor-not-allowed"
          )}
        >
          <IconBolt size={20} />
          Start Discovery
        </Button>
      </div>
    </div>
  )
}