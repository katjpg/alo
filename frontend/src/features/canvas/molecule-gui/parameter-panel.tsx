"use client"

import React from 'react'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { PropertySlider } from './property-slider'
import { PropertySelector } from '../components/property-selector'
import { IconBolt } from '@tabler/icons-react'
import { cn } from '@/lib/utils'
import type { MoleculeEditorState } from '@/hooks/use-molecule-editor'

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
  getProgress
}: ParameterPanelProps) {
  const handleStartDiscovery = () => {
    // Start discovery process
    // TODO: Implement discovery logic
  }

  return (
    <div className="flex flex-col h-full overflow-hidden">
      {/* Header */}
      <div className="flex-shrink-0 p-3 border-b border-border bg-gray-50">
        <h2 className="text-base font-semibold text-gray-900 truncate">Optimization Parameters</h2>
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
          </div>
        </div>

        {/* Properties Section */}
        <div className="p-3">
          <h3 className="text-base font-semibold text-gray-900 mb-3">Target Properties</h3>
          
          {/* Property Selector Dropdown */}
          <div className="mb-4">
            <PropertySelector
              selectedProperties={Object.entries(properties)
                .filter(([_, prop]) => prop.enabled)
                .map(([key, _]) => {
                  // Map state keys to PropertySelector IDs
                  const keyMap: Record<string, string> = {
                    'bbbp': 'BBBP',
                    'hia': 'HIA',
                    'mutag': 'Mutag',
                    'drd2': 'DRD2',
                    'plogP': 'plogP',
                    'qed': 'QED'
                  }
                  return keyMap[key] || key
                })}
              onPropertyChange={(propertyId, checked) => {
                // Map PropertySelector IDs to state keys
                const idMap: Record<string, keyof MoleculeEditorState['properties']> = {
                  'BBBP': 'bbbp',
                  'HIA': 'hia',
                  'Mutag': 'mutag',
                  'DRD2': 'drd2',
                  'plogP': 'plogP',
                  'QED': 'qed'
                }
                const propKey = idMap[propertyId]
                if (propKey) {
                  toggleProperty(propKey)
                }
              }}
            />
          </div>
          
          {/* Selected Properties with Sliders */}
          <div className="space-y-4">
            {Object.entries(properties)
              .filter(([_, prop]) => prop.enabled)
              .map(([key, prop]) => (
                <div key={key} className="space-y-2">
                  <div className="flex items-center justify-between">
                    <span className="text-sm font-medium text-gray-700">{prop.label}</span>
                    <span className="text-sm font-semibold text-gray-900">
                      {prop.value.toFixed(prop.step < 1 ? 1 : 0)}{prop.unit}
                    </span>
                  </div>
                  <PropertySlider
                    label=""
                    value={prop.value}
                    min={prop.min}
                    max={prop.max}
                    step={prop.step}
                    onChange={(value) => updateProperty(key as keyof MoleculeEditorState['properties'], value)}
                  />
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