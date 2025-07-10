"use client"

import { useEffect, useRef, useState } from 'react'
import { KetcherHeader } from './ketcher-header'
import type { MoleculeEditorState } from '@/hooks/use-molecule-editor'
import '@/styles/ketcher.css'

interface KetcherEditorProps {
  isValidated: boolean
  selectedRegion: MoleculeEditorState['selectedRegion']
  regions: MoleculeEditorState['regions']
  mockSelections: MoleculeEditorState['mockSelections']
  selectRegion: (region: 'core' | 'r1') => void
  defineRegion: (x: number, y: number) => void
  onKetcherInit?: (ketcher: any) => void
  smiles?: string
}

export function KetcherEditor({
  isValidated,
  selectedRegion,
  regions,
  mockSelections,
  selectRegion,
  defineRegion,
  onKetcherInit,
  smiles
}: KetcherEditorProps) {
  const containerRef = useRef<HTMLDivElement>(null)
  const [isReady, setIsReady] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [Editor, setEditor] = useState<any>(null)
  const [provider, setProvider] = useState<any>(null)

  // Initialize Ketcher
  useEffect(() => {
    const loadKetcher = async () => {
      try {
        // Dynamically import Ketcher components
        const [{ Editor: KetcherEditor }, { StandaloneStructServiceProvider }] = await Promise.all([
          import('ketcher-react'),
          import('ketcher-standalone')
        ])

        // Cast to any to avoid TypeScript issues (following the example)
        const structProvider = new (StandaloneStructServiceProvider as any)()
        
        setEditor(() => KetcherEditor)
        setProvider(structProvider)
        setIsReady(true)
      } catch (err) {
        console.error('Failed to load Ketcher:', err)
        setError('Failed to load molecular editor')
      }
    }

    loadKetcher()
  }, [])

  const handleKetcherInit = (ketcher: any) => {
    console.log('Ketcher initialized')
    
    // Store globally for debugging
    // @ts-ignore
    window.ketcher = ketcher
    
    // Call the parent callback
    if (onKetcherInit) {
      onKetcherInit(ketcher)
    }
  }

  // Ketcher button configuration  
  const buttonsConfig = {
    // Hide reaction tools
    'reaction-arrows': { hidden: true },
    'reaction-arrow-open-angle': { hidden: true },
    'reaction-arrow-filled-triangle': { hidden: true },
    'reaction-arrow-filled-bow': { hidden: true },
    'reaction-arrow-dashed-open-angle': { hidden: true },
    'reaction-arrow-failed': { hidden: true },
    'reaction-arrow-both-ends-filled-triangle': { hidden: true },
    'reaction-arrow-equilibrium-filled-half-bow': { hidden: true },
    'reaction-arrow-equilibrium-filled-triangle': { hidden: true },
    'reaction-arrow-equilibrium-open-angle': { hidden: true },
    'reaction-arrow-unbalanced-equilibrium-filled-half-bow': { hidden: true },
    'reaction-arrow-unbalanced-equilibrium-open-half-angle': { hidden: true },
    'reaction-arrow-unbalanced-equilibrium-large-filled-half-bow': { hidden: true },
    'reaction-arrow-unbalanced-equilibrium-filled-half-triangle': { hidden: true },
    'reaction-plus': { hidden: true },
    'reaction-mapping-tools': { hidden: true },
    'reaction-map': { hidden: true },
    'reaction-unmap': { hidden: true },
    'reaction-automap': { hidden: true },
    // Hide R-group tools
    'rgroup': { hidden: true },
    'rgroup-label': { hidden: true },
    'rgroup-fragment': { hidden: true },
    'rgroup-attpoints': { hidden: true },
    // Hide shape tools
    'shape': { hidden: true },
    'shape-ellipse': { hidden: true },
    'shape-rectangle': { hidden: true },
    'shape-line': { hidden: true },
    // Hide other tools
    'text': { hidden: true },
    'enhanced-stereo': { hidden: true },
    'sgroup': { hidden: true },
    'sgroup-data': { hidden: true }
  }

  return (
    <>
      <KetcherHeader
        selectedRegion={selectedRegion}
        regions={regions}
        selectRegion={selectRegion}
        isValidated={isValidated}
      />
      
      <div 
        ref={containerRef}
        className="flex-1 relative overflow-hidden ketcher-container" 
        style={{ minHeight: '500px' }}
      >
        {error ? (
          <div className="ketcher-error">
            <p className="text-red-600 mb-2">{error}</p>
            <button
              onClick={() => window.location.reload()}
              className="px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700"
            >
              Reload Page
            </button>
          </div>
        ) : !isReady ? (
          <div className="ketcher-loading">
            <div className="text-center">
              <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-600 mx-auto mb-4"></div>
              <p className="text-gray-600">Loading molecular editor...</p>
            </div>
          </div>
        ) : (
          <>
            {Editor && provider && (
              <Editor
                staticResourcesUrl={typeof window !== 'undefined' ? window.location.origin : ''}
                structServiceProvider={provider}
                errorHandler={(message: string) => {
                  // Only log errors we care about
                  if (!message.includes('KetcherLogger') && 
                      !message.includes('SVGLength') &&
                      !message.includes('relative length')) {
                    console.error('Ketcher error:', message)
                  }
                }}
                onInit={handleKetcherInit}
                buttons={buttonsConfig}
              />
            )}
            
            {/* Overlay for mock selection when a region is selected */}
            {selectedRegion && (
              <div 
                className="absolute inset-0 cursor-crosshair z-10"
                onClick={(e) => {
                  const rect = e.currentTarget.getBoundingClientRect()
                  const x = e.clientX - rect.left
                  const y = e.clientY - rect.top
                  defineRegion({ x, y })
                }}
                style={{ pointerEvents: selectedRegion ? 'auto' : 'none' }}
              />
            )}
            
            {selectedRegion && (
              <div className="absolute bottom-3 left-1/2 transform -translate-x-1/2 bg-gray-900 text-white px-4 py-2 rounded-md text-sm whitespace-nowrap z-10 pointer-events-none">
                Click to select atoms for {selectedRegion}
              </div>
            )}
            
            {/* Mock selections visualization */}
            {mockSelections.map((selection, index) => (
              <div
                key={index}
                className="absolute border-2 border-dashed pointer-events-none"
                style={{
                  left: `${selection.x}px`,
                  top: `${selection.y}px`,
                  width: `${selection.width}px`,
                  height: `${selection.height}px`,
                  borderColor: selection.region === 'core' ? '#e74c3c' : '#f39c12',
                  backgroundColor: selection.region === 'core' ? 'rgba(231, 76, 60, 0.1)' : 'rgba(243, 156, 18, 0.1)'
                }}
              />
            ))}
          </>
        )}
      </div>
    </>
  )
}