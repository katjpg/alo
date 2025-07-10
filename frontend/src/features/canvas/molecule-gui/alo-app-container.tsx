"use client"

import { useMoleculeEditor } from '@/hooks/use-molecule-editor'
import { KetcherEditorWrapper } from './ketcher-editor-wrapper'
import { ParameterPanel } from './parameter-panel'

export function ALOAppContainer() {
  const moleculeEditor = useMoleculeEditor()

  return (
    <div className="flex flex-col md:flex-row h-full w-full overflow-hidden" style={{ height: '600px' }}>
      {/* Left Panel - Molecule Editor */}
      <div className="h-[65%] md:h-full md:flex-1 bg-white border-b md:border-b-0 md:border-r border-border flex flex-col overflow-hidden" style={{ minHeight: '400px' }}>
        <KetcherEditorWrapper 
          {...moleculeEditor} 
          onKetcherInit={moleculeEditor.handleKetcherInit}
          smiles={moleculeEditor.smiles}
        />
      </div>

      {/* Right Panel - Parameters */}
      <div className="h-[35%] md:h-full w-full md:w-[300px] lg:w-[320px] xl:w-[360px] bg-white flex flex-col overflow-hidden">
        <ParameterPanel {...moleculeEditor} />
      </div>
    </div>
  )
}