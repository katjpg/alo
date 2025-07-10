"use client"

import dynamic from 'next/dynamic'
import { Suspense } from 'react'
import type { MoleculeEditorState } from '@/hooks/use-molecule-editor'

const KetcherEditor = dynamic(
  () => import('./ketcher-editor').then(mod => mod.KetcherEditor),
  { 
    ssr: false,
    loading: () => (
      <div className="flex-1 flex items-center justify-center">
        <div className="text-center">
          <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-600 mx-auto mb-4"></div>
          <p className="text-gray-600">Loading molecular editor...</p>
        </div>
      </div>
    )
  }
)

interface KetcherEditorWrapperProps {
  isValidated: boolean
  selectedRegion: MoleculeEditorState['selectedRegion']
  regions: MoleculeEditorState['regions']
  mockSelections: MoleculeEditorState['mockSelections']
  selectRegion: (region: 'core' | 'r1') => void
  defineRegion: (x: number, y: number) => void
  onKetcherInit?: (ketcher: any) => void
  smiles?: string
}

export function KetcherEditorWrapper(props: KetcherEditorWrapperProps) {
  return (
    <Suspense fallback={
      <div className="flex-1 flex items-center justify-center">
        <div className="text-center">
          <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-600 mx-auto mb-4"></div>
          <p className="text-gray-600">Loading molecular editor...</p>
        </div>
      </div>
    }>
      <KetcherEditor {...props} />
    </Suspense>
  )
}