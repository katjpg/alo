"use client"

import { useState, useCallback } from 'react'

export interface PropertyConfig {
  value: number
  enabled: boolean
  min: number
  max: number
  step: number
  unit?: string
  label: string
  tooltip: string
}

export interface MoleculeEditorState {
  smiles: string
  isValidated: boolean
  selectedRegion: 'core' | 'r1' | null
  regions: {
    core: string | null
    r1: string | null
  }
  properties: {
    bbbp: PropertyConfig
    hia: PropertyConfig
    mutag: PropertyConfig
    drd2: PropertyConfig
    plogP: PropertyConfig
    qed: PropertyConfig
  }
  mockSelections: Array<{
    x: number
    y: number
    width: number
    height: number
    region: string
  }>
  ketcherInstance: any | null
}

export function useMoleculeEditor() {
  const [smiles, setSmiles] = useState('')
  const [isValidated, setIsValidated] = useState(false)
  const [isValidating, setIsValidating] = useState(false)
  const [selectedRegion, setSelectedRegion] = useState<'core' | 'r1' | null>(null)
  const [regions, setRegions] = useState<MoleculeEditorState['regions']>({
    core: null,
    r1: null
  })
  const [debugMode] = useState(() => {
    // Debug mode via URL parameter / localStorage
    if (typeof window !== 'undefined') {
      const urlParams = new URLSearchParams(window.location.search)
      return urlParams.get('debug') === 'true' || localStorage.getItem('ketcherDebug') === 'true'
    }
    return false
  })
  const [properties, setProperties] = useState<MoleculeEditorState['properties']>({
    bbbp: {
      value: 0.8,
      enabled: true,
      min: 0,
      max: 1,
      step: 0.01,
      label: 'BBBP',
      tooltip: 'Blood-Brain Barrier Permeability'
    },
    hia: {
      value: 0.9,
      enabled: true,
      min: 0,
      max: 1,
      step: 0.01,
      label: 'HIA',
      tooltip: 'Human Intestinal Absorption'
    },
    mutag: {
      value: 0.1,
      enabled: false,
      min: 0,
      max: 1,
      step: 0.01,
      label: 'Mutag',
      tooltip: 'Mutagenicity'
    },
    drd2: {
      value: 0.7,
      enabled: true,
      min: 0,
      max: 1,
      step: 0.01,
      label: 'DRD2',
      tooltip: 'Dopamine Receptor D2 Activity'
    },
    plogP: {
      value: 2.5,
      enabled: false,
      min: -2,
      max: 7,
      step: 0.1,
      label: 'plogP',
      tooltip: 'Predicted Lipophilicity'
    },
    qed: {
      value: 0.7,
      enabled: true,
      min: 0,
      max: 1,
      step: 0.01,
      label: 'QED',
      tooltip: 'Quantitative Estimate of Drug-likeness'
    }
  })
  const [mockSelections, setMockSelections] = useState<MoleculeEditorState['mockSelections']>([])
  const [ketcherInstance, setKetcherInstance] = useState<any | null>(null)

  const validateMolecule = useCallback(async () => {
    // Check if we have a SMILES string to validate
    const smilesString = smiles.trim()
    if (!smilesString) {
      console.log('No SMILES to validate')
      return false
    }
    
    // Check if Ketcher instance exists
    if (!ketcherInstance) {
      console.warn('Cannot validate molecule: Ketcher not initialized')
      return false
    }
    
    setIsValidating(true)
    
    try {
      console.log('Setting molecule in Ketcher:', smilesString)
      
      // Simply call setMolecule - trust that it works
      await ketcherInstance.setMolecule(smilesString)
      
      // Wait a bit for rendering
      await new Promise(resolve => setTimeout(resolve, 300))
      
      // Mark as validated
      setIsValidated(true)
      setIsValidating(false)
      return true
      
    } catch (error) {
      console.error('Error setting molecule:', error)
      setIsValidating(false)
      return false
    }
  }, [smiles, ketcherInstance])

  const selectRegion = useCallback((region: 'core' | 'r1') => {
    // Disabled for now - regions functionality not in use
    return;
    
    // Original code remains below...
    if (!isValidated) return
    setSelectedRegion(region)
  }, [isValidated])

  const defineRegion = useCallback((atomIds: number[] | { x: number; y: number }) => {
    if (!selectedRegion) return

    // Handle both old mock selection and new real atom selection
    if (Array.isArray(atomIds)) {
      // Real atom selection from Ketcher
      setRegions(prev => ({
        ...prev,
        [selectedRegion]: `Selected atoms: ${atomIds.join(', ')}`
      }))
      
      // Store atom IDs for backend processing
      if (ketcherInstance) {
        // You can use this data to send to backend for fragment analysis
        console.log(`Region ${selectedRegion} atoms:`, atomIds)
        
        // Example: Get more details about the selected atoms
        try {
          const struct = ketcherInstance.editor.struct()
          const atomLabels = atomIds.map((atomId: number) => {
            const atom = struct.atoms.get(atomId)
            return atom ? `${atomId}:${atom.label}` : `${atomId}:?`
          })
          console.log(`Region ${selectedRegion} atom details:`, atomLabels.join(', '))
        } catch (e) {
          console.warn('Could not get atom details:', e)
        }
      }
    } else {
      // Legacy mock selection (for backward compatibility)
      const selection = {
        x: atomIds.x - 30,
        y: atomIds.y - 30,
        width: 60,
        height: 60,
        region: selectedRegion
      }

      setMockSelections(prev => [...prev.filter(s => s.region !== selectedRegion), selection])
      setRegions(prev => ({
        ...prev,
        [selectedRegion]: `Selected atoms near (${Math.round(atomIds.x)}, ${Math.round(atomIds.y)})`
      }))
    }
    
    setSelectedRegion(null)
  }, [selectedRegion, ketcherInstance])

  const removeRegion = useCallback((region: 'core' | 'r1') => {
    setRegions(prev => ({ ...prev, [region]: null }))
    setMockSelections(prev => prev.filter(s => s.region !== region))
  }, [])

  const updateProperty = useCallback((property: keyof MoleculeEditorState['properties'], value: number) => {
    setProperties(prev => ({
      ...prev,
      [property]: {
        ...prev[property],
        value
      }
    }))
  }, [])

  const toggleProperty = useCallback((property: keyof MoleculeEditorState['properties']) => {
    setProperties(prev => ({
      ...prev,
      [property]: {
        ...prev[property],
        enabled: !prev[property].enabled
      }
    }))
  }, [])

  const isComplete = useCallback(() => {
    return isValidated
  }, [isValidated])

  const getProgress = useCallback(() => {
    const steps = [false]
    if (isValidated) steps[0] = true
    return steps
  }, [isValidated])

  const reset = useCallback(() => {
    setSmiles('')
    setIsValidated(false)
    setSelectedRegion(null)
    setRegions({ core: null, r1: null })
    setMockSelections([])
    // Don't reset ketcherReady - keep the editor initialized
  }, [])

  const handleKetcherInit = useCallback((ketcher: any) => {
    console.log('Ketcher init callback received in hook')
    
    // Store the instance
    setKetcherInstance(ketcher)
    
    // Store the ketcher instance globally for debugging
    // @ts-ignore
    window.ketcherInstance = ketcher
  }, [])

  const getSelectionInfo = useCallback(() => {
    if (!ketcherInstance || !ketcherInstance.editor) {
      console.warn('Ketcher editor not available')
      return null
    }

    try {
      // Get the current selection
      const selection = ketcherInstance.editor.selection()
      
      if (!selection || (!selection.atoms || selection.atoms.length === 0)) {
        console.log('No atoms selected')
        return { atoms: [], bonds: [], summary: { atomCount: 0, bondCount: 0 } }
      }

      // Get the molecular structure
      const struct = ketcherInstance.editor.struct()
      
      const selectionInfo = {
        atoms: [] as any[],
        bonds: [] as any[],
        summary: {
          atomCount: 0,
          bondCount: 0,
          hasOtherEntities: false
        }
      }

      // Process selected atoms
      if (selection.atoms && selection.atoms.length > 0) {
        console.log(`Found ${selection.atoms.length} selected atoms:`)
        
        selection.atoms.forEach((atomId: number) => {
          const atom = struct.atoms.get(atomId)
          if (atom) {
            const atomInfo = {
              id: atomId,
              label: atom.label,
              position: {
                x: atom.pp?.x || 0,
                y: atom.pp?.y || 0
              },
              fragment: atom.fragment,
              charge: atom.charge || 0,
              isotope: atom.isotope || null,
              radical: atom.radical || 0,
              implicitH: atom.implicitH,
              neighbors: [] as any[]
            }

            // Get neighbor information
            try {
              const neighbors = struct.atomGetNeighbors(atomId)
              atomInfo.neighbors = neighbors.map((nei: any) => ({
                atomId: nei.aid,
                bondId: nei.bid
              }))
            } catch (e) {
              console.warn(`Could not get neighbors for atom ${atomId}`)
            }

            selectionInfo.atoms.push(atomInfo)
            console.log(`  Atom ${atomId}: ${atom.label} at (${atomInfo.position.x.toFixed(2)}, ${atomInfo.position.y.toFixed(2)})`)
          }
        })
      }

      // Process selected bonds
      if (selection.bonds && selection.bonds.length > 0) {
        console.log(`Found ${selection.bonds.length} selected bonds:`)
        
        selection.bonds.forEach((bondId: number) => {
          const bond = struct.bonds.get(bondId)
          if (bond) {
            const bondInfo = {
              id: bondId,
              begin: bond.begin,
              end: bond.end,
              type: bond.type,
              stereo: bond.stereo || 0,
              topology: bond.topology || 0,
              reactingCenterStatus: bond.reactingCenterStatus || 0
            }
            selectionInfo.bonds.push(bondInfo)
            console.log(`  Bond ${bondId}: connects atoms ${bond.begin} and ${bond.end}`)
          }
        })
      }

      // Update summary
      selectionInfo.summary.atomCount = selectionInfo.atoms.length
      selectionInfo.summary.bondCount = selectionInfo.bonds.length
      selectionInfo.summary.hasOtherEntities = !!(
        selection.texts?.length || 
        selection.rxnArrows?.length || 
        selection.rxnPluses?.length || 
        selection.simpleObjects?.length
      )

      console.log('Selection Summary:', selectionInfo.summary)
      console.log('Full Selection Data:', selectionInfo)

      return selectionInfo
    } catch (error) {
      console.error('Error getting selection info:', error)
      return null
    }
  }, [ketcherInstance])

  const getSelectionAsSmiles = useCallback(async () => {
    if (!ketcherInstance || !ketcherInstance.editor) {
      console.warn('Ketcher editor not available')
      return null
    }

    try {
      const selection = ketcherInstance.editor.selection()
      
      if (!selection || !selection.atoms || selection.atoms.length === 0) {
        console.log('No atoms selected')
        return null
      }

      // Get the full molecule SMILES
      const fullMoleculeSmiles = await ketcherInstance.getSmiles()
      
      // Get the molecular structure
      const struct = ketcherInstance.editor.struct()
      
      // Create a set of selected atoms for quick lookup
      const selectedAtomSet = new Set(selection.atoms)
      
      // Find all bonds that connect selected atoms
      const selectedBonds = [] as number[]
      const attachmentPoints = [] as any[]
      
      struct.bonds.forEach((bond: any, bondId: number) => {
        const beginSelected = selectedAtomSet.has(bond.begin)
        const endSelected = selectedAtomSet.has(bond.end)
        
        if (beginSelected && endSelected) {
          // Both atoms are selected - this bond is part of the fragment
          selectedBonds.push(bondId)
        } else if (beginSelected || endSelected) {
          // One atom is selected, one is not - this is an attachment point
          attachmentPoints.push({
            atomIndex: beginSelected ? bond.begin : bond.end,
            externalAtomIndex: beginSelected ? bond.end : bond.begin
          })
        }
      })

      // Create atom map with essential information
      const atomMap = selection.atoms.map((atomId: number) => {
        const atom = struct.atoms.get(atomId)
        const neighbors = struct.atomGetNeighbors(atomId)
        return {
          index: atomId,
          element: atom?.label || 'C',
          neighbors: neighbors.map((n: any) => n.aid)
        }
      })

      // Create bond map
      const bondMap = selectedBonds.map((bondId: number) => {
        const bond = struct.bonds.get(bondId)
        return {
          index: bondId,
          atoms: [bond?.begin, bond?.end],
          type: bond?.type || 1 // 1=single, 2=double, 3=triple, 4=aromatic
        }
      })

      // Return clean JSON structure
      return {
        fullMoleculeSmiles,
        selectedAtomIndices: selection.atoms,
        selectedBondIndices: selectedBonds,
        atomMap,
        bondMap,
        attachmentPoints
      }
    } catch (error) {
      console.error('Error converting selection to SMILES:', error)
      return null
    }
  }, [ketcherInstance])

  return {
    // State
    smiles,
    isValidated,
    isValidating,
    selectedRegion,
    regions,
    properties,
    mockSelections,
    ketcherInstance,
    debugMode,
    
    // Actions
    setSmiles,
    validateMolecule,
    selectRegion,
    defineRegion,
    removeRegion,
    updateProperty,
    toggleProperty,
    handleKetcherInit,
    getSelectionInfo,
    getSelectionAsSmiles,
    
    // Computed
    isComplete,
    getProgress,
    reset
  }
}