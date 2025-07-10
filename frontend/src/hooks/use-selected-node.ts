"use client"

import { create } from 'zustand'
import type { ThoughtTreeNode } from '@/features/canvas/thought-tree/data/mock-data'

interface SelectedNodeState {
  selectedNode: ThoughtTreeNode | null
  setSelectedNode: (node: ThoughtTreeNode | null) => void
  clearSelection: () => void
}

export const useSelectedNode = create<SelectedNodeState>((set) => ({
  selectedNode: null,
  setSelectedNode: (node) => set({ selectedNode: node }),
  clearSelection: () => set({ selectedNode: null }),
}))