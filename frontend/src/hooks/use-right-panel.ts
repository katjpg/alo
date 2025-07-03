"use client"

import { create } from 'zustand'
import { persist } from 'zustand/middleware'

interface RightPanelState {
  isOpen: boolean
  toggle: () => void
  setIsOpen: (open: boolean) => void
}

export const useRightPanelState = create<RightPanelState>()(
  persist(
    (set) => ({
      isOpen: true,
      toggle: () => set((state) => ({ isOpen: !state.isOpen })),
      setIsOpen: (open: boolean) => set({ isOpen: open }),
    }),
    {
      name: 'alo-right-panel-state',
    }
  )
)