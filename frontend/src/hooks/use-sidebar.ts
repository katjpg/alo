"use client"

import { create } from 'zustand'
import { persist } from 'zustand/middleware'

interface SidebarState {
  isOpen: boolean
  isHydrated: boolean
  toggle: () => void
  setIsOpen: (open: boolean) => void
  setHydrated: (hydrated: boolean) => void
}

export const useSidebarState = create<SidebarState>()(
  persist(
    (set) => ({
      isOpen: true,
      isHydrated: false,
      toggle: () => set((state) => ({ isOpen: !state.isOpen })),
      setIsOpen: (open: boolean) => set({ isOpen: open }),
      setHydrated: (hydrated: boolean) => set({ isHydrated: hydrated }),
    }),
    {
      name: 'alo-sidebar-state',
      partialize: (state) => ({ isOpen: state.isOpen }),
    }
  )
)