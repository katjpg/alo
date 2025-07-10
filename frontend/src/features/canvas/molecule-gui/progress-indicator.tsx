"use client"

import { cn } from '@/lib/utils'

interface ProgressIndicatorProps {
  steps: boolean[]
}

export function ProgressIndicator({ steps }: ProgressIndicatorProps) {
  return (
    <div className="flex justify-center gap-2 py-4">
      {steps.map((completed, index) => (
        <div
          key={index}
          className={cn(
            "w-2 h-2 rounded-full transition-all duration-300",
            completed 
              ? "bg-green-500 scale-125" 
              : "bg-gray-300"
          )}
        />
      ))}
    </div>
  )
}