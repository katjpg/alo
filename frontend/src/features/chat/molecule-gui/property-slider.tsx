"use client"

import { Slider } from '@/components/ui/slider'
import { cn } from '@/lib/utils'

interface PropertySliderProps {
  label: string
  value: number
  min: number
  max: number
  step: number
  unit?: string
  tooltip?: string
  onChange: (value: number) => void
  rangeLabels?: [string, string]
}

export function PropertySlider({
  label,
  value,
  min,
  max,
  step,
  unit = '',
  tooltip,
  onChange,
  rangeLabels
}: PropertySliderProps) {
  return (
    <div className="space-y-2">
      <div className="flex items-center justify-between">
        {label && (
          <label className="text-sm font-medium text-gray-700">
            {label}
            {tooltip && (
              <span className="ml-1 text-xs text-gray-500">ⓘ</span>
            )}
          </label>
        )}
        <span className={cn(
          "text-sm font-semibold text-gray-900 min-w-[45px] text-right",
          !label && "ml-auto"
        )}>
          {value.toFixed(step < 1 ? 1 : 0)}{unit}
        </span>
      </div>
      
      <Slider
        value={[value]}
        onValueChange={([newValue]) => onChange(newValue)}
        min={min}
        max={max}
        step={step}
        className="w-full"
      />
      
      <div className="flex justify-between text-xs text-gray-500">
        <span>{rangeLabels ? rangeLabels[0] : min}</span>
        <span>{rangeLabels ? rangeLabels[1] : max}</span>
      </div>
    </div>
  )
}