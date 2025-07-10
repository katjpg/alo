"use client"

import { Badge } from "@/components/ui/badge"
import { 
  IconX, 
  IconBrain, 
  IconStethoscope, 
  IconAlertTriangle, 
  IconMoodSmile, 
  IconDroplet, 
  IconShield 
} from "@tabler/icons-react"
import { admetProperties } from "./property-selector"

interface PropertyBadgesProps {
  selectedProperties: string[]
  onRemoveProperty: (propertyId: string) => void
}

const propertyIcons = {
  BBBP: IconBrain,        // Blood-Brain Barrier - brain icon
  HIA: IconStethoscope,   // Human Intestinal Absorption - medical icon
  Mutag: IconAlertTriangle, // Mutagenicity - warning icon
  DRD2: IconMoodSmile,    // Dopamine Receptor - mood/psychiatric icon
  plogP: IconDroplet,     // Lipophilicity - droplet icon
  QED: IconShield         // Drug-likeness - shield/quality icon
} as const

export function PropertyBadges({ selectedProperties, onRemoveProperty }: PropertyBadgesProps) {
  if (selectedProperties.length === 0) {
    return null
  }

  return (
    <div className="flex items-center gap-2 flex-wrap mt-3">
      {selectedProperties.map((propertyId) => {
        const property = admetProperties.find(p => p.id === propertyId)
        const IconComponent = propertyIcons[propertyId as keyof typeof propertyIcons]
        
        if (!property || !IconComponent) return null

        return (
          <Badge
            key={propertyId}
            variant="outline"
            className="rounded-full pr-1 gap-1.5 cursor-default"
          >
            <IconComponent className="h-4 w-4" />
            {property.name}
            <button
              onClick={() => onRemoveProperty(propertyId)}
              className="ml-1 hover:bg-black/10 dark:hover:bg-white/10 rounded-full p-0.5 transition-colors"
            >
              <IconX className="h-3 w-3" />
            </button>
          </Badge>
        )
      })}
    </div>
  )
}