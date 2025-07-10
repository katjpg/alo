"use client"

import { Button } from "@/components/ui/button"
import {
  DropdownMenu,
  DropdownMenuCheckboxItem,
  DropdownMenuContent,
  DropdownMenuLabel,
  DropdownMenuSeparator,
  DropdownMenuTrigger,
} from "@/components/ui/dropdown-menu"
import { IconChevronDown } from "@tabler/icons-react"

export interface ADMETProperty {
  id: string
  name: string
  fullName: string
  description: string
}

export const admetProperties: ADMETProperty[] = [
  {
    id: "BBBP",
    name: "BBBP",
    fullName: "Blood-Brain Barrier Permeability",
    description: "Predicts CNS drug penetration ability"
  },
  {
    id: "HIA",
    name: "HIA", 
    fullName: "Human Intestinal Absorption",
    description: "Estimates oral bioavailability potential"
  },
  {
    id: "Mutag",
    name: "Mutag",
    fullName: "Mutagenicity",
    description: "Screens for DNA damage and cancer risk"
  },
  {
    id: "DRD2",
    name: "DRD2",
    fullName: "Dopamine Receptor D2 Inhibition", 
    description: "Predicts antipsychotic activity potential"
  },
  {
    id: "plogP",
    name: "plogP",
    fullName: "Penalized LogP",
    description: "Optimized lipophilicity for drug-likeness"
  },
  {
    id: "QED",
    name: "QED",
    fullName: "Quantitative Estimation of Drug-likeness",
    description: "Overall pharmaceutical attractiveness score"
  }
]

interface PropertySelectorProps {
  selectedProperties: string[]
  onPropertyChange: (propertyId: string, checked: boolean) => void
}

export function PropertySelector({ selectedProperties, onPropertyChange }: PropertySelectorProps) {
  return (
    <DropdownMenu>
      <DropdownMenuTrigger asChild>
        <Button variant="outline" className="gap-2">
          Select Properties
          <IconChevronDown className="h-4 w-4" />
        </Button>
      </DropdownMenuTrigger>
      <DropdownMenuContent 
        className="w-80 h-60" 
        align="start" 
        side="bottom"
        sideOffset={4}
      >
        <DropdownMenuLabel>ADMET Properties</DropdownMenuLabel>
        <DropdownMenuSeparator />
        {admetProperties.map((property) => (
          <DropdownMenuCheckboxItem
            key={property.id}
            checked={selectedProperties.includes(property.id)}
            onCheckedChange={(checked) => onPropertyChange(property.id, checked)}
            onSelect={(e) => e.preventDefault()}
            className="flex flex-col items-start py-3"
          >
            <div className="font-medium text-sm">{property.name}</div>
            <div className="text-sm text-muted-foreground">
              {property.fullName}
            </div>
            <div className="text-sm text-muted-foreground leading-tight">
              {property.description}
            </div>
          </DropdownMenuCheckboxItem>
        ))}
      </DropdownMenuContent>
    </DropdownMenu>
  )
}
