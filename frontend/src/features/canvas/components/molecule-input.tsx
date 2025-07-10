"use client"

import { useState } from 'react'
import { IconPaperclip, IconSend, IconFileText } from '@tabler/icons-react'
import { Button } from '@/components/ui/button'
import { Textarea } from '@/components/ui/textarea'
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from '@/components/ui/tooltip'
import { cn } from '@/lib/utils'
import { PropertySelector } from './property-selector'
import { PropertyBadges } from './property-badges'

interface MolecularInputProps {
  value: string
  onChange: (value: string) => void
  onSend: () => void
  onFileUpload: (file: File) => void
  disabled?: boolean
  mode?: 'initial' | 'chat'
}

export function MolecularInput({
  value,
  onChange,
  onSend,
  onFileUpload,
  disabled = false,
  mode = 'initial'
}: MolecularInputProps) {
  const [isDragOver, setIsDragOver] = useState(false)
  const [uploadedFile, setUploadedFile] = useState<File | null>(null)
  const [selectedProperties, setSelectedProperties] = useState<string[]>([])

  const handlePropertyChange = (propertyId: string, checked: boolean) => {
    if (checked) {
      setSelectedProperties(prev => [...prev, propertyId])
    } else {
      setSelectedProperties(prev => prev.filter(id => id !== propertyId))
    }
  }

  const handleRemoveProperty = (propertyId: string) => {
    setSelectedProperties(prev => prev.filter(id => id !== propertyId))
  }

  const handleFileInput = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0]
    if (file) {
      setUploadedFile(file)
      onFileUpload(file)
    }
  }

  const handleDragOver = (e: React.DragEvent) => {
    e.preventDefault()
    setIsDragOver(true)
  }

  const handleDragLeave = (e: React.DragEvent) => {
    e.preventDefault()
    setIsDragOver(false)
  }

  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault()
    setIsDragOver(false)
    
    const file = e.dataTransfer.files[0]
    if (file && (file.name.endsWith('.mol') || file.name.endsWith('.sdf'))) {
      setUploadedFile(file)
      onFileUpload(file)
    }
  }

  const handleKeyPress = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault()
      onSend()
    }
  }

  const hasContent = value.trim() || uploadedFile || selectedProperties.length > 0

  return (
    <div className="w-full max-w-2xl">
      <div
        className={cn(
          "rounded-xl bg-background transition-all",
          isDragOver
            ? "border border-primary bg-primary/5"
            : "border border-border",
          disabled && "opacity-50 pointer-events-none"
        )}
      >
        {/* SMILES Input Section */}
        <div
          className="relative"
          onDragOver={handleDragOver}
          onDragLeave={handleDragLeave}
          onDrop={handleDrop}
        >
          <div className="flex items-start">
            <Textarea
              value={value}
              onChange={(e: React.ChangeEvent<HTMLTextAreaElement>) => onChange(e.target.value)}
              onKeyDown={handleKeyPress}
              placeholder={uploadedFile ? "File uploaded successfully!" : mode === 'chat' ? "Ask about molecules, properties, or optimizations..." : "Enter your initial ligand and select your desired properties..."}
              disabled={disabled || !!uploadedFile}
              className={cn(
                "flex-1 border-0 bg-transparent px-4 focus-visible:ring-0 focus-visible:ring-offset-0 resize-none shadow-none outline-none",
                mode === 'chat' 
                  ? "text-sm min-h-[40px] placeholder:text-sm pt-2 pb-2" 
                  : "text-lg min-h-[80px] placeholder:text-base pt-4 pb-6"
              )}
            />
          </div>
        </div>

        {/* Property Selection and Start Discovery Section */}
        <div className="p-4 space-y-3">
          <div className="flex items-center justify-between gap-4">
            {/* Property selector only in initial mode */}
            {mode === 'initial' && (
              <PropertySelector
                selectedProperties={selectedProperties}
                onPropertyChange={handlePropertyChange}
              />
            )}
            
            <div className={cn(
              "flex items-center gap-2",
              mode === 'chat' && "w-full justify-end"
            )}>
              {/* File Upload - only in initial mode */}
              {mode === 'initial' && (
                <>
                  {/* File Upload Indicator */}
                  {uploadedFile && (
                    <div className="flex items-center gap-2 px-3 py-2 border border-border rounded-lg bg-background">
                      <IconFileText className="h-4 w-4 text-muted-foreground" />
                      <span className="text-sm text-muted-foreground">{uploadedFile.name}</span>
                      <button
                        onClick={() => {
                          setUploadedFile(null)
                          onChange('')
                        }}
                        className="text-foreground hover:text-foreground/80 ml-1"
                      >
                        ×
                      </button>
                    </div>
                  )}
                  
                  <input
                    type="file"
                    accept=".mol,.sdf"
                    onChange={handleFileInput}
                    className="hidden"
                    id="file-upload"
                  />
                  <TooltipProvider>
                    <Tooltip>
                      <TooltipTrigger asChild>
                        <label htmlFor="file-upload">
                          <Button
                            variant="outline"
                            size="sm"
                            asChild
                            className="h-8 w-8 hover:bg-muted cursor-pointer"
                          >
                            <span>
                              <IconPaperclip className="h-4 w-4" />
                            </span>
                          </Button>
                        </label>
                      </TooltipTrigger>
                      <TooltipContent>
                        <p>Upload file</p>
                      </TooltipContent>
                    </Tooltip>
                  </TooltipProvider>
                </>
              )}
              
              {/* Send Button - different styles for each mode */}
              {mode === 'initial' ? (
                <Button
                  onClick={onSend}
                  disabled={!hasContent || disabled}
                  className="bg-blue-600 hover:bg-blue-700 text-white gap-2"
                  size="sm"
                >
                  Start Discovery
                  <IconSend className="h-4 w-4" />
                </Button>
              ) : (
                <Button
                  onClick={onSend}
                  disabled={!value.trim() || disabled}
                  className="bg-blue-600 hover:bg-blue-700 text-white h-8 w-8"
                  size="sm"
                >
                  <IconSend className="h-4 w-4" />
                </Button>
              )}
            </div>
          </div>
          
          {/* Property badges only in initial mode */}
          {mode === 'initial' && (
            <PropertyBadges
              selectedProperties={selectedProperties}
              onRemoveProperty={handleRemoveProperty}
            />
          )}
        </div>
      </div>

      {/* File Type Helper - only in initial mode */}
      {mode === 'initial' && (
        <div className="mt-2 text-center text-xs text-muted-foreground">
          Supported formats: SMILES strings, .mol, .sdf files
        </div>
      )}
    </div>
  )
}