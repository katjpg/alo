"use client"

import { useState } from 'react'
import { MolecularInput } from './molecule-input'

export function ChatInput() {
  const [inputValue, setInputValue] = useState('')
  const [isLoading, setIsLoading] = useState(false)

  const handleFileUpload = (file: File) => {
    // TODO: File upload functionality will be implemented later
    console.log('File uploaded:', file.name)
    setInputValue('') // Clear text input when file is uploaded
  }

  const handleSend = () => {
    if (!inputValue.trim()) return
    
    setIsLoading(true)
    // TODO: Send functionality will be implemented later
    console.log('Sending:', inputValue)
    
    // Simulate loading
    setTimeout(() => {
      setInputValue('')
      setIsLoading(false)
    }, 1000)
  }

  return (
    <MolecularInput
      value={inputValue}
      onChange={setInputValue}
      onSend={handleSend}
      onFileUpload={handleFileUpload}
      disabled={isLoading}
      mode="chat"
    />
  )
}