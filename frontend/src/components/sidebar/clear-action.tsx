'use client'

import { IconTrash } from '@tabler/icons-react'
import { Button } from '@/components/ui/button'
import {
  AlertDialog,
  AlertDialogAction,
  AlertDialogCancel,
  AlertDialogContent,
  AlertDialogDescription,
  AlertDialogFooter,
  AlertDialogHeader,
  AlertDialogTitle,
  AlertDialogTrigger,
} from '@/components/ui/alert-dialog'

interface ClearHistoryActionProps {
  empty: boolean
}

export function ClearAction({ empty }: ClearHistoryActionProps) {
  const handleClearHistory = async () => {
    try {
      // Implement clear history logic here
      console.log('Clearing discovery history...')
      
      // Dispatch event to refresh history
      window.dispatchEvent(new CustomEvent('discovery-history-updated'))
    } catch (error) {
      console.error('Failed to clear history:', error)
    }
  }

  if (empty) {
    return null
  }

  return (
    <AlertDialog>
      <AlertDialogTrigger asChild>
        <Button 
          variant="ghost" 
          size="sm" 
          className="h-6 w-6 p-0 hover:bg-destructive/10"
        >
          <IconTrash className="size-3" />
        </Button>
      </AlertDialogTrigger>
      <AlertDialogContent>
        <AlertDialogHeader>
          <AlertDialogTitle>Clear Discovery History</AlertDialogTitle>
          <AlertDialogDescription>
            Are you sure you want to clear all discovery sessions? This action cannot be undone.
          </AlertDialogDescription>
        </AlertDialogHeader>
        <AlertDialogFooter>
          <AlertDialogCancel>Cancel</AlertDialogCancel>
          <AlertDialogAction 
            onClick={handleClearHistory}
            className="bg-destructive text-destructive-foreground hover:bg-destructive/90"
          >
            Clear History
          </AlertDialogAction>
        </AlertDialogFooter>
      </AlertDialogContent>
    </AlertDialog>
  )
}