import { SidebarOverlay } from "@/components/layout/sidebar-overlay"

interface ChatLayoutProps {
  children: React.ReactNode
}

export default async function ChatLayout({ children }: ChatLayoutProps) {
  return (
    <>
      <SidebarOverlay />
      {children}
    </>
  )
}