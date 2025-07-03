import { cookies } from "next/headers"
import { ChatPanels } from "@/features/chat/components/chat-panels"
import { SidebarOverlay } from "@/components/layout/sidebar-overlay"

interface ChatLayoutProps {
  children: React.ReactNode
}

export default async function ChatLayout({ children }: ChatLayoutProps) {
  const cookieStore = await cookies()
  const layout = cookieStore.get("react-resizable-panels:layout:chat")
  const defaultLayout = layout ? JSON.parse(layout.value) : [26, 52, 22]

  return (
    <>
      <SidebarOverlay />
      <ChatPanels defaultLayout={defaultLayout}>
        {children}
      </ChatPanels>
    </>
  )
}