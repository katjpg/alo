import { ChatInput } from '../features/chat/components/chat-input'

export default function HomePage() {
  return (
    <div className="flex h-screen w-full items-center justify-center bg-background">
      <div className="w-full max-w-3xl px-8">
        <ChatInput />
      </div>
    </div>
  )
}