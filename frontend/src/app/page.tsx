import { ALOAppContainer } from '../features/canvas/molecule-gui/alo-app-container'
import { SidebarOverlay } from '@/components/layout/sidebar-overlay'
import 'ketcher-react/dist/index.css'
import '@/styles/ketcher.css'

export default function HomePage() {
  return (
    <>
      <SidebarOverlay />
      <div className="min-h-screen w-full bg-gray-50 py-8 px-4">
        <div className="max-w-7xl mx-auto space-y-8">
          {/* Header */}
          <div className="text-center">
            <h1 className="text-4xl font-bold text-gray-900 mb-2">Welcome to ALO</h1>
            <p className="text-lg text-gray-600">AI-based Ligand Optimization</p>
          </div>
          
          {/* Molecule Editor Container */}
          <div className="bg-white border border-border overflow-hidden" style={{ maxHeight: '600px' }}>
            <ALOAppContainer />
          </div>
        </div>
      </div>
    </>
  )
}