# ALO Interface

> ALO (AI-based Ligand Optimization) is a Human-AI collaborative tool designed to support drug discovery research. It combines fragment-based drug design with AI-driven molecular generation to enable iterative refinement and systematic exploration of molecular optimization pathways.

## Table of Contents
- [ALO Interface](#alo-interface)
  - [Table of Contents](#table-of-contents)
  - [Intro](#intro)
  - [Packages Used](#packages-used)
  - [Getting Started](#getting-started)
  - [Project Structure](#project-structure)
  - [Navigation](#navigation)
    - [Root (`/`)](#root-)
    - [Chat (`/chat`)](#chat-chat)
  - [Features](#features)
    - [`/chat` System](#chat-system)
      - [Primitive Components](#primitive-components)
      - [Thought Tree](#thought-tree)
    - [State Management](#state-management)

## Intro

The ALO interface is a Next.js 15 application that includes visual thought trees and iterative chat interactions. Features a three-panel layout system.. 

## Packages Used

| Package | Version | Purpose |
|---------|---------|---------|
| Next.js | 15.3.4 | React framework with App Router |
| React | 19.0.0 | UI library |
| TypeScript | 5.x | Type safety |
| @xyflow/react | 12.8.1 | Thought tree visualization |
| shadcn/ui | Latest | Component library |
| zustand | 5.0.5 | State management |
| react-resizable-panels | 2.2.0 | Resizable layout panels |
| Tailwind CSS | 4.x | Styling framework |

## Getting Started

```bash
# Install dependencies
pnpm install

# Run development server
pnpm run dev


```

## Project Structure

```
src/
├── app/                           # Next.js App Router
│   ├── api/                       # API routes
│   ├── chat/                      # Chat interface pages
│   │   └── [id]/                  # Dynamic chat sessions
│   ├── layout.tsx                 # Root layout
│   └── page.tsx                   # Landing page
│
├── components/                    # Reusable UI components
│   ├── layout/                    # Layout components
│   ├── sidebar/                   # Navigation components
│   └── ui/                        # Base UI components
│
├── features/                      # Feature-specific modules
│   ├── artifacts/                 # Saved results (TODO)
│   ├── chat/                      # Chat system
│   │   ├── components/            # Chat UI elements
│   │   ├── info-panel/            # Details panel (TODO)
│   │   └── thought-tree/          # Visualization
│   │       ├── edges/             # Edge types
│   │       ├── nodes/             # Node components
│   │       └── tree-hierarchy.tsx # Main tree
│   └── threads/                   # Conversation history (TODO)
│
├── hooks/                         # Custom React hooks
├── lib/                           # Utility functions
└── styles/                        # Global styles
```

## Navigation

### Root (`/`)
```tsx
<div className="flex h-screen">
  {/* Collapsible Sidebar - 60px collapsed, 240px expanded */}
  <AppSidebar />
  
  {/* Main Content Area */}
  <main className="flex-1 p-8">
    {/* Welcome Message */}
    <h1>Welcome to ALO</h1>
    
    {/* Large Molecule Input Form */}
    <MoleculeInput mode="initial" />
    
    {/* Start Discovery Button */}
  </main>
</div>
```

### Chat (`/chat`)
```tsx
<ChatPanels>
  {/* Left Panel - 26% width (resizable) */}
  <div>
    <ChatMessages />
    <MoleculeInput mode="chat" />
  </div>
  
  {/* Middle Panel - 52% width (resizable) */}
  <TreeHierarchy />
  
  {/* Right Panel - 280px fixed (toggleable) */}
  <InfoPanel />
</ChatPanels>
```

## Features

### `/chat` System

Provides an integrated environment for molecular discovery through AI collaboration.

#### Primitive Components

1. **MoleculeInput** (`molecule-input.tsx`)
   - SMILES string input with validation
   - File upload (.mol, .sdf) with drag-and-drop
   - Property selection badges
   - Responsive design with mode-specific styling

2. **ChatInput** (`chat-input.tsx`)
   - Textarea with auto-resize
   - Send button with loading states
   - Keyboard shortcuts (Cmd/Ctrl + Enter)
   - Placeholder guidance text

3. **ChatPanels** (`chat-panels.tsx`)
   - Three-panel resizable layout using react-resizable-panels
   - Persistent panel sizes via cookies
   - Mobile-responsive tab navigation

#### Thought Tree

The thought tree (`thought-tree/`) visualizes the molecular optimization workflow using ReactFlow:

**Node Types:**
- **Initial Molecule**: Starting compound with properties
- **Generated Molecules**: AI-suggested variations
- **Edited Molecules**: User-modified structures
- **Cluster Nodes**: Grouped similar molecules

**Edge Types:**
- **Decomposition**: Breaking down molecules (animated dashed lines)
- **Generation**: Creating new variants
- **Editing**: User modifications

**Features:**
- Interactive zoom and pan
- MiniMap for navigation
- Custom node designs with molecular previews
- Hierarchical layout algorithm
- Real-time updates during optimization

### State Management

The application uses Zustand for persistent state management across sidebar, panel sizes, and session data.