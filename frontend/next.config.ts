import type { NextConfig } from "next";

const nextConfig: NextConfig = {
  reactStrictMode: false, // Temporarily disabled to prevent Ketcher double-initialization issues
  webpack: (config, { isServer }) => {
    // Handle Ketcher's node dependencies
    if (!isServer) {
      config.resolve.fallback = {
        ...config.resolve.fallback,
        fs: false,
        path: false,
        crypto: false,
        buffer: false,
        stream: false,
      };
    }

    return config;
  },
  // Transpile ketcher packages
  transpilePackages: ["ketcher-core", "ketcher-react", "ketcher-standalone"],
};

export default nextConfig;
