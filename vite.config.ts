import { defineConfig } from 'vite'
import { readFileSync } from 'node:fs'
import react from '@vitejs/plugin-react'
import tailwindcss from '@tailwindcss/vite'

// https://vite.dev/config/
export default defineConfig({
  base: './',
  plugins: [
    {
      name: 'app-source-template',
      transformIndexHtml: {
        order: 'pre',
        handler: () => readFileSync(new URL('./index.template.html', import.meta.url), 'utf8'),
      },
    },
    react(),
    tailwindcss(),
  ],
})
