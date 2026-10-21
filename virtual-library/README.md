# Biblioteca virtual

Una estantería personal en 3D construida a partir de la
[Virtual Library Build Guide](https://github.com/carollia99/virtual-library-guide/blob/main/virtual_library_build_guide.md):

- Un riel horizontal infinito de libros reales en 3D, con la tapa verdadera y un lomo
  del color muestreado de esa tapa.
- Al pasar el mouse el libro se adelanta y aparece una tarjeta con título, autor, año,
  encuadernación, puntuación y géneros.
- Al hacer clic el libro sale de la estantería y gira hasta mostrar la tapa, con la
  estantería desenfocada detrás (← / → para pasar de libro, Esc para devolverlo).
- Búsqueda en lenguaje natural ("What are you looking for?") con IA, más filtros por género.
- Botón opcional **Recommend a book** para que las visitas te recomienden libros, que
  aparecen en una segunda estantería ("Recommended to me").

Stack: TanStack Start + React 19 + Vite + Tailwind CSS v4 + TypeScript.

## Puesta en marcha

Requiere Node 22.18 o superior.

```bash
cd virtual-library
npm install
cp .env.example .env   # opcional: claves de IA y Supabase
npm run preview        # compila y abre el sitio en http://localhost:3000
```

`npm run preview` es lo que conviene para **ver y usar** el sitio: compila la versión de
producción, que es la rápida. `npm run dev` sirve para editar el código (se recarga solo
con cada cambio), pero es bastante más lento, sobre todo al cargar.

## 1. Cargar tu biblioteca de Goodreads

La estantería empieza vacía: el sitio no inventa libros.

1. En Goodreads: **My Books → Import and export → Export Library**. Esperá el mail y
   descargá `goodreads_library_export.csv`.
2. Copialo dentro de `virtual-library/` (está en `.gitignore`, no se sube al repo).
3. Corré (poné tu email después de `--contact`):

   ```bash
   npm run import:goodreads -- goodreads_library_export.csv --contact tu@email.com
   ```

Open Library permite 1 consulta por segundo a quien no se identifica, y 3 por segundo
si la consulta lleva un email de contacto. Por eso `--contact` hace la importación
unas 3 veces más rápida: unos 2–3 minutos para 300 libros, contra unos 7 sin el email.
El email solo viaja en el encabezado de las consultas a Open Library. Al empezar, el
script te dice cuánto va a tardar, y mientras corre muestra cuánto falta.

Lo que ya buscó queda guardado en `.cache/` (no se sube al repo). Si lo cortás con
Ctrl+C, o si algún libro falla porque Open Library no respondió, corré el mismo comando
otra vez: retoma desde donde quedó y completa solo lo que falta. Volver a importar
después de agregar libros en Goodreads tarda segundos. Con `--refresh` se vuelve a
buscar todo de cero.

El script (`scripts/import-goodreads.ts`) reescribe `src/data/books.ts`:

- Solo toma filas con *Exclusive Shelf* `read` o `currently-reading`, de la más reciente
  a la más vieja.
- Busca cada libro en Open Library (por ISBN y por título/autor) para obtener la tapa,
  un resumen corto y los temas.
- Muestrea el borde izquierdo de la tapa para el color del lomo (`spine`), elige un
  acento saturado (`band`) y el color de la tinta (`ink`).
- Deriva encuadernación, tamaño, inclinación, profundidad y desgaste a partir de la
  cantidad de páginas y de un hash de la clave de Open Library (siempre da lo mismo).
- Asigna géneros solo de este conjunto: Fiction, Nonfiction, Sci-Fi,
  Mystery & Thriller, Fantasy, Romance. Usa tus estantes de Goodreads y los temas de
  Open Library, en inglés o en castellano. Los botones de filtro se generan solos.
- Quita del título la serie y la edición ("The Name of the Wind (The Kingkiller
  Chronicle, #1)" → "The Name of the Wind"; "Sapiens (Spanish Edition)" → "Sapiens").

Si preferís no usar el script, también podés pasarle el CSV a un asistente de IA junto
con el prompt de la sección 4 de la guía; el resultado es el mismo archivo.

## 2. Búsqueda con IA (opcional)

La búsqueda corre en el servidor (`src/lib/librarySearch.functions.ts` /
`librarySearch.server.ts`), así que la clave nunca llega al navegador. En `.env`:

| Variable | Uso |
|---|---|
| `ANTHROPIC_API_KEY` | Usa Claude (por defecto `claude-haiku-5-5`, cambiable con `ANTHROPIC_MODEL`). Tiene prioridad. |
| `OPENAI_API_KEY` (+ `OPENAI_BASE_URL`, `OPENAI_MODEL`) | Cualquier API compatible con OpenAI: OpenAI, Gemini, Groq, OpenRouter, Lovable. |

Sin ninguna clave la búsqueda sigue funcionando, con coincidencia de palabras clave en
el navegador (misma interfaz, ranking más simple).

## 3. Estantería de recomendaciones (opcional)

1. Creá un proyecto en [Supabase](https://supabase.com) (o cualquier Postgres).
2. Ejecutá `supabase/migrations/0001_recommendations.sql` en el editor SQL. Los `GRANT`
   y las políticas de RLS permiten leer e insertar con límites de longitud, pero no
   editar ni borrar.
3. Poné `VITE_SUPABASE_URL` y `VITE_SUPABASE_ANON_KEY` en `.env`.

Sin esas variables el botón **Recommend a book** y la segunda estantería no se muestran.

## Personalizar

- **Velocidad de las animaciones**: todos los tiempos están juntos en
  `src/components/motion.ts`, en milisegundos. Por ejemplo `HOVER_MS` (el libro que se
  adelanta al pasar el mouse) o `PULL_MS` (el libro que sale de la estantería). Un
  número más chico hace la animación más rápida.

- **Título**: el texto en `src/components/TypedTitle.tsx` y el `head()` de
  `src/routes/index.tsx` (título y metadatos `og:`/`twitter:`).
- **Sensación del hover**: `HOVER_PULL` (96) y `HOVER_LIFT` (-26) en `BookSpine.tsx`.
- **Paleta**: los tokens oklch de `:root` en `src/styles.css`; los componentes no
  tienen colores fijos.

## Publicar

Es un build estándar de Vite + Nitro: `npm run build` genera un servidor Node en
`.output/` que corre en Railway, Render, Fly, un VPS o Docker. Para Vercel, Netlify o
Cloudflare, ver la sección de deployment de TanStack Start (se cambia el preset o el
plugin en `vite.config.ts`). Antes de publicar, revisá que el título y los metadatos
`og:`/`twitter:` sean los tuyos.

## Estructura

```
src/
  data/books.ts                 tipo Book + el array `books` (generado)
  components/
    bookFaces.ts                COVER_W y tipografías del lomo
    TypedTitle.tsx              "Welcome to my library" tipeado
    BookSpine.tsx               el libro 3D (lomo, tapas, páginas) + tarjeta flotante
    Shelf.tsx                   riel infinito, curvatura, rueda/arrastre/teclado
    BookDetail.tsx              el libro fuera de la estantería
    LibraryFilter.tsx           búsqueda con IA + géneros
    RecommendBookDialog.tsx     modal para recomendar
  hooks/useRecommendations.ts   lectura/inserción en Supabase
  lib/
    bookPhysics.ts              propiedades físicas y paleta (compartido con el script)
    openLibrary.ts              búsqueda en Open Library, paleta de tapas, buildBook
    keywordSearch.ts            búsqueda sin IA
    librarySearch.*.ts          búsqueda con IA en el servidor
  routes/__root.tsx, index.tsx  documento y página principal
scripts/import-goodreads.ts     CSV de Goodreads → src/data/books.ts
supabase/migrations/            SQL de la tabla de recomendaciones
```
