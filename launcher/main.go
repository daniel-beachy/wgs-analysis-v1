// Genome Dashboard launcher.
//
// One small, dependency-free binary per platform (macOS universal, Windows x64). It embeds the
// built viewer, finds the analysis results (wgs-data/releases) by walking up from wherever the
// binary lives, serves both on 127.0.0.1, and opens the default browser. It exits on its own
// about a minute after the last browser tab closes (the page sends a heartbeat).
package main

import (
	"embed"
	"errors"
	"flag"
	"fmt"
	"io/fs"
	"log"
	"net"
	"net/http"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"sync/atomic"
	"time"
)

//go:embed all:app
var appFS embed.FS

const preferredPort = 8787

var version = "dev"

func main() {
	dataFlag := flag.String("data", "", "path to wgs-data/releases (default: search upward from this program)")
	port := flag.Int("port", preferredPort, "port to try first")
	noBrowser := flag.Bool("no-browser", false, "do not open a browser")
	idle := flag.Duration("idle", 75*time.Second, "exit this long after the last browser heartbeat (0 = never)")
	flag.Parse()
	setupLog()

	releases, err := findReleases(*dataFlag)
	if err != nil {
		log.Printf("warning: %v (the dashboard will explain how to fix this)", err)
	} else {
		log.Printf("data: %s", releases)
	}

	ln, url, existing := listen(*port)
	if existing != "" {
		log.Printf("already running at %s; opening it", existing)
		if !*noBrowser {
			openBrowser(existing)
		}
		return
	}

	var last atomic.Int64
	last.Store(time.Now().Unix())
	mux := http.NewServeMux()
	mux.HandleFunc("/api/ping", func(w http.ResponseWriter, _ *http.Request) { fmt.Fprint(w, "genome-dashboard ", version) })
	mux.HandleFunc("/api/heartbeat", func(w http.ResponseWriter, _ *http.Request) { last.Store(time.Now().Unix()); w.WriteHeader(204) })
	mux.HandleFunc("/api/bye", func(w http.ResponseWriter, _ *http.Request) {
		// Closing one tab (or reloading) sends "bye"; give the page a few seconds to re-announce itself.
		last.Store(time.Now().Add(-*idle).Add(15 * time.Second).Unix())
		w.WriteHeader(204)
	})
	if releases != "" {
		mux.Handle("/data/", http.StripPrefix("/data/", dataHandler(releases)))
	}
	sub, _ := fs.Sub(appFS, "app")
	mux.Handle("/", noStore(http.FileServerFS(sub)))

	srv := &http.Server{Handler: mux, ReadHeaderTimeout: 10 * time.Second}
	go func() {
		if err := srv.Serve(ln); err != nil && !errors.Is(err, http.ErrServerClosed) {
			log.Fatal(err)
		}
	}()
	log.Printf("serving %s", url)
	if !*noBrowser {
		openBrowser(url)
	}
	if *idle == 0 {
		select {}
	}
	for range time.Tick(5 * time.Second) {
		if time.Since(time.Unix(last.Load(), 0)) > *idle {
			log.Printf("no open dashboard tabs for %s; exiting", *idle)
			return
		}
	}
}

// listen binds 127.0.0.1 on the preferred port, or reports an already-running instance, or falls back to any free port.
func listen(port int) (net.Listener, string, string) {
	addr := fmt.Sprintf("127.0.0.1:%d", port)
	if ln, err := net.Listen("tcp", addr); err == nil {
		return ln, "http://" + addr + "/", ""
	}
	c := http.Client{Timeout: time.Second}
	if r, err := c.Get("http://" + addr + "/api/ping"); err == nil {
		r.Body.Close()
		if r.StatusCode == 200 {
			return nil, "", "http://" + addr + "/"
		}
	}
	ln, err := net.Listen("tcp", "127.0.0.1:0")
	if err != nil {
		log.Fatal(err)
	}
	return ln, fmt.Sprintf("http://%s/", ln.Addr()), ""
}

// findReleases looks for wgs-data/releases/index.json in the flag, $WGS_RELEASES, then every ancestor of the executable.
func findReleases(flagVal string) (string, error) {
	var cands []string
	for _, v := range []string{flagVal, os.Getenv("WGS_RELEASES")} {
		if v != "" {
			cands = append(cands, v, filepath.Join(v, "releases"), filepath.Join(v, "wgs-data", "releases"))
		}
	}
	if exe, err := os.Executable(); err == nil {
		if r, err := filepath.EvalSymlinks(exe); err == nil {
			exe = r
		}
		for dir := filepath.Dir(exe); ; dir = filepath.Dir(dir) {
			cands = append(cands, filepath.Join(dir, "wgs-data", "releases"))
			if filepath.Dir(dir) == dir {
				break
			}
		}
	}
	if wd, err := os.Getwd(); err == nil {
		cands = append(cands, filepath.Join(wd, "wgs-data", "releases"))
	}
	for _, c := range cands {
		if st, err := os.Stat(filepath.Join(c, "index.json")); err == nil && !st.IsDir() {
			abs, _ := filepath.Abs(c)
			return abs, nil
		}
	}
	return "", errors.New("could not find wgs-data/releases/index.json next to or above this program; use --data PATH")
}

// dataHandler serves release files read-only, with HTTP Range support (DuckDB reads Parquet in ranges).
func dataHandler(root string) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		name := filepath.Base(r.URL.Path)
		if strings.HasPrefix(name, "._") || strings.HasPrefix(name, ".") || strings.Contains(r.URL.Path, "..") {
			http.NotFound(w, r)
			return
		}
		p := filepath.Join(root, filepath.FromSlash(strings.TrimPrefix(r.URL.Path, "/")))
		f, err := os.Open(p)
		if err != nil {
			http.NotFound(w, r)
			return
		}
		defer f.Close()
		st, err := f.Stat()
		if err != nil || st.IsDir() {
			http.NotFound(w, r)
			return
		}
		if strings.HasSuffix(p, ".json") {
			w.Header().Set("Content-Type", "application/json")
			w.Header().Set("Cache-Control", "no-cache")
		} else if strings.HasSuffix(p, ".html") {
			// Tool reports (e.g. PharmCAT) opened in their own tab; content-addressed like other objects.
			w.Header().Set("Content-Type", "text/html; charset=utf-8")
			w.Header().Set("Cache-Control", "public, max-age=31536000, immutable")
		} else {
			w.Header().Set("Content-Type", "application/octet-stream")
			// Objects are content-addressed and never change.
			w.Header().Set("Cache-Control", "public, max-age=31536000, immutable")
		}
		http.ServeContent(w, r, name, st.ModTime(), f)
	})
}

func noStore(h http.Handler) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		if r.URL.Path == "/" || strings.HasSuffix(r.URL.Path, ".html") {
			w.Header().Set("Cache-Control", "no-cache")
		}
		h.ServeHTTP(w, r)
	})
}

func openBrowser(url string) {
	var cmd *exec.Cmd
	switch runtime.GOOS {
	case "darwin":
		cmd = exec.Command("open", url)
	case "windows":
		cmd = exec.Command("rundll32", "url.dll,FileProtocolHandler", url)
	default:
		cmd = exec.Command("xdg-open", url)
	}
	if err := cmd.Start(); err != nil {
		log.Printf("open %s in your browser (%v)", url, err)
	}
}

// setupLog writes to stderr and, since GUI launches have no console, to a log file in the temp dir.
func setupLog() {
	log.SetFlags(log.Ltime)
	if f, err := os.OpenFile(filepath.Join(os.TempDir(), "genome-dashboard.log"), os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0o644); err == nil {
		log.SetOutput(multi{os.Stderr, f})
	}
}

type multi []*os.File

func (m multi) Write(p []byte) (int, error) {
	for _, f := range m {
		f.Write(p)
	}
	return len(p), nil
}
