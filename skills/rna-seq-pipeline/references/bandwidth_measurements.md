| Method                              | Throughput     | Status      |
|-------------------------------------|----------------|-------------|
| 1× curl direct (no proxy)           | 1.09 MB/s      | bad         |
| 1× curl FTP protocol (no proxy)     | 0.13 MB/s      | terrible    |
| 1× curl HTTPS via sing-box proxy    | 1.0–5.2 MB/s   | bursty      |
| **4× parallel curl HTTPS via proxy**| **4.7–5.9 MB/s**| **best**   |
| 8× parallel curl HTTPS via proxy    | 5.4 MB/s       | no benefit  |
| NCBI S3 endpoint (sci-hub style)    | 0.13 MB/s      | terrible    |

Tested 2026-05-19 on CNS cluster:
- proxy: sing-box VLESS+Reality at 127.0.0.1:10809
- ENA host: ftp.sra.ebi.ac.uk (HTTPS)
- target file: ENA SRR21846941_1.fastq.gz (~1.9 GB)

Conclusion: 4 concurrent connections saturate the proxy. Going beyond 4 is
counterproductive (more file handles, identical aggregate throughput).

Implications:
- 100 GB / 5 MB/s ≈ 5.5 hours
- For larger jobs (1+ TB), this proxy is a bottleneck — need faster route or
  to run from a host with whitelisted IP (institutional network).
