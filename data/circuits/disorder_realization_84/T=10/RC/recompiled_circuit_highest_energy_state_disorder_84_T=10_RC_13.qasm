OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.52047211) q[0];
sx q[0];
rz(-2.1055129) q[0];
sx q[0];
rz(2.5007024) q[0];
rz(2.9042397) q[1];
sx q[1];
rz(-0.64164716) q[1];
sx q[1];
rz(0.084903804) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2239421) q[0];
sx q[0];
rz(-0.71086796) q[0];
sx q[0];
rz(-1.720233) q[0];
rz(-0.058481599) q[2];
sx q[2];
rz(-1.6019434) q[2];
sx q[2];
rz(-2.8498788) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0247268) q[1];
sx q[1];
rz(-1.6849396) q[1];
sx q[1];
rz(-2.8495339) q[1];
x q[2];
rz(-0.39438168) q[3];
sx q[3];
rz(-2.1254) q[3];
sx q[3];
rz(1.3485874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.25195965) q[2];
sx q[2];
rz(-1.2026938) q[2];
sx q[2];
rz(-1.6331875) q[2];
rz(-0.83186045) q[3];
sx q[3];
rz(-0.32604495) q[3];
sx q[3];
rz(1.3445725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2107596) q[0];
sx q[0];
rz(-2.5063214) q[0];
sx q[0];
rz(-2.8765836) q[0];
rz(-2.3459072) q[1];
sx q[1];
rz(-0.34514752) q[1];
sx q[1];
rz(1.4211242) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6832631) q[0];
sx q[0];
rz(-2.4340672) q[0];
sx q[0];
rz(1.9842771) q[0];
rz(-pi) q[1];
rz(-0.46799731) q[2];
sx q[2];
rz(-2.1350522) q[2];
sx q[2];
rz(2.6411438) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2441229) q[1];
sx q[1];
rz(-0.25551968) q[1];
sx q[1];
rz(-1.3974415) q[1];
rz(-pi) q[2];
rz(-0.029985241) q[3];
sx q[3];
rz(-1.4032149) q[3];
sx q[3];
rz(-2.1927137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3694094) q[2];
sx q[2];
rz(-0.99516827) q[2];
sx q[2];
rz(2.4824202) q[2];
rz(2.2199421) q[3];
sx q[3];
rz(-2.0989336) q[3];
sx q[3];
rz(1.5825745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7121861) q[0];
sx q[0];
rz(-2.5262008) q[0];
sx q[0];
rz(-1.8581101) q[0];
rz(0.42690024) q[1];
sx q[1];
rz(-2.5476397) q[1];
sx q[1];
rz(1.9740419) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99770633) q[0];
sx q[0];
rz(-2.6397815) q[0];
sx q[0];
rz(0.83117475) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3857791) q[2];
sx q[2];
rz(-0.67833704) q[2];
sx q[2];
rz(2.8752022) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.75558799) q[1];
sx q[1];
rz(-0.30522871) q[1];
sx q[1];
rz(-2.8237923) q[1];
rz(-pi) q[2];
rz(-0.57691388) q[3];
sx q[3];
rz(-2.8261999) q[3];
sx q[3];
rz(-1.9384991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9350962) q[2];
sx q[2];
rz(-1.4859716) q[2];
sx q[2];
rz(-0.11779724) q[2];
rz(-1.3055034) q[3];
sx q[3];
rz(-0.85956231) q[3];
sx q[3];
rz(-0.99087244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8369668) q[0];
sx q[0];
rz(-2.0898297) q[0];
sx q[0];
rz(-0.032531746) q[0];
rz(0.98371983) q[1];
sx q[1];
rz(-2.3301221) q[1];
sx q[1];
rz(-2.5987015) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27730091) q[0];
sx q[0];
rz(-1.8435209) q[0];
sx q[0];
rz(-0.69820531) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2603733) q[2];
sx q[2];
rz(-0.32698787) q[2];
sx q[2];
rz(2.5782158) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.429814) q[1];
sx q[1];
rz(-0.79718726) q[1];
sx q[1];
rz(-1.2120257) q[1];
rz(-pi) q[2];
rz(1.4063605) q[3];
sx q[3];
rz(-2.3963266) q[3];
sx q[3];
rz(1.1370873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5797609) q[2];
sx q[2];
rz(-1.1299955) q[2];
sx q[2];
rz(-0.96060166) q[2];
rz(0.97551712) q[3];
sx q[3];
rz(-2.2701008) q[3];
sx q[3];
rz(0.52198854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3112711) q[0];
sx q[0];
rz(-2.4595342) q[0];
sx q[0];
rz(-0.3325381) q[0];
rz(0.63811103) q[1];
sx q[1];
rz(-2.5159914) q[1];
sx q[1];
rz(0.45613751) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7903731) q[0];
sx q[0];
rz(-0.14650061) q[0];
sx q[0];
rz(1.1266493) q[0];
rz(0.78807414) q[2];
sx q[2];
rz(-1.2196793) q[2];
sx q[2];
rz(1.9478363) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3553473) q[1];
sx q[1];
rz(-1.4817217) q[1];
sx q[1];
rz(-2.1554355) q[1];
x q[2];
rz(-2.2182369) q[3];
sx q[3];
rz(-1.277194) q[3];
sx q[3];
rz(-1.3577611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0629603) q[2];
sx q[2];
rz(-0.76138622) q[2];
sx q[2];
rz(2.6893993) q[2];
rz(-0.26442987) q[3];
sx q[3];
rz(-0.14614883) q[3];
sx q[3];
rz(1.9047846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0662769) q[0];
sx q[0];
rz(-0.83106581) q[0];
sx q[0];
rz(0.7178632) q[0];
rz(-1.5526937) q[1];
sx q[1];
rz(-1.6009067) q[1];
sx q[1];
rz(2.9343361) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.358905) q[0];
sx q[0];
rz(-2.0126632) q[0];
sx q[0];
rz(1.601041) q[0];
rz(-pi) q[1];
rz(-1.7278094) q[2];
sx q[2];
rz(-2.7141389) q[2];
sx q[2];
rz(-2.4326716) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8718827) q[1];
sx q[1];
rz(-1.1445395) q[1];
sx q[1];
rz(2.0987418) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8317811) q[3];
sx q[3];
rz(-1.6768322) q[3];
sx q[3];
rz(0.98766726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.76102105) q[2];
sx q[2];
rz(-1.0552768) q[2];
sx q[2];
rz(-2.5804248) q[2];
rz(1.0593972) q[3];
sx q[3];
rz(-1.1166409) q[3];
sx q[3];
rz(1.4580956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2137432) q[0];
sx q[0];
rz(-2.1708467) q[0];
sx q[0];
rz(0.84681502) q[0];
rz(2.1532586) q[1];
sx q[1];
rz(-1.7653932) q[1];
sx q[1];
rz(-0.83289897) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7306379) q[0];
sx q[0];
rz(-0.22408501) q[0];
sx q[0];
rz(-2.6574294) q[0];
rz(-1.49176) q[2];
sx q[2];
rz(-1.7124868) q[2];
sx q[2];
rz(2.2159383) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3551402) q[1];
sx q[1];
rz(-1.9398098) q[1];
sx q[1];
rz(-2.2166316) q[1];
rz(0.049120589) q[3];
sx q[3];
rz(-1.4572772) q[3];
sx q[3];
rz(0.037179557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46574584) q[2];
sx q[2];
rz(-2.4387359) q[2];
sx q[2];
rz(-2.3715026) q[2];
rz(-0.13104023) q[3];
sx q[3];
rz(-1.4821056) q[3];
sx q[3];
rz(-1.8946764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7939664) q[0];
sx q[0];
rz(-0.93526953) q[0];
sx q[0];
rz(-2.6499709) q[0];
rz(-2.8288016) q[1];
sx q[1];
rz(-0.28952315) q[1];
sx q[1];
rz(3.0591931) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9461878) q[0];
sx q[0];
rz(-3.051149) q[0];
sx q[0];
rz(1.7121332) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2213242) q[2];
sx q[2];
rz(-0.81020497) q[2];
sx q[2];
rz(-1.4820463) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1090549) q[1];
sx q[1];
rz(-1.3961642) q[1];
sx q[1];
rz(-0.1638549) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8776603) q[3];
sx q[3];
rz(-0.32684775) q[3];
sx q[3];
rz(-0.89034789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4440492) q[2];
sx q[2];
rz(-1.9335582) q[2];
sx q[2];
rz(3.1235798) q[2];
rz(2.2703914) q[3];
sx q[3];
rz(-2.8049991) q[3];
sx q[3];
rz(0.66756162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1910601) q[0];
sx q[0];
rz(-0.55452269) q[0];
sx q[0];
rz(-0.46448034) q[0];
rz(0.49098268) q[1];
sx q[1];
rz(-1.6698488) q[1];
sx q[1];
rz(-1.6741265) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6086836) q[0];
sx q[0];
rz(-0.94588806) q[0];
sx q[0];
rz(1.5708357) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50602127) q[2];
sx q[2];
rz(-1.3014761) q[2];
sx q[2];
rz(-1.276818) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58337823) q[1];
sx q[1];
rz(-1.5102856) q[1];
sx q[1];
rz(1.3223745) q[1];
x q[2];
rz(-0.72148607) q[3];
sx q[3];
rz(-0.90264796) q[3];
sx q[3];
rz(-2.9381882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.014650194) q[2];
sx q[2];
rz(-1.6501959) q[2];
sx q[2];
rz(-1.8102185) q[2];
rz(2.2018382) q[3];
sx q[3];
rz(-2.0846114) q[3];
sx q[3];
rz(-1.1872928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2347655) q[0];
sx q[0];
rz(-2.6563788) q[0];
sx q[0];
rz(1.0377129) q[0];
rz(-2.0872133) q[1];
sx q[1];
rz(-1.8260006) q[1];
sx q[1];
rz(2.0920928) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9572311) q[0];
sx q[0];
rz(-1.5210694) q[0];
sx q[0];
rz(3.141409) q[0];
rz(0.56042258) q[2];
sx q[2];
rz(-1.9935521) q[2];
sx q[2];
rz(-2.1566856) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8250268) q[1];
sx q[1];
rz(-2.5399643) q[1];
sx q[1];
rz(-1.5161995) q[1];
rz(0.97941937) q[3];
sx q[3];
rz(-2.0775177) q[3];
sx q[3];
rz(3.0924071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1513169) q[2];
sx q[2];
rz(-1.8345202) q[2];
sx q[2];
rz(1.6952197) q[2];
rz(-0.15898786) q[3];
sx q[3];
rz(-1.9828601) q[3];
sx q[3];
rz(0.75001636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0010407) q[0];
sx q[0];
rz(-0.97698553) q[0];
sx q[0];
rz(-2.3401596) q[0];
rz(-1.4871545) q[1];
sx q[1];
rz(-0.14688891) q[1];
sx q[1];
rz(2.6085703) q[1];
rz(0.86511378) q[2];
sx q[2];
rz(-0.8556753) q[2];
sx q[2];
rz(-0.52134261) q[2];
rz(0.63538649) q[3];
sx q[3];
rz(-1.0839331) q[3];
sx q[3];
rz(-1.4270368) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
