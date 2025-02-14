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
rz(-2.6211205) q[0];
sx q[0];
rz(-1.0360798) q[0];
sx q[0];
rz(-2.5007024) q[0];
rz(2.9042397) q[1];
sx q[1];
rz(-0.64164716) q[1];
sx q[1];
rz(-3.0566888) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42006832) q[0];
sx q[0];
rz(-0.86948538) q[0];
sx q[0];
rz(3.0140956) q[0];
rz(-pi) q[1];
rz(0.058481599) q[2];
sx q[2];
rz(-1.5396492) q[2];
sx q[2];
rz(-2.8498788) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7218923) q[1];
sx q[1];
rz(-1.2806935) q[1];
sx q[1];
rz(1.6899403) q[1];
rz(-0.97982429) q[3];
sx q[3];
rz(-1.9035467) q[3];
sx q[3];
rz(2.7036288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.25195965) q[2];
sx q[2];
rz(-1.2026938) q[2];
sx q[2];
rz(1.6331875) q[2];
rz(-0.83186045) q[3];
sx q[3];
rz(-0.32604495) q[3];
sx q[3];
rz(-1.7970201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2107596) q[0];
sx q[0];
rz(-2.5063214) q[0];
sx q[0];
rz(2.8765836) q[0];
rz(2.3459072) q[1];
sx q[1];
rz(-0.34514752) q[1];
sx q[1];
rz(-1.4211242) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9819337) q[0];
sx q[0];
rz(-0.93330416) q[0];
sx q[0];
rz(-2.8106014) q[0];
rz(-2.1876343) q[2];
sx q[2];
rz(-1.9618615) q[2];
sx q[2];
rz(-0.80634889) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7184192) q[1];
sx q[1];
rz(-1.8224026) q[1];
sx q[1];
rz(-3.0965641) q[1];
rz(-pi) q[2];
x q[2];
rz(0.029985241) q[3];
sx q[3];
rz(-1.4032149) q[3];
sx q[3];
rz(-0.94887892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77218324) q[2];
sx q[2];
rz(-0.99516827) q[2];
sx q[2];
rz(0.65917242) q[2];
rz(-2.2199421) q[3];
sx q[3];
rz(-1.042659) q[3];
sx q[3];
rz(1.5825745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7121861) q[0];
sx q[0];
rz(-2.5262008) q[0];
sx q[0];
rz(-1.2834826) q[0];
rz(-0.42690024) q[1];
sx q[1];
rz(-0.59395298) q[1];
sx q[1];
rz(1.9740419) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8030464) q[0];
sx q[0];
rz(-1.9340705) q[0];
sx q[0];
rz(2.7873894) q[0];
x q[1];
rz(-2.3857791) q[2];
sx q[2];
rz(-0.67833704) q[2];
sx q[2];
rz(-0.2663904) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.42347149) q[1];
sx q[1];
rz(-1.8602861) q[1];
sx q[1];
rz(-1.4726588) q[1];
x q[2];
rz(1.7469206) q[3];
sx q[3];
rz(-1.8338037) q[3];
sx q[3];
rz(0.60282487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20649642) q[2];
sx q[2];
rz(-1.4859716) q[2];
sx q[2];
rz(0.11779724) q[2];
rz(-1.3055034) q[3];
sx q[3];
rz(-0.85956231) q[3];
sx q[3];
rz(-0.99087244) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30462581) q[0];
sx q[0];
rz(-2.0898297) q[0];
sx q[0];
rz(3.1090609) q[0];
rz(-2.1578728) q[1];
sx q[1];
rz(-0.81147057) q[1];
sx q[1];
rz(-0.54289114) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98292339) q[0];
sx q[0];
rz(-0.74115935) q[0];
sx q[0];
rz(-0.41038402) q[0];
rz(-pi) q[1];
rz(-3.0383598) q[2];
sx q[2];
rz(-1.8816173) q[2];
sx q[2];
rz(0.23676591) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0263152) q[1];
sx q[1];
rz(-1.8247073) q[1];
sx q[1];
rz(-0.8064958) q[1];
rz(0.832295) q[3];
sx q[3];
rz(-1.4595539) q[3];
sx q[3];
rz(0.55505607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5618318) q[2];
sx q[2];
rz(-1.1299955) q[2];
sx q[2];
rz(2.180991) q[2];
rz(2.1660755) q[3];
sx q[3];
rz(-0.87149182) q[3];
sx q[3];
rz(0.52198854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8303216) q[0];
sx q[0];
rz(-2.4595342) q[0];
sx q[0];
rz(-2.8090546) q[0];
rz(2.5034816) q[1];
sx q[1];
rz(-2.5159914) q[1];
sx q[1];
rz(-0.45613751) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77959427) q[0];
sx q[0];
rz(-1.5080305) q[0];
sx q[0];
rz(-1.7032575) q[0];
rz(-pi) q[1];
x q[1];
rz(2.664723) q[2];
sx q[2];
rz(-0.84706351) q[2];
sx q[2];
rz(-0.70666955) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0597555) q[1];
sx q[1];
rz(-2.55099) q[1];
sx q[1];
rz(1.4103622) q[1];
x q[2];
rz(-1.1061323) q[3];
sx q[3];
rz(-2.439508) q[3];
sx q[3];
rz(2.9891069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0629603) q[2];
sx q[2];
rz(-2.3802064) q[2];
sx q[2];
rz(2.6893993) q[2];
rz(-0.26442987) q[3];
sx q[3];
rz(-2.9954438) q[3];
sx q[3];
rz(1.2368081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0753157) q[0];
sx q[0];
rz(-2.3105268) q[0];
sx q[0];
rz(0.7178632) q[0];
rz(1.588899) q[1];
sx q[1];
rz(-1.6009067) q[1];
sx q[1];
rz(-0.20725651) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.358905) q[0];
sx q[0];
rz(-1.1289294) q[0];
sx q[0];
rz(-1.601041) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4137832) q[2];
sx q[2];
rz(-2.7141389) q[2];
sx q[2];
rz(0.70892109) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0644857) q[1];
sx q[1];
rz(-2.047309) q[1];
sx q[1];
rz(2.6576192) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8317811) q[3];
sx q[3];
rz(-1.4647604) q[3];
sx q[3];
rz(-2.1539254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.76102105) q[2];
sx q[2];
rz(-2.0863159) q[2];
sx q[2];
rz(2.5804248) q[2];
rz(1.0593972) q[3];
sx q[3];
rz(-2.0249517) q[3];
sx q[3];
rz(-1.4580956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92784944) q[0];
sx q[0];
rz(-0.97074592) q[0];
sx q[0];
rz(-0.84681502) q[0];
rz(0.98833409) q[1];
sx q[1];
rz(-1.3761995) q[1];
sx q[1];
rz(-0.83289897) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5079437) q[0];
sx q[0];
rz(-1.6744153) q[0];
sx q[0];
rz(-2.9425462) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9994644) q[2];
sx q[2];
rz(-1.4925537) q[2];
sx q[2];
rz(-2.4852666) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.78645241) q[1];
sx q[1];
rz(-1.9398098) q[1];
sx q[1];
rz(2.2166316) q[1];
rz(-pi) q[2];
rz(-1.6844514) q[3];
sx q[3];
rz(-1.6196005) q[3];
sx q[3];
rz(1.5391853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.46574584) q[2];
sx q[2];
rz(-0.70285672) q[2];
sx q[2];
rz(0.77009002) q[2];
rz(-3.0105524) q[3];
sx q[3];
rz(-1.4821056) q[3];
sx q[3];
rz(1.8946764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34762621) q[0];
sx q[0];
rz(-0.93526953) q[0];
sx q[0];
rz(-2.6499709) q[0];
rz(-2.8288016) q[1];
sx q[1];
rz(-0.28952315) q[1];
sx q[1];
rz(-0.082399592) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0534957) q[0];
sx q[0];
rz(-1.6603357) q[0];
sx q[0];
rz(-0.012774668) q[0];
x q[1];
rz(0.56677197) q[2];
sx q[2];
rz(-0.9563891) q[2];
sx q[2];
rz(-0.64726171) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8696599) q[1];
sx q[1];
rz(-2.9026976) q[1];
sx q[1];
rz(0.82456692) q[1];
rz(-pi) q[2];
rz(-1.8833722) q[3];
sx q[3];
rz(-1.4736611) q[3];
sx q[3];
rz(2.7526906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6975434) q[2];
sx q[2];
rz(-1.9335582) q[2];
sx q[2];
rz(3.1235798) q[2];
rz(0.87120122) q[3];
sx q[3];
rz(-2.8049991) q[3];
sx q[3];
rz(2.474031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1910601) q[0];
sx q[0];
rz(-0.55452269) q[0];
sx q[0];
rz(2.6771123) q[0];
rz(2.65061) q[1];
sx q[1];
rz(-1.6698488) q[1];
sx q[1];
rz(1.6741265) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53290908) q[0];
sx q[0];
rz(-2.1957046) q[0];
sx q[0];
rz(1.5708357) q[0];
x q[1];
rz(-0.51767434) q[2];
sx q[2];
rz(-0.56768226) q[2];
sx q[2];
rz(2.9879251) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.58337823) q[1];
sx q[1];
rz(-1.5102856) q[1];
sx q[1];
rz(1.8192181) q[1];
rz(-pi) q[2];
rz(0.72148607) q[3];
sx q[3];
rz(-2.2389447) q[3];
sx q[3];
rz(0.20340445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1269425) q[2];
sx q[2];
rz(-1.6501959) q[2];
sx q[2];
rz(-1.3313741) q[2];
rz(-0.93975449) q[3];
sx q[3];
rz(-1.0569812) q[3];
sx q[3];
rz(1.1872928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.0494999) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3864256) q[0];
sx q[0];
rz(-1.5706129) q[0];
sx q[0];
rz(1.6205233) q[0];
rz(2.0590354) q[2];
sx q[2];
rz(-1.0647213) q[2];
sx q[2];
rz(-0.33389865) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2503719) q[1];
sx q[1];
rz(-0.97019056) q[1];
sx q[1];
rz(-3.1041464) q[1];
rz(-0.97941937) q[3];
sx q[3];
rz(-1.064075) q[3];
sx q[3];
rz(-0.049185569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.99027571) q[2];
sx q[2];
rz(-1.3070725) q[2];
sx q[2];
rz(-1.6952197) q[2];
rz(-0.15898786) q[3];
sx q[3];
rz(-1.9828601) q[3];
sx q[3];
rz(0.75001636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1405519) q[0];
sx q[0];
rz(-0.97698553) q[0];
sx q[0];
rz(-2.3401596) q[0];
rz(1.4871545) q[1];
sx q[1];
rz(-2.9947037) q[1];
sx q[1];
rz(-0.53302232) q[1];
rz(-2.2904446) q[2];
sx q[2];
rz(-1.0590886) q[2];
sx q[2];
rz(1.5589489) q[2];
rz(0.728353) q[3];
sx q[3];
rz(-0.77941685) q[3];
sx q[3];
rz(0.70913915) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
