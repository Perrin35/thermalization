OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4274974) q[0];
sx q[0];
rz(-0.56199718) q[0];
sx q[0];
rz(-2.9105817) q[0];
rz(0.24569874) q[1];
sx q[1];
rz(-0.45431554) q[1];
sx q[1];
rz(-1.8543724) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79345771) q[0];
sx q[0];
rz(-1.0890111) q[0];
sx q[0];
rz(2.8346377) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76164772) q[2];
sx q[2];
rz(-0.33312329) q[2];
sx q[2];
rz(1.3093349) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8637064) q[1];
sx q[1];
rz(-1.7369629) q[1];
sx q[1];
rz(-1.6707129) q[1];
rz(-pi) q[2];
rz(-1.3955529) q[3];
sx q[3];
rz(-1.6449494) q[3];
sx q[3];
rz(1.7046622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40453688) q[2];
sx q[2];
rz(-2.4344567) q[2];
sx q[2];
rz(-2.7810968) q[2];
rz(-2.0170085) q[3];
sx q[3];
rz(-2.093061) q[3];
sx q[3];
rz(1.8726965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26038134) q[0];
sx q[0];
rz(-0.17558782) q[0];
sx q[0];
rz(-0.65688175) q[0];
rz(2.8086713) q[1];
sx q[1];
rz(-2.0491144) q[1];
sx q[1];
rz(-0.93516707) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2774076) q[0];
sx q[0];
rz(-1.5404697) q[0];
sx q[0];
rz(-0.10459374) q[0];
rz(-0.25721154) q[2];
sx q[2];
rz(-1.4818483) q[2];
sx q[2];
rz(-2.6572029) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6015687) q[1];
sx q[1];
rz(-1.717713) q[1];
sx q[1];
rz(1.6107035) q[1];
rz(1.3627429) q[3];
sx q[3];
rz(-1.2223772) q[3];
sx q[3];
rz(1.5294242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8344581) q[2];
sx q[2];
rz(-1.6259401) q[2];
sx q[2];
rz(-1.2163986) q[2];
rz(-2.6136716) q[3];
sx q[3];
rz(-1.0390176) q[3];
sx q[3];
rz(-1.886604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5288178) q[0];
sx q[0];
rz(-2.3154494) q[0];
sx q[0];
rz(-0.58498996) q[0];
rz(-3.0168369) q[1];
sx q[1];
rz(-2.545732) q[1];
sx q[1];
rz(-1.6927208) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.751469) q[0];
sx q[0];
rz(-1.5678143) q[0];
sx q[0];
rz(-0.0026767038) q[0];
rz(-pi) q[1];
rz(2.9563006) q[2];
sx q[2];
rz(-0.073436471) q[2];
sx q[2];
rz(-2.9543608) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6022012) q[1];
sx q[1];
rz(-0.81967205) q[1];
sx q[1];
rz(-0.96266268) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65421974) q[3];
sx q[3];
rz(-2.7252135) q[3];
sx q[3];
rz(-0.49921303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2584194) q[2];
sx q[2];
rz(-2.1397739) q[2];
sx q[2];
rz(-1.6955356) q[2];
rz(-2.3840733) q[3];
sx q[3];
rz(-1.5816553) q[3];
sx q[3];
rz(-2.3022046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7739173) q[0];
sx q[0];
rz(-2.2125419) q[0];
sx q[0];
rz(-1.0876592) q[0];
rz(-2.7663973) q[1];
sx q[1];
rz(-2.0162069) q[1];
sx q[1];
rz(-0.96022022) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69273573) q[0];
sx q[0];
rz(-1.6105284) q[0];
sx q[0];
rz(1.3055152) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4713418) q[2];
sx q[2];
rz(-1.5410454) q[2];
sx q[2];
rz(-3.0029675) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.5615446) q[1];
sx q[1];
rz(-2.5001038) q[1];
sx q[1];
rz(-2.3381691) q[1];
x q[2];
rz(0.26953617) q[3];
sx q[3];
rz(-1.663066) q[3];
sx q[3];
rz(2.5470981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.20597657) q[2];
sx q[2];
rz(-1.874186) q[2];
sx q[2];
rz(1.4975632) q[2];
rz(-0.86132541) q[3];
sx q[3];
rz(-0.70278168) q[3];
sx q[3];
rz(2.6845045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.977026) q[0];
sx q[0];
rz(-0.67104665) q[0];
sx q[0];
rz(-0.60428756) q[0];
rz(1.1445649) q[1];
sx q[1];
rz(-1.6555758) q[1];
sx q[1];
rz(0.42246517) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0351663) q[0];
sx q[0];
rz(-1.6521593) q[0];
sx q[0];
rz(1.401713) q[0];
rz(2.1429135) q[2];
sx q[2];
rz(-1.0829751) q[2];
sx q[2];
rz(-1.7833379) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3389272) q[1];
sx q[1];
rz(-2.4484854) q[1];
sx q[1];
rz(0.1130123) q[1];
x q[2];
rz(-2.7606008) q[3];
sx q[3];
rz(-1.89011) q[3];
sx q[3];
rz(-1.3407941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3182688) q[2];
sx q[2];
rz(-0.6508998) q[2];
sx q[2];
rz(0.65625119) q[2];
rz(-1.6107791) q[3];
sx q[3];
rz(-1.2698413) q[3];
sx q[3];
rz(0.58073616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73606473) q[0];
sx q[0];
rz(-1.0117714) q[0];
sx q[0];
rz(2.5166125) q[0];
rz(-0.75195733) q[1];
sx q[1];
rz(-2.0729005) q[1];
sx q[1];
rz(-1.6065074) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84941712) q[0];
sx q[0];
rz(-1.2818953) q[0];
sx q[0];
rz(1.6221415) q[0];
rz(-0.97246031) q[2];
sx q[2];
rz(-2.1876799) q[2];
sx q[2];
rz(-1.0765558) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2481164) q[1];
sx q[1];
rz(-1.7991589) q[1];
sx q[1];
rz(-0.62947692) q[1];
rz(-pi) q[2];
rz(2.1702437) q[3];
sx q[3];
rz(-2.0347715) q[3];
sx q[3];
rz(2.6148207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9525166) q[2];
sx q[2];
rz(-1.357888) q[2];
sx q[2];
rz(-2.2124186) q[2];
rz(2.3164228) q[3];
sx q[3];
rz(-1.630183) q[3];
sx q[3];
rz(-2.0391803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55649844) q[0];
sx q[0];
rz(-2.3346021) q[0];
sx q[0];
rz(1.7671385) q[0];
rz(0.51042405) q[1];
sx q[1];
rz(-0.7904895) q[1];
sx q[1];
rz(-0.62320954) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0334367) q[0];
sx q[0];
rz(-1.2386453) q[0];
sx q[0];
rz(2.9812814) q[0];
rz(-0.34131949) q[2];
sx q[2];
rz(-2.6067197) q[2];
sx q[2];
rz(0.90082263) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0859061) q[1];
sx q[1];
rz(-1.4594363) q[1];
sx q[1];
rz(-1.423905) q[1];
x q[2];
rz(1.3680063) q[3];
sx q[3];
rz(-0.8335127) q[3];
sx q[3];
rz(1.1356419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0509384) q[2];
sx q[2];
rz(-0.81523681) q[2];
sx q[2];
rz(-1.9742924) q[2];
rz(-3.1189611) q[3];
sx q[3];
rz(-0.52112094) q[3];
sx q[3];
rz(-1.1289271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8959494) q[0];
sx q[0];
rz(-2.0063945) q[0];
sx q[0];
rz(-2.1242712) q[0];
rz(1.7474878) q[1];
sx q[1];
rz(-2.930495) q[1];
sx q[1];
rz(-0.62754935) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.548259) q[0];
sx q[0];
rz(-0.98932668) q[0];
sx q[0];
rz(0.47971804) q[0];
rz(-2.4623221) q[2];
sx q[2];
rz(-1.5049266) q[2];
sx q[2];
rz(0.010347376) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1422954) q[1];
sx q[1];
rz(-0.31041708) q[1];
sx q[1];
rz(0.69119549) q[1];
x q[2];
rz(2.5379337) q[3];
sx q[3];
rz(-0.55760477) q[3];
sx q[3];
rz(1.738036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.62577406) q[2];
sx q[2];
rz(-0.63483441) q[2];
sx q[2];
rz(-2.7308357) q[2];
rz(-0.050203236) q[3];
sx q[3];
rz(-1.8886731) q[3];
sx q[3];
rz(-3.0661809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.582616) q[0];
sx q[0];
rz(-2.2447383) q[0];
sx q[0];
rz(-2.8726752) q[0];
rz(0.31002632) q[1];
sx q[1];
rz(-1.0870442) q[1];
sx q[1];
rz(0.1300098) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1300505) q[0];
sx q[0];
rz(-0.39549144) q[0];
sx q[0];
rz(0.96590913) q[0];
rz(-pi) q[1];
rz(-2.347888) q[2];
sx q[2];
rz(-0.48173258) q[2];
sx q[2];
rz(-2.6793753) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8925608) q[1];
sx q[1];
rz(-1.7317803) q[1];
sx q[1];
rz(2.575483) q[1];
x q[2];
rz(-0.9762398) q[3];
sx q[3];
rz(-0.77139445) q[3];
sx q[3];
rz(2.6736128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.93398634) q[2];
sx q[2];
rz(-2.3748368) q[2];
sx q[2];
rz(-0.096573528) q[2];
rz(1.5038331) q[3];
sx q[3];
rz(-1.8042754) q[3];
sx q[3];
rz(2.3418929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12829517) q[0];
sx q[0];
rz(-0.18823637) q[0];
sx q[0];
rz(-2.6085594) q[0];
rz(0.036272613) q[1];
sx q[1];
rz(-2.3590922) q[1];
sx q[1];
rz(1.689555) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8743458) q[0];
sx q[0];
rz(-1.1935992) q[0];
sx q[0];
rz(-2.1826571) q[0];
rz(-1.6797941) q[2];
sx q[2];
rz(-1.8882635) q[2];
sx q[2];
rz(0.017680971) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2382792) q[1];
sx q[1];
rz(-2.7809445) q[1];
sx q[1];
rz(-0.084899501) q[1];
rz(-pi) q[2];
rz(-2.6759869) q[3];
sx q[3];
rz(-1.4911663) q[3];
sx q[3];
rz(-1.3543983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71296802) q[2];
sx q[2];
rz(-2.4098101) q[2];
sx q[2];
rz(-3.1089605) q[2];
rz(2.2964358) q[3];
sx q[3];
rz(-1.9149575) q[3];
sx q[3];
rz(-2.0613861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3451155) q[0];
sx q[0];
rz(-2.4892172) q[0];
sx q[0];
rz(-0.31443483) q[0];
rz(0.79700094) q[1];
sx q[1];
rz(-0.92756699) q[1];
sx q[1];
rz(-3.0753593) q[1];
rz(2.6673139) q[2];
sx q[2];
rz(-2.5534292) q[2];
sx q[2];
rz(-3.0234887) q[2];
rz(1.6975523) q[3];
sx q[3];
rz(-0.76382617) q[3];
sx q[3];
rz(2.7504117) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
