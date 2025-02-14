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
rz(-0.44218818) q[0];
sx q[0];
rz(-1.1531885) q[0];
sx q[0];
rz(-0.12864223) q[0];
rz(1.6755942) q[1];
sx q[1];
rz(-2.0425551) q[1];
sx q[1];
rz(1.0817945) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3386782) q[0];
sx q[0];
rz(-0.071597425) q[0];
sx q[0];
rz(2.6971571) q[0];
rz(-pi) q[1];
rz(0.84613447) q[2];
sx q[2];
rz(-2.6577009) q[2];
sx q[2];
rz(-2.0563375) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9777367) q[1];
sx q[1];
rz(-2.9675936) q[1];
sx q[1];
rz(-1.15088) q[1];
rz(-pi) q[2];
rz(-1.2461189) q[3];
sx q[3];
rz(-2.8522308) q[3];
sx q[3];
rz(-2.112767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.094540207) q[2];
sx q[2];
rz(-2.4461942) q[2];
sx q[2];
rz(2.409234) q[2];
rz(0.098946027) q[3];
sx q[3];
rz(-2.1061335) q[3];
sx q[3];
rz(-1.8100479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.052213) q[0];
sx q[0];
rz(-2.3591924) q[0];
sx q[0];
rz(3.1001477) q[0];
rz(-0.6913569) q[1];
sx q[1];
rz(-1.3203011) q[1];
sx q[1];
rz(0.88821214) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23032886) q[0];
sx q[0];
rz(-2.3155876) q[0];
sx q[0];
rz(0.86358304) q[0];
rz(-pi) q[1];
rz(0.67866831) q[2];
sx q[2];
rz(-1.1801085) q[2];
sx q[2];
rz(-0.44579577) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6710818) q[1];
sx q[1];
rz(-2.8948519) q[1];
sx q[1];
rz(-0.17784878) q[1];
x q[2];
rz(0.4732186) q[3];
sx q[3];
rz(-0.56269533) q[3];
sx q[3];
rz(2.0317047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1516252) q[2];
sx q[2];
rz(-1.7513559) q[2];
sx q[2];
rz(-1.3803231) q[2];
rz(3.0917998) q[3];
sx q[3];
rz(-0.68792611) q[3];
sx q[3];
rz(0.5602347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77408537) q[0];
sx q[0];
rz(-2.9893576) q[0];
sx q[0];
rz(1.2138858) q[0];
rz(1.8269352) q[1];
sx q[1];
rz(-2.1036802) q[1];
sx q[1];
rz(1.5710057) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9938068) q[0];
sx q[0];
rz(-1.3249614) q[0];
sx q[0];
rz(0.39061201) q[0];
rz(-2.7309787) q[2];
sx q[2];
rz(-1.7143692) q[2];
sx q[2];
rz(-0.5856572) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27109001) q[1];
sx q[1];
rz(-1.6047641) q[1];
sx q[1];
rz(1.4064596) q[1];
x q[2];
rz(2.5475471) q[3];
sx q[3];
rz(-2.0875153) q[3];
sx q[3];
rz(1.3498211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7870002) q[2];
sx q[2];
rz(-0.74942333) q[2];
sx q[2];
rz(-0.3053537) q[2];
rz(2.2954156) q[3];
sx q[3];
rz(-1.2140707) q[3];
sx q[3];
rz(-0.86887104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2884035) q[0];
sx q[0];
rz(-1.0193634) q[0];
sx q[0];
rz(0.97002059) q[0];
rz(2.8963529) q[1];
sx q[1];
rz(-1.0673362) q[1];
sx q[1];
rz(1.5018357) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63642913) q[0];
sx q[0];
rz(-1.7254524) q[0];
sx q[0];
rz(1.4419492) q[0];
x q[1];
rz(-1.3972598) q[2];
sx q[2];
rz(-1.3385941) q[2];
sx q[2];
rz(0.47877889) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8322595) q[1];
sx q[1];
rz(-0.24411476) q[1];
sx q[1];
rz(-0.67863087) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7795622) q[3];
sx q[3];
rz(-2.3811582) q[3];
sx q[3];
rz(2.8575051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.60845) q[2];
sx q[2];
rz(-1.5773982) q[2];
sx q[2];
rz(-1.8428141) q[2];
rz(-0.53457824) q[3];
sx q[3];
rz(-0.60408533) q[3];
sx q[3];
rz(0.66876137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0534441) q[0];
sx q[0];
rz(-1.3548594) q[0];
sx q[0];
rz(2.163929) q[0];
rz(0.68556249) q[1];
sx q[1];
rz(-0.79427636) q[1];
sx q[1];
rz(-0.96644863) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9622877) q[0];
sx q[0];
rz(-2.0317269) q[0];
sx q[0];
rz(2.3910752) q[0];
x q[1];
rz(0.62018779) q[2];
sx q[2];
rz(-2.475707) q[2];
sx q[2];
rz(0.87845368) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9453166) q[1];
sx q[1];
rz(-1.4320201) q[1];
sx q[1];
rz(0.86244418) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2174986) q[3];
sx q[3];
rz(-1.0477001) q[3];
sx q[3];
rz(-2.9016665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55296772) q[2];
sx q[2];
rz(-1.3543465) q[2];
sx q[2];
rz(-2.8589613) q[2];
rz(0.96453729) q[3];
sx q[3];
rz(-1.5805565) q[3];
sx q[3];
rz(0.36981043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.033217) q[0];
sx q[0];
rz(-2.7300457) q[0];
sx q[0];
rz(-2.0630398) q[0];
rz(-0.82303965) q[1];
sx q[1];
rz(-2.0185202) q[1];
sx q[1];
rz(-0.80026921) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6155933) q[0];
sx q[0];
rz(-1.5999927) q[0];
sx q[0];
rz(-3.133081) q[0];
rz(-0.79485017) q[2];
sx q[2];
rz(-0.42758105) q[2];
sx q[2];
rz(-2.7045369) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2507402) q[1];
sx q[1];
rz(-1.0167334) q[1];
sx q[1];
rz(2.3046369) q[1];
x q[2];
rz(-2.0144281) q[3];
sx q[3];
rz(-0.95329696) q[3];
sx q[3];
rz(1.3408183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0018953) q[2];
sx q[2];
rz(-1.8997833) q[2];
sx q[2];
rz(-0.099420698) q[2];
rz(-2.5771778) q[3];
sx q[3];
rz(-1.5921009) q[3];
sx q[3];
rz(-2.2216589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1426706) q[0];
sx q[0];
rz(-2.8746334) q[0];
sx q[0];
rz(0.024918407) q[0];
rz(1.092356) q[1];
sx q[1];
rz(-1.9905636) q[1];
sx q[1];
rz(-2.5999462) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9077907) q[0];
sx q[0];
rz(-0.72392091) q[0];
sx q[0];
rz(2.1320599) q[0];
rz(-1.0850026) q[2];
sx q[2];
rz(-2.5143904) q[2];
sx q[2];
rz(1.021334) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4389651) q[1];
sx q[1];
rz(-2.781233) q[1];
sx q[1];
rz(1.9907336) q[1];
x q[2];
rz(-2.1775167) q[3];
sx q[3];
rz(-1.4243444) q[3];
sx q[3];
rz(1.7646542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.8267374) q[2];
sx q[2];
rz(-1.22437) q[2];
sx q[2];
rz(1.1401736) q[2];
rz(1.4402116) q[3];
sx q[3];
rz(-2.7482016) q[3];
sx q[3];
rz(-0.17471974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88406554) q[0];
sx q[0];
rz(-1.1648357) q[0];
sx q[0];
rz(-0.87673941) q[0];
rz(2.9675617) q[1];
sx q[1];
rz(-2.048025) q[1];
sx q[1];
rz(1.7736951) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48094952) q[0];
sx q[0];
rz(-2.5202453) q[0];
sx q[0];
rz(2.0578007) q[0];
x q[1];
rz(2.9139745) q[2];
sx q[2];
rz(-1.2151067) q[2];
sx q[2];
rz(1.5700036) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2217283) q[1];
sx q[1];
rz(-2.113453) q[1];
sx q[1];
rz(-2.7046287) q[1];
rz(0.075966751) q[3];
sx q[3];
rz(-1.6923762) q[3];
sx q[3];
rz(1.0978804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6690663) q[2];
sx q[2];
rz(-1.8146699) q[2];
sx q[2];
rz(-2.9533022) q[2];
rz(-1.3663728) q[3];
sx q[3];
rz(-1.3683189) q[3];
sx q[3];
rz(1.202047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5240204) q[0];
sx q[0];
rz(-3.0078648) q[0];
sx q[0];
rz(2.4753841) q[0];
rz(-1.1353525) q[1];
sx q[1];
rz(-2.1609781) q[1];
sx q[1];
rz(1.1599783) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2839417) q[0];
sx q[0];
rz(-0.91383024) q[0];
sx q[0];
rz(0.82842555) q[0];
rz(-pi) q[1];
rz(0.23187821) q[2];
sx q[2];
rz(-2.3740951) q[2];
sx q[2];
rz(2.9434443) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4593606) q[1];
sx q[1];
rz(-2.5256093) q[1];
sx q[1];
rz(-1.5840696) q[1];
rz(-pi) q[2];
rz(2.8256666) q[3];
sx q[3];
rz(-2.1135309) q[3];
sx q[3];
rz(2.3019615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3975415) q[2];
sx q[2];
rz(-0.56637374) q[2];
sx q[2];
rz(0.19691021) q[2];
rz(-2.2714553) q[3];
sx q[3];
rz(-1.0715485) q[3];
sx q[3];
rz(-1.2675233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63721913) q[0];
sx q[0];
rz(-1.0727896) q[0];
sx q[0];
rz(3.1043501) q[0];
rz(1.0875018) q[1];
sx q[1];
rz(-1.0481513) q[1];
sx q[1];
rz(0.16669272) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.113362) q[0];
sx q[0];
rz(-2.2489002) q[0];
sx q[0];
rz(1.9166758) q[0];
x q[1];
rz(-1.8718821) q[2];
sx q[2];
rz(-1.1043806) q[2];
sx q[2];
rz(0.25638858) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6065258) q[1];
sx q[1];
rz(-1.6589612) q[1];
sx q[1];
rz(2.8463581) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14504542) q[3];
sx q[3];
rz(-1.5022455) q[3];
sx q[3];
rz(-1.1762432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.4166261) q[2];
sx q[2];
rz(-2.5278957) q[2];
sx q[2];
rz(1.7566768) q[2];
rz(-2.7409605) q[3];
sx q[3];
rz(-1.8568361) q[3];
sx q[3];
rz(-2.1304776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8818125) q[0];
sx q[0];
rz(-2.0406944) q[0];
sx q[0];
rz(2.5115321) q[0];
rz(-0.13314816) q[1];
sx q[1];
rz(-1.1590191) q[1];
sx q[1];
rz(-1.0551183) q[1];
rz(-1.3721977) q[2];
sx q[2];
rz(-0.19459859) q[2];
sx q[2];
rz(2.9698402) q[2];
rz(-1.4478441) q[3];
sx q[3];
rz(-2.3364441) q[3];
sx q[3];
rz(0.91329109) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
