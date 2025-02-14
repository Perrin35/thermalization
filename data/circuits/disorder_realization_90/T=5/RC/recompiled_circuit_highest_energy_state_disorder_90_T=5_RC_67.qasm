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
rz(-0.65547216) q[0];
sx q[0];
rz(4.0937427) q[0];
sx q[0];
rz(9.5185315) q[0];
rz(-1.7733511) q[1];
sx q[1];
rz(-1.5612839) q[1];
sx q[1];
rz(0.66722792) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8225559) q[0];
sx q[0];
rz(-1.0780436) q[0];
sx q[0];
rz(2.7125554) q[0];
x q[1];
rz(-2.480443) q[2];
sx q[2];
rz(-2.8183658) q[2];
sx q[2];
rz(1.8913392) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.88244438) q[1];
sx q[1];
rz(-1.8516774) q[1];
sx q[1];
rz(0.85390635) q[1];
rz(0.85140075) q[3];
sx q[3];
rz(-0.82653294) q[3];
sx q[3];
rz(1.1917308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61624709) q[2];
sx q[2];
rz(-1.5364001) q[2];
sx q[2];
rz(3.1180535) q[2];
rz(2.8871138) q[3];
sx q[3];
rz(-1.9813709) q[3];
sx q[3];
rz(-0.75828534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94754058) q[0];
sx q[0];
rz(-1.7451311) q[0];
sx q[0];
rz(-2.5260455) q[0];
rz(-2.5415892) q[1];
sx q[1];
rz(-1.871385) q[1];
sx q[1];
rz(2.792865) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.00086) q[0];
sx q[0];
rz(-1.6475186) q[0];
sx q[0];
rz(-0.64312913) q[0];
rz(-pi) q[1];
rz(-2.8519657) q[2];
sx q[2];
rz(-1.8995325) q[2];
sx q[2];
rz(1.5266974) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.80417577) q[1];
sx q[1];
rz(-1.2722208) q[1];
sx q[1];
rz(2.0470624) q[1];
rz(-pi) q[2];
rz(2.8564151) q[3];
sx q[3];
rz(-2.0582407) q[3];
sx q[3];
rz(0.91966682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.49150026) q[2];
sx q[2];
rz(-1.8365752) q[2];
sx q[2];
rz(2.4864062) q[2];
rz(0.16264597) q[3];
sx q[3];
rz(-0.9442257) q[3];
sx q[3];
rz(-2.9338525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.1732037) q[0];
sx q[0];
rz(-2.0252616) q[0];
sx q[0];
rz(-3.047347) q[0];
rz(1.8415797) q[1];
sx q[1];
rz(-0.899122) q[1];
sx q[1];
rz(-1.1945486) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0604977) q[0];
sx q[0];
rz(-1.0733685) q[0];
sx q[0];
rz(1.917385) q[0];
x q[1];
rz(3.0897806) q[2];
sx q[2];
rz(-1.2070884) q[2];
sx q[2];
rz(1.5075945) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0654203) q[1];
sx q[1];
rz(-1.3143871) q[1];
sx q[1];
rz(-0.045412046) q[1];
rz(-pi) q[2];
rz(-0.39154176) q[3];
sx q[3];
rz(-0.80219183) q[3];
sx q[3];
rz(-1.4456309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.82075787) q[2];
sx q[2];
rz(-1.6944378) q[2];
sx q[2];
rz(2.7336332) q[2];
rz(-1.3784493) q[3];
sx q[3];
rz(-2.2429376) q[3];
sx q[3];
rz(0.11817008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5648062) q[0];
sx q[0];
rz(-1.7498359) q[0];
sx q[0];
rz(2.3396662) q[0];
rz(0.39464125) q[1];
sx q[1];
rz(-1.0486187) q[1];
sx q[1];
rz(-3.1065497) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5749767) q[0];
sx q[0];
rz(-2.1144951) q[0];
sx q[0];
rz(0.4369785) q[0];
x q[1];
rz(-1.7614683) q[2];
sx q[2];
rz(-1.9004656) q[2];
sx q[2];
rz(-0.53527385) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9600507) q[1];
sx q[1];
rz(-1.6142577) q[1];
sx q[1];
rz(3.0994013) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48864103) q[3];
sx q[3];
rz(-2.7436769) q[3];
sx q[3];
rz(-0.15109381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0617712) q[2];
sx q[2];
rz(-2.0786736) q[2];
sx q[2];
rz(1.4351832) q[2];
rz(1.4052514) q[3];
sx q[3];
rz(-1.0900213) q[3];
sx q[3];
rz(-2.0315571) q[3];
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
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46471304) q[0];
sx q[0];
rz(-1.1486624) q[0];
sx q[0];
rz(-1.1418463) q[0];
rz(2.6308718) q[1];
sx q[1];
rz(-1.2341713) q[1];
sx q[1];
rz(1.7150735) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079581408) q[0];
sx q[0];
rz(-2.1163762) q[0];
sx q[0];
rz(2.1929492) q[0];
rz(-pi) q[1];
rz(1.8829338) q[2];
sx q[2];
rz(-1.3545879) q[2];
sx q[2];
rz(2.5463111) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1514476) q[1];
sx q[1];
rz(-1.6832608) q[1];
sx q[1];
rz(-0.97655762) q[1];
rz(-pi) q[2];
rz(-0.1673836) q[3];
sx q[3];
rz(-1.1684844) q[3];
sx q[3];
rz(1.2237807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7412173) q[2];
sx q[2];
rz(-1.3072689) q[2];
sx q[2];
rz(2.884414) q[2];
rz(-0.87279618) q[3];
sx q[3];
rz(-1.8200487) q[3];
sx q[3];
rz(2.0040373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
rz(-2.0208825) q[0];
sx q[0];
rz(-0.45419422) q[0];
sx q[0];
rz(1.0783476) q[0];
rz(0.9067761) q[1];
sx q[1];
rz(-2.4557143) q[1];
sx q[1];
rz(3.0996941) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2611194) q[0];
sx q[0];
rz(-0.47396892) q[0];
sx q[0];
rz(1.4920477) q[0];
rz(-pi) q[1];
rz(-2.2380377) q[2];
sx q[2];
rz(-1.6674526) q[2];
sx q[2];
rz(-2.5943611) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1160526) q[1];
sx q[1];
rz(-1.2092822) q[1];
sx q[1];
rz(-0.29360969) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.83591977) q[3];
sx q[3];
rz(-1.984388) q[3];
sx q[3];
rz(1.6294162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28356734) q[2];
sx q[2];
rz(-2.4391386) q[2];
sx q[2];
rz(-2.2502327) q[2];
rz(2.4274965) q[3];
sx q[3];
rz(-0.73914206) q[3];
sx q[3];
rz(-1.063063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5080268) q[0];
sx q[0];
rz(-1.897568) q[0];
sx q[0];
rz(2.4666393) q[0];
rz(1.5114816) q[1];
sx q[1];
rz(-2.5210896) q[1];
sx q[1];
rz(2.564548) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9723608) q[0];
sx q[0];
rz(-1.0857538) q[0];
sx q[0];
rz(-0.26974704) q[0];
rz(-pi) q[1];
rz(-2.9222832) q[2];
sx q[2];
rz(-2.6990934) q[2];
sx q[2];
rz(-2.2377917) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7865078) q[1];
sx q[1];
rz(-1.4889281) q[1];
sx q[1];
rz(-2.2360389) q[1];
x q[2];
rz(-0.61935987) q[3];
sx q[3];
rz(-0.97638789) q[3];
sx q[3];
rz(-2.847282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6881037) q[2];
sx q[2];
rz(-2.0192396) q[2];
sx q[2];
rz(-2.2171059) q[2];
rz(1.2201355) q[3];
sx q[3];
rz(-2.852738) q[3];
sx q[3];
rz(0.060062241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7500551) q[0];
sx q[0];
rz(-2.4704762) q[0];
sx q[0];
rz(1.3439939) q[0];
rz(-1.9186107) q[1];
sx q[1];
rz(-1.5822625) q[1];
sx q[1];
rz(1.5498243) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5179199) q[0];
sx q[0];
rz(-2.8123283) q[0];
sx q[0];
rz(-0.98401208) q[0];
rz(0.37495592) q[2];
sx q[2];
rz(-1.4658615) q[2];
sx q[2];
rz(-2.8255759) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.87639522) q[1];
sx q[1];
rz(-1.9361456) q[1];
sx q[1];
rz(-0.72138924) q[1];
rz(-pi) q[2];
rz(0.76390169) q[3];
sx q[3];
rz(-1.502486) q[3];
sx q[3];
rz(0.98136434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.082077114) q[2];
sx q[2];
rz(-2.1166708) q[2];
sx q[2];
rz(-0.22641851) q[2];
rz(0.039479937) q[3];
sx q[3];
rz(-1.4791146) q[3];
sx q[3];
rz(1.1994908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74462849) q[0];
sx q[0];
rz(-1.9292984) q[0];
sx q[0];
rz(0.82897559) q[0];
rz(0.12116155) q[1];
sx q[1];
rz(-2.1794901) q[1];
sx q[1];
rz(1.4519579) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71428107) q[0];
sx q[0];
rz(-1.4330919) q[0];
sx q[0];
rz(0.72231035) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0214861) q[2];
sx q[2];
rz(-2.3419437) q[2];
sx q[2];
rz(0.19935184) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.87178225) q[1];
sx q[1];
rz(-2.2024558) q[1];
sx q[1];
rz(-0.58806432) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97254689) q[3];
sx q[3];
rz(-2.498811) q[3];
sx q[3];
rz(-1.9721667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5539603) q[2];
sx q[2];
rz(-0.82531896) q[2];
sx q[2];
rz(2.2770605) q[2];
rz(0.65166059) q[3];
sx q[3];
rz(-1.0385907) q[3];
sx q[3];
rz(-0.017875044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.637735) q[0];
sx q[0];
rz(-2.2561769) q[0];
sx q[0];
rz(-0.46022415) q[0];
rz(-1.3453311) q[1];
sx q[1];
rz(-0.63667744) q[1];
sx q[1];
rz(1.3705137) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11281989) q[0];
sx q[0];
rz(-1.1826853) q[0];
sx q[0];
rz(2.81243) q[0];
x q[1];
rz(2.7659594) q[2];
sx q[2];
rz(-2.4690095) q[2];
sx q[2];
rz(3.1194558) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.94356905) q[1];
sx q[1];
rz(-2.7881099) q[1];
sx q[1];
rz(2.3193053) q[1];
rz(-pi) q[2];
rz(1.8200582) q[3];
sx q[3];
rz(-2.9279997) q[3];
sx q[3];
rz(0.64208618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3687849) q[2];
sx q[2];
rz(-1.3489172) q[2];
sx q[2];
rz(0.14300145) q[2];
rz(0.16608206) q[3];
sx q[3];
rz(-0.69286418) q[3];
sx q[3];
rz(2.3486229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3731257) q[0];
sx q[0];
rz(-1.65092) q[0];
sx q[0];
rz(-1.7478818) q[0];
rz(1.1557747) q[1];
sx q[1];
rz(-1.118569) q[1];
sx q[1];
rz(-2.7973693) q[1];
rz(2.406698) q[2];
sx q[2];
rz(-1.1558888) q[2];
sx q[2];
rz(0.75239858) q[2];
rz(2.803346) q[3];
sx q[3];
rz(-2.1460642) q[3];
sx q[3];
rz(1.2301302) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
