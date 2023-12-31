OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0796354) q[0];
sx q[0];
rz(-2.893653) q[0];
sx q[0];
rz(0.82013446) q[0];
rz(0.040854383) q[1];
sx q[1];
rz(-2.3526469) q[1];
sx q[1];
rz(0.087892428) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.983842) q[0];
sx q[0];
rz(-2.0294667) q[0];
sx q[0];
rz(1.0755324) q[0];
rz(0.90859969) q[2];
sx q[2];
rz(-2.5052921) q[2];
sx q[2];
rz(1.5696021) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.81401285) q[1];
sx q[1];
rz(-2.6373133) q[1];
sx q[1];
rz(3.1370647) q[1];
x q[2];
rz(3.1357364) q[3];
sx q[3];
rz(-1.8294888) q[3];
sx q[3];
rz(-0.70219016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1408046) q[2];
sx q[2];
rz(-1.6892097) q[2];
sx q[2];
rz(0.56935707) q[2];
rz(1.5287483) q[3];
sx q[3];
rz(-2.5053535) q[3];
sx q[3];
rz(1.7830085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5434718) q[0];
sx q[0];
rz(-0.16508979) q[0];
sx q[0];
rz(-2.5879481) q[0];
rz(1.9042632) q[1];
sx q[1];
rz(-1.7747223) q[1];
sx q[1];
rz(-1.233261) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20576661) q[0];
sx q[0];
rz(-2.2201949) q[0];
sx q[0];
rz(-2.4539024) q[0];
rz(-pi) q[1];
rz(-2.4670062) q[2];
sx q[2];
rz(-1.470675) q[2];
sx q[2];
rz(-1.9763725) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.711449) q[1];
sx q[1];
rz(-1.5866382) q[1];
sx q[1];
rz(1.364691) q[1];
x q[2];
rz(2.4741461) q[3];
sx q[3];
rz(-1.7597886) q[3];
sx q[3];
rz(0.28039704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.0040940293) q[2];
sx q[2];
rz(-1.6823021) q[2];
sx q[2];
rz(3.1393576) q[2];
rz(-0.8301174) q[3];
sx q[3];
rz(-2.3486962) q[3];
sx q[3];
rz(1.8921651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24901351) q[0];
sx q[0];
rz(-0.61616388) q[0];
sx q[0];
rz(-2.2900443) q[0];
rz(-0.7710723) q[1];
sx q[1];
rz(-1.4665736) q[1];
sx q[1];
rz(3.070389) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.203513) q[0];
sx q[0];
rz(-2.8305614) q[0];
sx q[0];
rz(-2.2631893) q[0];
rz(-2.0729012) q[2];
sx q[2];
rz(-0.50902589) q[2];
sx q[2];
rz(0.44391649) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.60656602) q[1];
sx q[1];
rz(-2.2177794) q[1];
sx q[1];
rz(1.4024629) q[1];
rz(-pi) q[2];
rz(-0.72910587) q[3];
sx q[3];
rz(-2.2359072) q[3];
sx q[3];
rz(-0.8641181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8292134) q[2];
sx q[2];
rz(-0.92386121) q[2];
sx q[2];
rz(2.6181347) q[2];
rz(-2.8097025) q[3];
sx q[3];
rz(-3.0811946) q[3];
sx q[3];
rz(-2.6383242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4098542) q[0];
sx q[0];
rz(-2.8926352) q[0];
sx q[0];
rz(2.3160034) q[0];
rz(-1.1666974) q[1];
sx q[1];
rz(-2.2749133) q[1];
sx q[1];
rz(-1.694214) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3576413) q[0];
sx q[0];
rz(-0.28342993) q[0];
sx q[0];
rz(-0.95143239) q[0];
x q[1];
rz(-1.5201735) q[2];
sx q[2];
rz(-0.17855893) q[2];
sx q[2];
rz(-2.6176917) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.64926681) q[1];
sx q[1];
rz(-2.4267303) q[1];
sx q[1];
rz(-1.3267172) q[1];
rz(-1.3284239) q[3];
sx q[3];
rz(-1.216785) q[3];
sx q[3];
rz(-0.44204516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6772785) q[2];
sx q[2];
rz(-1.7769287) q[2];
sx q[2];
rz(0.14492598) q[2];
rz(-1.011301) q[3];
sx q[3];
rz(-2.5144808) q[3];
sx q[3];
rz(-3.1269126) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1458364) q[0];
sx q[0];
rz(-3.0881112) q[0];
sx q[0];
rz(2.4575535) q[0];
rz(1.9513291) q[1];
sx q[1];
rz(-0.66527706) q[1];
sx q[1];
rz(1.9285944) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4764403) q[0];
sx q[0];
rz(-1.9422918) q[0];
sx q[0];
rz(-1.2907589) q[0];
rz(-2.2239387) q[2];
sx q[2];
rz(-2.0954164) q[2];
sx q[2];
rz(1.4337002) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2781093) q[1];
sx q[1];
rz(-0.28243318) q[1];
sx q[1];
rz(-0.43538283) q[1];
rz(-pi) q[2];
rz(-0.9807068) q[3];
sx q[3];
rz(-2.8299035) q[3];
sx q[3];
rz(-0.89524549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6325536) q[2];
sx q[2];
rz(-1.8703439) q[2];
sx q[2];
rz(2.5115013) q[2];
rz(-0.57224327) q[3];
sx q[3];
rz(-2.4961491) q[3];
sx q[3];
rz(2.2563289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0260139) q[0];
sx q[0];
rz(-0.84340874) q[0];
sx q[0];
rz(0.28636006) q[0];
rz(-0.59965602) q[1];
sx q[1];
rz(-1.1279794) q[1];
sx q[1];
rz(0.58475959) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8557381) q[0];
sx q[0];
rz(-1.2585088) q[0];
sx q[0];
rz(-0.62494846) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5445968) q[2];
sx q[2];
rz(-0.95653406) q[2];
sx q[2];
rz(0.40528497) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6905745) q[1];
sx q[1];
rz(-1.283053) q[1];
sx q[1];
rz(-2.2698127) q[1];
rz(0.16634511) q[3];
sx q[3];
rz(-1.5208941) q[3];
sx q[3];
rz(0.14453105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.16453234) q[2];
sx q[2];
rz(-0.88447276) q[2];
sx q[2];
rz(1.5131081) q[2];
rz(-0.96308723) q[3];
sx q[3];
rz(-1.3685127) q[3];
sx q[3];
rz(0.62455463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0528089) q[0];
sx q[0];
rz(-1.4316906) q[0];
sx q[0];
rz(0.54876304) q[0];
rz(0.051963003) q[1];
sx q[1];
rz(-2.6497662) q[1];
sx q[1];
rz(2.5351977) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081189922) q[0];
sx q[0];
rz(-2.5572544) q[0];
sx q[0];
rz(1.1330182) q[0];
x q[1];
rz(-0.90784448) q[2];
sx q[2];
rz(-1.2693451) q[2];
sx q[2];
rz(-0.19336685) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.28601521) q[1];
sx q[1];
rz(-2.3708323) q[1];
sx q[1];
rz(-1.7325749) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2783979) q[3];
sx q[3];
rz(-0.51897012) q[3];
sx q[3];
rz(0.88245813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9397883) q[2];
sx q[2];
rz(-2.1494892) q[2];
sx q[2];
rz(-2.6331804) q[2];
rz(-2.0243747) q[3];
sx q[3];
rz(-0.17337392) q[3];
sx q[3];
rz(-1.4846876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8023119) q[0];
sx q[0];
rz(-2.7957714) q[0];
sx q[0];
rz(2.3876277) q[0];
rz(0.94999653) q[1];
sx q[1];
rz(-1.6597304) q[1];
sx q[1];
rz(2.9139013) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65093016) q[0];
sx q[0];
rz(-0.41914808) q[0];
sx q[0];
rz(1.7463669) q[0];
rz(-pi) q[1];
rz(-0.44616206) q[2];
sx q[2];
rz(-1.237962) q[2];
sx q[2];
rz(1.7762426) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.258192) q[1];
sx q[1];
rz(-1.0419783) q[1];
sx q[1];
rz(1.8316815) q[1];
x q[2];
rz(2.5826107) q[3];
sx q[3];
rz(-1.535927) q[3];
sx q[3];
rz(-2.1502286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0630539) q[2];
sx q[2];
rz(-0.20919122) q[2];
sx q[2];
rz(-0.83855808) q[2];
rz(-0.36455425) q[3];
sx q[3];
rz(-1.7611327) q[3];
sx q[3];
rz(-2.3013039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39695981) q[0];
sx q[0];
rz(-1.3752221) q[0];
sx q[0];
rz(-2.3378085) q[0];
rz(-1.0546168) q[1];
sx q[1];
rz(-2.5024253) q[1];
sx q[1];
rz(1.1484336) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20439416) q[0];
sx q[0];
rz(-3.0657401) q[0];
sx q[0];
rz(1.2501459) q[0];
rz(-pi) q[1];
rz(1.0461651) q[2];
sx q[2];
rz(-1.2534007) q[2];
sx q[2];
rz(1.5898926) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.25363898) q[1];
sx q[1];
rz(-1.1432853) q[1];
sx q[1];
rz(1.0484139) q[1];
x q[2];
rz(-3.0461379) q[3];
sx q[3];
rz(-2.7881552) q[3];
sx q[3];
rz(-2.8638338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.049872963) q[2];
sx q[2];
rz(-2.1240978) q[2];
sx q[2];
rz(-0.77825528) q[2];
rz(-0.19566472) q[3];
sx q[3];
rz(-1.0396495) q[3];
sx q[3];
rz(-1.0218609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8675999) q[0];
sx q[0];
rz(-0.99853981) q[0];
sx q[0];
rz(2.8883949) q[0];
rz(-0.7397488) q[1];
sx q[1];
rz(-0.7363798) q[1];
sx q[1];
rz(2.443312) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65933278) q[0];
sx q[0];
rz(-1.8730622) q[0];
sx q[0];
rz(1.4955826) q[0];
rz(0.024784879) q[2];
sx q[2];
rz(-1.4966655) q[2];
sx q[2];
rz(0.63400808) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.89430289) q[1];
sx q[1];
rz(-1.3548046) q[1];
sx q[1];
rz(-2.0386019) q[1];
rz(-pi) q[2];
rz(-1.2020338) q[3];
sx q[3];
rz(-0.59951111) q[3];
sx q[3];
rz(0.60839073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0239821) q[2];
sx q[2];
rz(-2.2587946) q[2];
sx q[2];
rz(2.3013766) q[2];
rz(-0.64030567) q[3];
sx q[3];
rz(-2.2655723) q[3];
sx q[3];
rz(-1.1673814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2587851) q[0];
sx q[0];
rz(-0.85492815) q[0];
sx q[0];
rz(1.7229765) q[0];
rz(0.029126833) q[1];
sx q[1];
rz(-0.11321414) q[1];
sx q[1];
rz(-1.3197457) q[1];
rz(-0.60224709) q[2];
sx q[2];
rz(-1.8613653) q[2];
sx q[2];
rz(-2.0801434) q[2];
rz(0.15016951) q[3];
sx q[3];
rz(-2.2554845) q[3];
sx q[3];
rz(3.1156202) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
