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
rz(-0.79987502) q[0];
sx q[0];
rz(-0.81698155) q[0];
sx q[0];
rz(2.6843827) q[0];
rz(-2.7297821) q[1];
sx q[1];
rz(-1.6219985) q[1];
sx q[1];
rz(2.5817459) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8053539) q[0];
sx q[0];
rz(-1.8913664) q[0];
sx q[0];
rz(-2.6692315) q[0];
rz(-pi) q[1];
rz(1.0798321) q[2];
sx q[2];
rz(-0.55779558) q[2];
sx q[2];
rz(-1.3809134) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.35750439) q[1];
sx q[1];
rz(-0.58114806) q[1];
sx q[1];
rz(0.52304348) q[1];
rz(-pi) q[2];
rz(-0.58497854) q[3];
sx q[3];
rz(-0.26250678) q[3];
sx q[3];
rz(2.7141904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.37602279) q[2];
sx q[2];
rz(-0.96161181) q[2];
sx q[2];
rz(1.8454856) q[2];
rz(2.5206595) q[3];
sx q[3];
rz(-0.66578484) q[3];
sx q[3];
rz(0.92500979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45944443) q[0];
sx q[0];
rz(-2.5542673) q[0];
sx q[0];
rz(0.71429724) q[0];
rz(0.722305) q[1];
sx q[1];
rz(-2.0562833) q[1];
sx q[1];
rz(-1.3246271) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2091529) q[0];
sx q[0];
rz(-2.5896833) q[0];
sx q[0];
rz(0.96262424) q[0];
rz(0.06109625) q[2];
sx q[2];
rz(-0.73373605) q[2];
sx q[2];
rz(0.82759418) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2683619) q[1];
sx q[1];
rz(-1.9385846) q[1];
sx q[1];
rz(1.7413543) q[1];
rz(0.12217317) q[3];
sx q[3];
rz(-2.2524912) q[3];
sx q[3];
rz(0.091191779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.033919949) q[2];
sx q[2];
rz(-1.4227957) q[2];
sx q[2];
rz(0.99204341) q[2];
rz(-2.3658559) q[3];
sx q[3];
rz(-0.78191596) q[3];
sx q[3];
rz(-1.0824664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7139605) q[0];
sx q[0];
rz(-1.3466703) q[0];
sx q[0];
rz(0.91208518) q[0];
rz(-2.7569547) q[1];
sx q[1];
rz(-2.1802528) q[1];
sx q[1];
rz(-2.6461163) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3867823) q[0];
sx q[0];
rz(-1.4982002) q[0];
sx q[0];
rz(3.0955546) q[0];
rz(-pi) q[1];
x q[1];
rz(2.936061) q[2];
sx q[2];
rz(-0.759911) q[2];
sx q[2];
rz(2.2874557) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.7621618) q[1];
sx q[1];
rz(-0.84745211) q[1];
sx q[1];
rz(-1.6382193) q[1];
rz(-0.94888249) q[3];
sx q[3];
rz(-1.5601741) q[3];
sx q[3];
rz(1.1158021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2428525) q[2];
sx q[2];
rz(-1.1246559) q[2];
sx q[2];
rz(0.43453547) q[2];
rz(-0.09856002) q[3];
sx q[3];
rz(-1.7193272) q[3];
sx q[3];
rz(-2.3677473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.3944405) q[0];
sx q[0];
rz(-0.23500615) q[0];
sx q[0];
rz(2.8727942) q[0];
rz(-2.0647743) q[1];
sx q[1];
rz(-1.4067168) q[1];
sx q[1];
rz(1.1928308) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6679316) q[0];
sx q[0];
rz(-2.8369378) q[0];
sx q[0];
rz(-0.61078914) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0054686) q[2];
sx q[2];
rz(-1.4413712) q[2];
sx q[2];
rz(0.65107513) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7320648) q[1];
sx q[1];
rz(-1.7842147) q[1];
sx q[1];
rz(1.7312538) q[1];
rz(-pi) q[2];
rz(-0.70130879) q[3];
sx q[3];
rz(-1.4150672) q[3];
sx q[3];
rz(0.45022545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8863135) q[2];
sx q[2];
rz(-0.89526075) q[2];
sx q[2];
rz(0.2992343) q[2];
rz(0.29087654) q[3];
sx q[3];
rz(-1.8894922) q[3];
sx q[3];
rz(0.65557426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3000779) q[0];
sx q[0];
rz(-0.97989196) q[0];
sx q[0];
rz(3.063524) q[0];
rz(0.95371753) q[1];
sx q[1];
rz(-0.51906145) q[1];
sx q[1];
rz(1.5740707) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6758976) q[0];
sx q[0];
rz(-0.90238304) q[0];
sx q[0];
rz(-0.23082478) q[0];
rz(-pi) q[1];
rz(0.6965397) q[2];
sx q[2];
rz(-1.3020421) q[2];
sx q[2];
rz(-2.0292676) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.138254) q[1];
sx q[1];
rz(-1.5383771) q[1];
sx q[1];
rz(-0.95862548) q[1];
rz(-pi) q[2];
rz(0.6930954) q[3];
sx q[3];
rz(-1.6247107) q[3];
sx q[3];
rz(-1.113232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8629525) q[2];
sx q[2];
rz(-0.43206698) q[2];
sx q[2];
rz(0.26818177) q[2];
rz(2.6523318) q[3];
sx q[3];
rz(-2.0475976) q[3];
sx q[3];
rz(1.6930273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0443403) q[0];
sx q[0];
rz(-3.1179929) q[0];
sx q[0];
rz(-2.6771255) q[0];
rz(-0.52344549) q[1];
sx q[1];
rz(-0.76952666) q[1];
sx q[1];
rz(-2.4109667) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3364612) q[0];
sx q[0];
rz(-1.4298273) q[0];
sx q[0];
rz(-0.8223429) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4695005) q[2];
sx q[2];
rz(-1.8037829) q[2];
sx q[2];
rz(1.9547878) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0981104) q[1];
sx q[1];
rz(-0.93096369) q[1];
sx q[1];
rz(-2.1049064) q[1];
rz(-pi) q[2];
rz(0.2517638) q[3];
sx q[3];
rz(-1.442896) q[3];
sx q[3];
rz(2.3826016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0940493) q[2];
sx q[2];
rz(-1.0796248) q[2];
sx q[2];
rz(-0.13761061) q[2];
rz(2.008647) q[3];
sx q[3];
rz(-1.9652941) q[3];
sx q[3];
rz(-1.2631811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9889744) q[0];
sx q[0];
rz(-1.512383) q[0];
sx q[0];
rz(0.033893943) q[0];
rz(-2.7821817) q[1];
sx q[1];
rz(-1.6172599) q[1];
sx q[1];
rz(-0.83438897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2914258) q[0];
sx q[0];
rz(-1.8803673) q[0];
sx q[0];
rz(-0.39637027) q[0];
x q[1];
rz(2.4530386) q[2];
sx q[2];
rz(-1.1821185) q[2];
sx q[2];
rz(1.9118903) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.31454489) q[1];
sx q[1];
rz(-0.27092182) q[1];
sx q[1];
rz(1.2071868) q[1];
x q[2];
rz(-3.062617) q[3];
sx q[3];
rz(-2.4684836) q[3];
sx q[3];
rz(2.1977294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.72839165) q[2];
sx q[2];
rz(-0.30343702) q[2];
sx q[2];
rz(-2.0835853) q[2];
rz(-0.47438619) q[3];
sx q[3];
rz(-2.3460903) q[3];
sx q[3];
rz(2.0856196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0742842) q[0];
sx q[0];
rz(-1.5165167) q[0];
sx q[0];
rz(1.3577331) q[0];
rz(-2.7108497) q[1];
sx q[1];
rz(-1.6669225) q[1];
sx q[1];
rz(0.46708435) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8170094) q[0];
sx q[0];
rz(-1.6139493) q[0];
sx q[0];
rz(-0.0014716455) q[0];
rz(0.70286669) q[2];
sx q[2];
rz(-1.4231668) q[2];
sx q[2];
rz(-0.63181782) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.39462806) q[1];
sx q[1];
rz(-2.10022) q[1];
sx q[1];
rz(1.5922597) q[1];
rz(-2.2086618) q[3];
sx q[3];
rz(-1.5278421) q[3];
sx q[3];
rz(1.651498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0908541) q[2];
sx q[2];
rz(-0.41768062) q[2];
sx q[2];
rz(2.0609071) q[2];
rz(-0.68197) q[3];
sx q[3];
rz(-2.3365648) q[3];
sx q[3];
rz(-2.4431156) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68898106) q[0];
sx q[0];
rz(-0.51849759) q[0];
sx q[0];
rz(-2.3338351) q[0];
rz(1.0001596) q[1];
sx q[1];
rz(-2.2554485) q[1];
sx q[1];
rz(2.3449786) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0611872) q[0];
sx q[0];
rz(-0.10594254) q[0];
sx q[0];
rz(1.9796014) q[0];
rz(-pi) q[1];
rz(0.32795017) q[2];
sx q[2];
rz(-2.4900576) q[2];
sx q[2];
rz(0.33199379) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.63846) q[1];
sx q[1];
rz(-1.1640932) q[1];
sx q[1];
rz(-0.15127123) q[1];
rz(-1.8851938) q[3];
sx q[3];
rz(-1.8125497) q[3];
sx q[3];
rz(-0.35255656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.54958582) q[2];
sx q[2];
rz(-2.0506115) q[2];
sx q[2];
rz(-2.8263212) q[2];
rz(-1.5677876) q[3];
sx q[3];
rz(-1.1476293) q[3];
sx q[3];
rz(1.1228336) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0386117) q[0];
sx q[0];
rz(-0.6655612) q[0];
sx q[0];
rz(-2.3200206) q[0];
rz(2.6577677) q[1];
sx q[1];
rz(-1.7211823) q[1];
sx q[1];
rz(0.22463591) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0517462) q[0];
sx q[0];
rz(-1.4205853) q[0];
sx q[0];
rz(1.2443022) q[0];
x q[1];
rz(-2.2251525) q[2];
sx q[2];
rz(-2.096837) q[2];
sx q[2];
rz(-1.9182084) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37525081) q[1];
sx q[1];
rz(-0.90000376) q[1];
sx q[1];
rz(-2.8743582) q[1];
rz(-pi) q[2];
rz(-1.3566689) q[3];
sx q[3];
rz(-1.6928634) q[3];
sx q[3];
rz(1.7904953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.75767526) q[2];
sx q[2];
rz(-3.0249247) q[2];
sx q[2];
rz(-3.0465916) q[2];
rz(1.654024) q[3];
sx q[3];
rz(-0.58949685) q[3];
sx q[3];
rz(-0.76475638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8944396) q[0];
sx q[0];
rz(-2.7338487) q[0];
sx q[0];
rz(-0.26640531) q[0];
rz(-15/(8*pi)) q[1];
sx q[1];
rz(-1.4529556) q[1];
sx q[1];
rz(-1.6642889) q[1];
rz(-2.2483027) q[2];
sx q[2];
rz(-1.2883452) q[2];
sx q[2];
rz(2.7569994) q[2];
rz(1.5726907) q[3];
sx q[3];
rz(-1.9463149) q[3];
sx q[3];
rz(-2.7494242) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
