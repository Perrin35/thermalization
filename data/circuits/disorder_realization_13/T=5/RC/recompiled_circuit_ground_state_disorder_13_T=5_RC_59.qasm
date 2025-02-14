OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.00081113022) q[0];
sx q[0];
rz(2.3164764) q[0];
sx q[0];
rz(9.0029132) q[0];
rz(2.3946664) q[1];
sx q[1];
rz(-0.59023017) q[1];
sx q[1];
rz(0.8210558) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29971545) q[0];
sx q[0];
rz(-1.4869191) q[0];
sx q[0];
rz(-1.5257349) q[0];
rz(0.73696359) q[2];
sx q[2];
rz(-0.30084494) q[2];
sx q[2];
rz(2.414051) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6118879) q[1];
sx q[1];
rz(-2.6571353) q[1];
sx q[1];
rz(2.2147708) q[1];
rz(-pi) q[2];
rz(1.2116634) q[3];
sx q[3];
rz(-0.36595038) q[3];
sx q[3];
rz(-2.8173878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7464298) q[2];
sx q[2];
rz(-0.91150993) q[2];
sx q[2];
rz(0.54473031) q[2];
rz(1.644545) q[3];
sx q[3];
rz(-2.0293197) q[3];
sx q[3];
rz(1.3397608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4561975) q[0];
sx q[0];
rz(-2.39769) q[0];
sx q[0];
rz(2.5681382) q[0];
rz(0.042512976) q[1];
sx q[1];
rz(-1.7492234) q[1];
sx q[1];
rz(-1.858985) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61717466) q[0];
sx q[0];
rz(-1.5502009) q[0];
sx q[0];
rz(-3.0153946) q[0];
rz(-pi) q[1];
rz(-1.8268711) q[2];
sx q[2];
rz(-1.8559858) q[2];
sx q[2];
rz(1.8023173) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.01067082) q[1];
sx q[1];
rz(-1.0464484) q[1];
sx q[1];
rz(-1.7337383) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76000472) q[3];
sx q[3];
rz(-2.3541321) q[3];
sx q[3];
rz(-1.2889047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.87806988) q[2];
sx q[2];
rz(-1.1493378) q[2];
sx q[2];
rz(2.6980706) q[2];
rz(-1.589132) q[3];
sx q[3];
rz(-0.3905206) q[3];
sx q[3];
rz(-2.922557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012101128) q[0];
sx q[0];
rz(-2.2624367) q[0];
sx q[0];
rz(2.7398859) q[0];
rz(1.7237639) q[1];
sx q[1];
rz(-2.6938853) q[1];
sx q[1];
rz(-1.1161425) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6594769) q[0];
sx q[0];
rz(-1.2667313) q[0];
sx q[0];
rz(2.461654) q[0];
x q[1];
rz(-1.7557451) q[2];
sx q[2];
rz(-0.70025899) q[2];
sx q[2];
rz(-1.8691693) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7771866) q[1];
sx q[1];
rz(-2.0677797) q[1];
sx q[1];
rz(-2.6650409) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0175727) q[3];
sx q[3];
rz(-0.55202602) q[3];
sx q[3];
rz(3.0025512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.22905722) q[2];
sx q[2];
rz(-2.8450862) q[2];
sx q[2];
rz(2.1336446) q[2];
rz(-1.6352765) q[3];
sx q[3];
rz(-1.1791891) q[3];
sx q[3];
rz(-0.87872163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85663831) q[0];
sx q[0];
rz(-3.1143739) q[0];
sx q[0];
rz(2.304049) q[0];
rz(1.2614999) q[1];
sx q[1];
rz(-1.3976169) q[1];
sx q[1];
rz(-1.8998442) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.094114) q[0];
sx q[0];
rz(-2.6015208) q[0];
sx q[0];
rz(0.54757) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9622039) q[2];
sx q[2];
rz(-1.3957011) q[2];
sx q[2];
rz(0.97543699) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.95835244) q[1];
sx q[1];
rz(-2.0236402) q[1];
sx q[1];
rz(2.0264105) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7697428) q[3];
sx q[3];
rz(-1.0222553) q[3];
sx q[3];
rz(-1.2814825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9813098) q[2];
sx q[2];
rz(-0.95429388) q[2];
sx q[2];
rz(2.1295638) q[2];
rz(-2.1435598) q[3];
sx q[3];
rz(-2.4055552) q[3];
sx q[3];
rz(1.4724822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-3.0618458) q[0];
sx q[0];
rz(-2.340402) q[0];
sx q[0];
rz(2.8629942) q[0];
rz(-2.8246236) q[1];
sx q[1];
rz(-1.7520889) q[1];
sx q[1];
rz(0.46151361) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7636895) q[0];
sx q[0];
rz(-2.9660241) q[0];
sx q[0];
rz(1.2922835) q[0];
rz(-pi) q[1];
rz(-0.95038173) q[2];
sx q[2];
rz(-1.2290406) q[2];
sx q[2];
rz(-0.9859964) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.52109226) q[1];
sx q[1];
rz(-2.8847155) q[1];
sx q[1];
rz(1.9476554) q[1];
x q[2];
rz(0.096624537) q[3];
sx q[3];
rz(-2.1139675) q[3];
sx q[3];
rz(-2.4718049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3658112) q[2];
sx q[2];
rz(-1.8041939) q[2];
sx q[2];
rz(-2.865045) q[2];
rz(-2.5413399) q[3];
sx q[3];
rz(-2.4252031) q[3];
sx q[3];
rz(0.53441179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9944331) q[0];
sx q[0];
rz(-2.5614547) q[0];
sx q[0];
rz(2.4603727) q[0];
rz(2.6874806) q[1];
sx q[1];
rz(-1.157016) q[1];
sx q[1];
rz(0.62201321) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9227064) q[0];
sx q[0];
rz(-1.4758631) q[0];
sx q[0];
rz(-0.19849284) q[0];
rz(-2.1900643) q[2];
sx q[2];
rz(-1.027809) q[2];
sx q[2];
rz(-1.1752626) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.71341638) q[1];
sx q[1];
rz(-0.82943664) q[1];
sx q[1];
rz(-2.1362638) q[1];
x q[2];
rz(1.9196904) q[3];
sx q[3];
rz(-1.8233646) q[3];
sx q[3];
rz(-0.080772922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5155718) q[2];
sx q[2];
rz(-0.61352789) q[2];
sx q[2];
rz(0.75508368) q[2];
rz(3.0355022) q[3];
sx q[3];
rz(-1.875149) q[3];
sx q[3];
rz(2.6809926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0945011) q[0];
sx q[0];
rz(-0.61151183) q[0];
sx q[0];
rz(3.0754572) q[0];
rz(-3.0905837) q[1];
sx q[1];
rz(-2.3936733) q[1];
sx q[1];
rz(1.8775108) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.105427) q[0];
sx q[0];
rz(-1.9308865) q[0];
sx q[0];
rz(-2.197916) q[0];
x q[1];
rz(-1.8813921) q[2];
sx q[2];
rz(-2.2618838) q[2];
sx q[2];
rz(-2.4533009) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7552774) q[1];
sx q[1];
rz(-1.1799954) q[1];
sx q[1];
rz(-0.74797191) q[1];
rz(-pi) q[2];
rz(2.407526) q[3];
sx q[3];
rz(-2.5841004) q[3];
sx q[3];
rz(-2.8898784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.85822004) q[2];
sx q[2];
rz(-1.1274575) q[2];
sx q[2];
rz(-0.14632012) q[2];
rz(-0.60394168) q[3];
sx q[3];
rz(-2.5950409) q[3];
sx q[3];
rz(1.7291791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42501763) q[0];
sx q[0];
rz(-1.1973493) q[0];
sx q[0];
rz(0.091751598) q[0];
rz(0.07235202) q[1];
sx q[1];
rz(-1.3183343) q[1];
sx q[1];
rz(2.4718557) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7726752) q[0];
sx q[0];
rz(-1.2784504) q[0];
sx q[0];
rz(1.7122853) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7013427) q[2];
sx q[2];
rz(-1.9770844) q[2];
sx q[2];
rz(-2.8216336) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.32423985) q[1];
sx q[1];
rz(-1.2659371) q[1];
sx q[1];
rz(-1.7225811) q[1];
rz(-pi) q[2];
rz(2.0132695) q[3];
sx q[3];
rz(-1.7024346) q[3];
sx q[3];
rz(1.7385755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.79121315) q[2];
sx q[2];
rz(-2.4243441) q[2];
sx q[2];
rz(-2.4199602) q[2];
rz(2.3676938) q[3];
sx q[3];
rz(-0.98267233) q[3];
sx q[3];
rz(-0.53988808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4197107) q[0];
sx q[0];
rz(-1.9419436) q[0];
sx q[0];
rz(-0.53959674) q[0];
rz(0.78060141) q[1];
sx q[1];
rz(-1.8949948) q[1];
sx q[1];
rz(-2.7194729) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1273859) q[0];
sx q[0];
rz(-0.84833586) q[0];
sx q[0];
rz(2.3425472) q[0];
rz(-1.5109748) q[2];
sx q[2];
rz(-2.0622396) q[2];
sx q[2];
rz(-0.37183842) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4735225) q[1];
sx q[1];
rz(-2.0501185) q[1];
sx q[1];
rz(0.042906656) q[1];
x q[2];
rz(-0.067764564) q[3];
sx q[3];
rz(-2.5518806) q[3];
sx q[3];
rz(-2.7763491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.098794) q[2];
sx q[2];
rz(-2.1197987) q[2];
sx q[2];
rz(0.54879028) q[2];
rz(0.15268606) q[3];
sx q[3];
rz(-2.1137674) q[3];
sx q[3];
rz(-2.4299183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.3720836) q[0];
sx q[0];
rz(-2.7476269) q[0];
sx q[0];
rz(0.60419303) q[0];
rz(-1.6744042) q[1];
sx q[1];
rz(-1.7908275) q[1];
sx q[1];
rz(1.9942572) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13191053) q[0];
sx q[0];
rz(-1.9846801) q[0];
sx q[0];
rz(2.4538246) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4060535) q[2];
sx q[2];
rz(-1.333964) q[2];
sx q[2];
rz(2.9273667) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5512451) q[1];
sx q[1];
rz(-0.4830803) q[1];
sx q[1];
rz(-0.03309588) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2828034) q[3];
sx q[3];
rz(-0.50903532) q[3];
sx q[3];
rz(-1.3648175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4597822) q[2];
sx q[2];
rz(-1.9777538) q[2];
sx q[2];
rz(2.5414844) q[2];
rz(2.4188304) q[3];
sx q[3];
rz(-2.1434651) q[3];
sx q[3];
rz(0.37989894) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73846524) q[0];
sx q[0];
rz(-1.5140139) q[0];
sx q[0];
rz(-1.4366666) q[0];
rz(-1.5557095) q[1];
sx q[1];
rz(-1.3974421) q[1];
sx q[1];
rz(-0.93346649) q[1];
rz(-0.49028291) q[2];
sx q[2];
rz(-0.96302196) q[2];
sx q[2];
rz(-1.8204126) q[2];
rz(0.71674552) q[3];
sx q[3];
rz(-2.8195753) q[3];
sx q[3];
rz(-1.4618518) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
