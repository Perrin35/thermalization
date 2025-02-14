OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.719425) q[0];
sx q[0];
rz(-0.72501215) q[0];
sx q[0];
rz(-2.3056735) q[0];
rz(-2.3140276) q[1];
sx q[1];
rz(-0.27048549) q[1];
sx q[1];
rz(-2.7343813) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0969047) q[0];
sx q[0];
rz(-0.36806199) q[0];
sx q[0];
rz(-0.9442455) q[0];
x q[1];
rz(0.23925067) q[2];
sx q[2];
rz(-1.9209338) q[2];
sx q[2];
rz(0.020933271) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5558124) q[1];
sx q[1];
rz(-1.7051053) q[1];
sx q[1];
rz(-2.7463169) q[1];
rz(-pi) q[2];
rz(0.61032774) q[3];
sx q[3];
rz(-0.84926987) q[3];
sx q[3];
rz(1.5278097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3621346) q[2];
sx q[2];
rz(-0.3388277) q[2];
sx q[2];
rz(0.96620488) q[2];
rz(-2.2400014) q[3];
sx q[3];
rz(-0.85351557) q[3];
sx q[3];
rz(2.0961659) q[3];
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
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40039429) q[0];
sx q[0];
rz(-0.61120954) q[0];
sx q[0];
rz(-3.0389767) q[0];
rz(-1.1688894) q[1];
sx q[1];
rz(-2.1714307) q[1];
sx q[1];
rz(2.6279367) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4054905) q[0];
sx q[0];
rz(-1.0879192) q[0];
sx q[0];
rz(-3.0391284) q[0];
x q[1];
rz(0.53698009) q[2];
sx q[2];
rz(-1.4479847) q[2];
sx q[2];
rz(1.7857743) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2741435) q[1];
sx q[1];
rz(-1.2701057) q[1];
sx q[1];
rz(-2.1656028) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30269514) q[3];
sx q[3];
rz(-1.154457) q[3];
sx q[3];
rz(-2.8973605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1932842) q[2];
sx q[2];
rz(-0.19364348) q[2];
sx q[2];
rz(-0.2085169) q[2];
rz(0.68566132) q[3];
sx q[3];
rz(-1.0582358) q[3];
sx q[3];
rz(0.16511495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0971646) q[0];
sx q[0];
rz(-1.2760289) q[0];
sx q[0];
rz(1.9097419) q[0];
rz(2.7992898) q[1];
sx q[1];
rz(-1.1016223) q[1];
sx q[1];
rz(1.5922348) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4195815) q[0];
sx q[0];
rz(-1.3886459) q[0];
sx q[0];
rz(0.31690449) q[0];
rz(-pi) q[1];
rz(1.4502656) q[2];
sx q[2];
rz(-1.4107586) q[2];
sx q[2];
rz(-1.1628422) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.07330896) q[1];
sx q[1];
rz(-1.7080293) q[1];
sx q[1];
rz(-0.84899606) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3089448) q[3];
sx q[3];
rz(-0.73048985) q[3];
sx q[3];
rz(-0.60527847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0478829) q[2];
sx q[2];
rz(-1.9423395) q[2];
sx q[2];
rz(-1.017978) q[2];
rz(-1.4114981) q[3];
sx q[3];
rz(-2.394815) q[3];
sx q[3];
rz(1.5311034) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8133076) q[0];
sx q[0];
rz(-1.7787373) q[0];
sx q[0];
rz(3.0174729) q[0];
rz(2.2165551) q[1];
sx q[1];
rz(-1.8412453) q[1];
sx q[1];
rz(-0.73840028) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3401854) q[0];
sx q[0];
rz(-0.75447318) q[0];
sx q[0];
rz(1.6615926) q[0];
rz(2.5321711) q[2];
sx q[2];
rz(-1.2487234) q[2];
sx q[2];
rz(-0.24064482) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.70917801) q[1];
sx q[1];
rz(-1.3644618) q[1];
sx q[1];
rz(-0.56992759) q[1];
x q[2];
rz(2.9876338) q[3];
sx q[3];
rz(-2.5601644) q[3];
sx q[3];
rz(-2.1179782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.063252123) q[2];
sx q[2];
rz(-1.1880778) q[2];
sx q[2];
rz(1.1032907) q[2];
rz(2.8830146) q[3];
sx q[3];
rz(-2.0662112) q[3];
sx q[3];
rz(-2.8598089) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2461808) q[0];
sx q[0];
rz(-0.87765944) q[0];
sx q[0];
rz(1.3377162) q[0];
rz(2.5216263) q[1];
sx q[1];
rz(-1.4620616) q[1];
sx q[1];
rz(1.5043129) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97457394) q[0];
sx q[0];
rz(-1.3822379) q[0];
sx q[0];
rz(1.2861286) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4392534) q[2];
sx q[2];
rz(-1.6026671) q[2];
sx q[2];
rz(-1.7447217) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6474939) q[1];
sx q[1];
rz(-1.955619) q[1];
sx q[1];
rz(0.61487756) q[1];
rz(-1.6880468) q[3];
sx q[3];
rz(-2.4053889) q[3];
sx q[3];
rz(0.29090009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6164246) q[2];
sx q[2];
rz(-1.1526356) q[2];
sx q[2];
rz(-1.1241414) q[2];
rz(-0.21640402) q[3];
sx q[3];
rz(-2.3085322) q[3];
sx q[3];
rz(-1.0696577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6935317) q[0];
sx q[0];
rz(-0.8140642) q[0];
sx q[0];
rz(2.6779209) q[0];
rz(-3.0060153) q[1];
sx q[1];
rz(-0.99628535) q[1];
sx q[1];
rz(-0.3840951) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039003619) q[0];
sx q[0];
rz(-1.4560934) q[0];
sx q[0];
rz(2.6713598) q[0];
x q[1];
rz(-0.21125085) q[2];
sx q[2];
rz(-2.3659035) q[2];
sx q[2];
rz(-0.83310177) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.64318618) q[1];
sx q[1];
rz(-1.2075042) q[1];
sx q[1];
rz(0.17885991) q[1];
rz(1.4469524) q[3];
sx q[3];
rz(-0.55168024) q[3];
sx q[3];
rz(1.1744896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42528459) q[2];
sx q[2];
rz(-0.3766489) q[2];
sx q[2];
rz(0.46871218) q[2];
rz(-2.5675755) q[3];
sx q[3];
rz(-1.470397) q[3];
sx q[3];
rz(0.65042692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18844093) q[0];
sx q[0];
rz(-1.3583536) q[0];
sx q[0];
rz(2.4515732) q[0];
rz(-0.72235876) q[1];
sx q[1];
rz(-1.1411618) q[1];
sx q[1];
rz(-0.77675995) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39333568) q[0];
sx q[0];
rz(-1.3384755) q[0];
sx q[0];
rz(1.6477702) q[0];
rz(0.58906196) q[2];
sx q[2];
rz(-1.2576332) q[2];
sx q[2];
rz(-0.89436603) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1550786) q[1];
sx q[1];
rz(-2.6342794) q[1];
sx q[1];
rz(0.55965565) q[1];
x q[2];
rz(-1.5942105) q[3];
sx q[3];
rz(-2.4613614) q[3];
sx q[3];
rz(-1.3772723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7338099) q[2];
sx q[2];
rz(-1.5747384) q[2];
sx q[2];
rz(0.053827914) q[2];
rz(-2.5194061) q[3];
sx q[3];
rz(-2.6572808) q[3];
sx q[3];
rz(-3.0622862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4902041) q[0];
sx q[0];
rz(-1.1505928) q[0];
sx q[0];
rz(1.7405317) q[0];
rz(0.98848629) q[1];
sx q[1];
rz(-1.8732312) q[1];
sx q[1];
rz(0.40774694) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0893129) q[0];
sx q[0];
rz(-1.3527217) q[0];
sx q[0];
rz(-0.41357354) q[0];
rz(-pi) q[1];
rz(-1.4125702) q[2];
sx q[2];
rz(-2.4605949) q[2];
sx q[2];
rz(-2.3883377) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0474075) q[1];
sx q[1];
rz(-2.2360206) q[1];
sx q[1];
rz(-2.1783504) q[1];
rz(-1.4325822) q[3];
sx q[3];
rz(-1.301389) q[3];
sx q[3];
rz(0.22002735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.72655788) q[2];
sx q[2];
rz(-1.1316391) q[2];
sx q[2];
rz(-2.4416907) q[2];
rz(-1.0837726) q[3];
sx q[3];
rz(-0.86728573) q[3];
sx q[3];
rz(-0.23610246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20339762) q[0];
sx q[0];
rz(-1.45881) q[0];
sx q[0];
rz(-2.7517125) q[0];
rz(-1.9436721) q[1];
sx q[1];
rz(-0.094466297) q[1];
sx q[1];
rz(-0.87497154) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1221473) q[0];
sx q[0];
rz(-2.5179389) q[0];
sx q[0];
rz(1.532293) q[0];
x q[1];
rz(1.0880555) q[2];
sx q[2];
rz(-0.74626479) q[2];
sx q[2];
rz(-1.3811228) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.81351133) q[1];
sx q[1];
rz(-1.2248779) q[1];
sx q[1];
rz(2.6807869) q[1];
x q[2];
rz(1.0040818) q[3];
sx q[3];
rz(-1.9329786) q[3];
sx q[3];
rz(0.099206533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3206869) q[2];
sx q[2];
rz(-2.4497538) q[2];
sx q[2];
rz(-0.64986491) q[2];
rz(-1.1167022) q[3];
sx q[3];
rz(-1.8442804) q[3];
sx q[3];
rz(-2.9445649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4138625) q[0];
sx q[0];
rz(-3.08857) q[0];
sx q[0];
rz(-2.9470288) q[0];
rz(-0.77231705) q[1];
sx q[1];
rz(-1.1338502) q[1];
sx q[1];
rz(2.9639941) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0227764) q[0];
sx q[0];
rz(-1.5714065) q[0];
sx q[0];
rz(1.5724036) q[0];
rz(-0.9382605) q[2];
sx q[2];
rz(-1.0192945) q[2];
sx q[2];
rz(0.86966627) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1298123) q[1];
sx q[1];
rz(-2.3200071) q[1];
sx q[1];
rz(2.473144) q[1];
rz(2.4836618) q[3];
sx q[3];
rz(-2.2773491) q[3];
sx q[3];
rz(-2.547675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.78581587) q[2];
sx q[2];
rz(-2.5444701) q[2];
sx q[2];
rz(-0.43006483) q[2];
rz(-1.1108584) q[3];
sx q[3];
rz(-1.6836124) q[3];
sx q[3];
rz(0.41657579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72660245) q[0];
sx q[0];
rz(-1.5530598) q[0];
sx q[0];
rz(-0.1834827) q[0];
rz(1.1732187) q[1];
sx q[1];
rz(-1.9216187) q[1];
sx q[1];
rz(0.12513195) q[1];
rz(2.0424615) q[2];
sx q[2];
rz(-0.57715125) q[2];
sx q[2];
rz(-1.1751529) q[2];
rz(-2.597025) q[3];
sx q[3];
rz(-0.93251317) q[3];
sx q[3];
rz(2.9208356) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
