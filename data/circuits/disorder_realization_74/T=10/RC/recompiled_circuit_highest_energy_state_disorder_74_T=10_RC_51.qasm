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
rz(1.9407152) q[0];
sx q[0];
rz(-2.241029) q[0];
sx q[0];
rz(0.23298921) q[0];
rz(-1.5268582) q[1];
sx q[1];
rz(-2.0695217) q[1];
sx q[1];
rz(1.1195247) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6132076) q[0];
sx q[0];
rz(-0.23477708) q[0];
sx q[0];
rz(1.6260207) q[0];
rz(2.5801611) q[2];
sx q[2];
rz(-1.5806287) q[2];
sx q[2];
rz(-1.9468284) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.2762766) q[1];
sx q[1];
rz(-0.36082339) q[1];
sx q[1];
rz(1.3555384) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2337716) q[3];
sx q[3];
rz(-2.153018) q[3];
sx q[3];
rz(-1.9523417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4423674) q[2];
sx q[2];
rz(-1.0590326) q[2];
sx q[2];
rz(0.11288682) q[2];
rz(-0.26116192) q[3];
sx q[3];
rz(-1.3449202) q[3];
sx q[3];
rz(-1.2507218) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79171044) q[0];
sx q[0];
rz(-3.0630906) q[0];
sx q[0];
rz(-3.0872524) q[0];
rz(0.21866523) q[1];
sx q[1];
rz(-1.668914) q[1];
sx q[1];
rz(-0.36453077) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2549627) q[0];
sx q[0];
rz(-2.1658362) q[0];
sx q[0];
rz(0.9592077) q[0];
rz(2.0979375) q[2];
sx q[2];
rz(-0.69913617) q[2];
sx q[2];
rz(1.8086036) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6682525) q[1];
sx q[1];
rz(-2.1917731) q[1];
sx q[1];
rz(1.4771858) q[1];
x q[2];
rz(1.0899312) q[3];
sx q[3];
rz(-0.37853795) q[3];
sx q[3];
rz(-0.35096452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9008122) q[2];
sx q[2];
rz(-1.6996926) q[2];
sx q[2];
rz(2.0330632) q[2];
rz(-0.19567868) q[3];
sx q[3];
rz(-1.9340197) q[3];
sx q[3];
rz(1.6746563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3942669) q[0];
sx q[0];
rz(-1.0402004) q[0];
sx q[0];
rz(1.0362097) q[0];
rz(0.58468435) q[1];
sx q[1];
rz(-1.4671289) q[1];
sx q[1];
rz(-2.5856957) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7202458) q[0];
sx q[0];
rz(-1.9860391) q[0];
sx q[0];
rz(0.74851997) q[0];
x q[1];
rz(3.1039882) q[2];
sx q[2];
rz(-1.6480664) q[2];
sx q[2];
rz(-1.8118878) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7207009) q[1];
sx q[1];
rz(-1.4399505) q[1];
sx q[1];
rz(-2.2761005) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49212891) q[3];
sx q[3];
rz(-0.50072296) q[3];
sx q[3];
rz(-2.2496417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7613775) q[2];
sx q[2];
rz(-0.20647241) q[2];
sx q[2];
rz(-1.9492487) q[2];
rz(0.10382593) q[3];
sx q[3];
rz(-1.1622279) q[3];
sx q[3];
rz(1.1473568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1059859) q[0];
sx q[0];
rz(-2.5150531) q[0];
sx q[0];
rz(-1.4242127) q[0];
rz(0.72319889) q[1];
sx q[1];
rz(-2.9056748) q[1];
sx q[1];
rz(-2.9211488) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0940277) q[0];
sx q[0];
rz(-1.936543) q[0];
sx q[0];
rz(2.0630552) q[0];
rz(-pi) q[1];
rz(-2.0218684) q[2];
sx q[2];
rz(-0.60904087) q[2];
sx q[2];
rz(0.90897564) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6666777) q[1];
sx q[1];
rz(-1.4581469) q[1];
sx q[1];
rz(1.9289609) q[1];
rz(0.8414874) q[3];
sx q[3];
rz(-0.53181767) q[3];
sx q[3];
rz(-1.8790661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.38349884) q[2];
sx q[2];
rz(-1.8269962) q[2];
sx q[2];
rz(-2.626075) q[2];
rz(-3.0564485) q[3];
sx q[3];
rz(-2.4862423) q[3];
sx q[3];
rz(2.3084124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4715217) q[0];
sx q[0];
rz(-0.18297289) q[0];
sx q[0];
rz(2.4772189) q[0];
rz(-0.5270671) q[1];
sx q[1];
rz(-2.0030237) q[1];
sx q[1];
rz(-1.1767496) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.697061) q[0];
sx q[0];
rz(-0.089655487) q[0];
sx q[0];
rz(-1.2517002) q[0];
rz(-pi) q[1];
rz(-2.2559153) q[2];
sx q[2];
rz(-1.9373477) q[2];
sx q[2];
rz(-1.4057856) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5866111) q[1];
sx q[1];
rz(-0.23906318) q[1];
sx q[1];
rz(-0.081940941) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.142435) q[3];
sx q[3];
rz(-1.7222341) q[3];
sx q[3];
rz(-3.0157523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.027017) q[2];
sx q[2];
rz(-0.21275529) q[2];
sx q[2];
rz(2.2034755) q[2];
rz(-2.0959057) q[3];
sx q[3];
rz(-0.18200471) q[3];
sx q[3];
rz(-0.81030455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3351347) q[0];
sx q[0];
rz(-1.9430176) q[0];
sx q[0];
rz(1.8916116) q[0];
rz(0.20420034) q[1];
sx q[1];
rz(-2.4649492) q[1];
sx q[1];
rz(-0.85320371) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4049339) q[0];
sx q[0];
rz(-2.6879426) q[0];
sx q[0];
rz(-2.4429295) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3647652) q[2];
sx q[2];
rz(-2.3558801) q[2];
sx q[2];
rz(-1.7091027) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0029162) q[1];
sx q[1];
rz(-2.2670548) q[1];
sx q[1];
rz(2.584444) q[1];
rz(3.0141287) q[3];
sx q[3];
rz(-1.9098305) q[3];
sx q[3];
rz(-1.6515428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.044067232) q[2];
sx q[2];
rz(-1.7996457) q[2];
sx q[2];
rz(1.9419144) q[2];
rz(-2.1357644) q[3];
sx q[3];
rz(-2.2483716) q[3];
sx q[3];
rz(-0.83542663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47578874) q[0];
sx q[0];
rz(-2.8732193) q[0];
sx q[0];
rz(-2.1589808) q[0];
rz(-1.3081374) q[1];
sx q[1];
rz(-1.2160701) q[1];
sx q[1];
rz(1.4391724) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2932201) q[0];
sx q[0];
rz(-1.2407899) q[0];
sx q[0];
rz(-2.7820935) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.077567049) q[2];
sx q[2];
rz(-2.0915389) q[2];
sx q[2];
rz(-2.5375053) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5251617) q[1];
sx q[1];
rz(-2.3724864) q[1];
sx q[1];
rz(2.0517212) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4286707) q[3];
sx q[3];
rz(-1.097659) q[3];
sx q[3];
rz(-0.17156916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.16284379) q[2];
sx q[2];
rz(-0.89357251) q[2];
sx q[2];
rz(-0.54235512) q[2];
rz(-0.98313037) q[3];
sx q[3];
rz(-1.977908) q[3];
sx q[3];
rz(-2.2014309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56383175) q[0];
sx q[0];
rz(-1.7153808) q[0];
sx q[0];
rz(-0.78052178) q[0];
rz(1.966194) q[1];
sx q[1];
rz(-0.54882097) q[1];
sx q[1];
rz(1.0343879) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9517453) q[0];
sx q[0];
rz(-1.7380413) q[0];
sx q[0];
rz(0.37519223) q[0];
rz(-pi) q[1];
rz(-1.7138359) q[2];
sx q[2];
rz(-1.6517706) q[2];
sx q[2];
rz(-1.9097999) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3745034) q[1];
sx q[1];
rz(-1.6853635) q[1];
sx q[1];
rz(0.40278816) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54776056) q[3];
sx q[3];
rz(-1.3410853) q[3];
sx q[3];
rz(-1.821777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.15022755) q[2];
sx q[2];
rz(-1.4069724) q[2];
sx q[2];
rz(-3.0312209) q[2];
rz(-1.8433833) q[3];
sx q[3];
rz(-1.7896264) q[3];
sx q[3];
rz(2.8492294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2503535) q[0];
sx q[0];
rz(-2.1284916) q[0];
sx q[0];
rz(-1.2497586) q[0];
rz(-0.091014422) q[1];
sx q[1];
rz(-2.4200771) q[1];
sx q[1];
rz(2.4299842) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2219951) q[0];
sx q[0];
rz(-1.564043) q[0];
sx q[0];
rz(2.3632977) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3468412) q[2];
sx q[2];
rz(-1.0074136) q[2];
sx q[2];
rz(-1.3842441) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6775637) q[1];
sx q[1];
rz(-2.0256126) q[1];
sx q[1];
rz(-1.9536665) q[1];
rz(-2.3859714) q[3];
sx q[3];
rz(-1.8241624) q[3];
sx q[3];
rz(-2.8112335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8083501) q[2];
sx q[2];
rz(-1.4834206) q[2];
sx q[2];
rz(1.9688155) q[2];
rz(-1.5277537) q[3];
sx q[3];
rz(-2.0634191) q[3];
sx q[3];
rz(0.095002256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.461819) q[0];
sx q[0];
rz(-3.0175896) q[0];
sx q[0];
rz(-1.5754196) q[0];
rz(2.7836986) q[1];
sx q[1];
rz(-1.0708258) q[1];
sx q[1];
rz(-2.2391052) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3751463) q[0];
sx q[0];
rz(-1.3129108) q[0];
sx q[0];
rz(-0.93574406) q[0];
rz(-pi) q[1];
rz(2.4975791) q[2];
sx q[2];
rz(-1.5884382) q[2];
sx q[2];
rz(2.0669075) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6186699) q[1];
sx q[1];
rz(-0.82397193) q[1];
sx q[1];
rz(2.0342779) q[1];
x q[2];
rz(-2.1211347) q[3];
sx q[3];
rz(-1.6452879) q[3];
sx q[3];
rz(1.0005152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0365399) q[2];
sx q[2];
rz(-0.5883216) q[2];
sx q[2];
rz(-1.6416637) q[2];
rz(-0.52122742) q[3];
sx q[3];
rz(-0.14385496) q[3];
sx q[3];
rz(2.1874793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70984107) q[0];
sx q[0];
rz(-0.7702282) q[0];
sx q[0];
rz(-1.7256398) q[0];
rz(0.00090986666) q[1];
sx q[1];
rz(-1.4716499) q[1];
sx q[1];
rz(1.6843527) q[1];
rz(1.4191237) q[2];
sx q[2];
rz(-2.5628452) q[2];
sx q[2];
rz(2.645523) q[2];
rz(-0.72003638) q[3];
sx q[3];
rz(-2.2514718) q[3];
sx q[3];
rz(2.9801647) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
