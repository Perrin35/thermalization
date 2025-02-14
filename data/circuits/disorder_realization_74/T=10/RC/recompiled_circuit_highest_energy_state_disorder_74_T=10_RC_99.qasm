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
rz(-2.9086034) q[0];
rz(-1.5268582) q[1];
sx q[1];
rz(-2.0695217) q[1];
sx q[1];
rz(1.1195247) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6699864) q[0];
sx q[0];
rz(-1.3363839) q[0];
sx q[0];
rz(3.1283911) q[0];
x q[1];
rz(1.5824116) q[2];
sx q[2];
rz(-2.1321974) q[2];
sx q[2];
rz(2.7593768) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0927375) q[1];
sx q[1];
rz(-1.6462781) q[1];
sx q[1];
rz(1.2176179) q[1];
rz(-pi) q[2];
rz(0.60910881) q[3];
sx q[3];
rz(-1.850633) q[3];
sx q[3];
rz(-0.5718872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.69922525) q[2];
sx q[2];
rz(-2.0825601) q[2];
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
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3498822) q[0];
sx q[0];
rz(-0.078502027) q[0];
sx q[0];
rz(-3.0872524) q[0];
rz(-2.9229274) q[1];
sx q[1];
rz(-1.668914) q[1];
sx q[1];
rz(2.7770619) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35848356) q[0];
sx q[0];
rz(-2.3158584) q[0];
sx q[0];
rz(0.70348212) q[0];
rz(-pi) q[1];
rz(-2.1991992) q[2];
sx q[2];
rz(-1.9004873) q[2];
sx q[2];
rz(-0.18131944) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.042885429) q[1];
sx q[1];
rz(-1.4946997) q[1];
sx q[1];
rz(2.5185381) q[1];
rz(-0.18192911) q[3];
sx q[3];
rz(-1.9046138) q[3];
sx q[3];
rz(0.16063375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.24078044) q[2];
sx q[2];
rz(-1.6996926) q[2];
sx q[2];
rz(-2.0330632) q[2];
rz(2.945914) q[3];
sx q[3];
rz(-1.207573) q[3];
sx q[3];
rz(1.4669363) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7473258) q[0];
sx q[0];
rz(-1.0402004) q[0];
sx q[0];
rz(-1.0362097) q[0];
rz(0.58468435) q[1];
sx q[1];
rz(-1.6744637) q[1];
sx q[1];
rz(-0.55589693) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7202458) q[0];
sx q[0];
rz(-1.1555536) q[0];
sx q[0];
rz(-0.74851997) q[0];
x q[1];
rz(-2.0228407) q[2];
sx q[2];
rz(-0.085918203) q[2];
sx q[2];
rz(0.8762067) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4208918) q[1];
sx q[1];
rz(-1.4399505) q[1];
sx q[1];
rz(2.2761005) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44938748) q[3];
sx q[3];
rz(-1.7996178) q[3];
sx q[3];
rz(-0.23923161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3802152) q[2];
sx q[2];
rz(-0.20647241) q[2];
sx q[2];
rz(-1.9492487) q[2];
rz(0.10382593) q[3];
sx q[3];
rz(-1.1622279) q[3];
sx q[3];
rz(-1.9942358) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.035606774) q[0];
sx q[0];
rz(-2.5150531) q[0];
sx q[0];
rz(-1.71738) q[0];
rz(-0.72319889) q[1];
sx q[1];
rz(-2.9056748) q[1];
sx q[1];
rz(-0.22044388) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0766692) q[0];
sx q[0];
rz(-0.60416302) q[0];
sx q[0];
rz(-2.2518065) q[0];
x q[1];
rz(1.1197243) q[2];
sx q[2];
rz(-0.60904087) q[2];
sx q[2];
rz(-2.232617) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4749149) q[1];
sx q[1];
rz(-1.4581469) q[1];
sx q[1];
rz(-1.2126318) q[1];
rz(-pi) q[2];
rz(-1.9842159) q[3];
sx q[3];
rz(-1.9154906) q[3];
sx q[3];
rz(-0.34806309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7580938) q[2];
sx q[2];
rz(-1.3145964) q[2];
sx q[2];
rz(0.51551762) q[2];
rz(-0.085144194) q[3];
sx q[3];
rz(-0.65535039) q[3];
sx q[3];
rz(-0.83318025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67007095) q[0];
sx q[0];
rz(-2.9586198) q[0];
sx q[0];
rz(0.66437379) q[0];
rz(-2.6145256) q[1];
sx q[1];
rz(-1.138569) q[1];
sx q[1];
rz(1.9648431) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.697061) q[0];
sx q[0];
rz(-0.089655487) q[0];
sx q[0];
rz(-1.8898925) q[0];
x q[1];
rz(2.2559153) q[2];
sx q[2];
rz(-1.204245) q[2];
sx q[2];
rz(1.7358071) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.47065319) q[1];
sx q[1];
rz(-1.8090418) q[1];
sx q[1];
rz(-1.5508503) q[1];
rz(-1.2958636) q[3];
sx q[3];
rz(-2.5524013) q[3];
sx q[3];
rz(-1.9269892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1145757) q[2];
sx q[2];
rz(-0.21275529) q[2];
sx q[2];
rz(2.2034755) q[2];
rz(2.0959057) q[3];
sx q[3];
rz(-0.18200471) q[3];
sx q[3];
rz(0.81030455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3351347) q[0];
sx q[0];
rz(-1.198575) q[0];
sx q[0];
rz(-1.2499811) q[0];
rz(-0.20420034) q[1];
sx q[1];
rz(-0.67664346) q[1];
sx q[1];
rz(-0.85320371) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6224504) q[0];
sx q[0];
rz(-1.8565489) q[0];
sx q[0];
rz(-2.7842872) q[0];
rz(-pi) q[1];
rz(-2.5218042) q[2];
sx q[2];
rz(-1.0519805) q[2];
sx q[2];
rz(-0.74558115) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1386765) q[1];
sx q[1];
rz(-0.87453784) q[1];
sx q[1];
rz(2.584444) q[1];
rz(-0.12746396) q[3];
sx q[3];
rz(-1.9098305) q[3];
sx q[3];
rz(1.4900498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.044067232) q[2];
sx q[2];
rz(-1.3419469) q[2];
sx q[2];
rz(-1.9419144) q[2];
rz(-2.1357644) q[3];
sx q[3];
rz(-2.2483716) q[3];
sx q[3];
rz(2.306166) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6658039) q[0];
sx q[0];
rz(-2.8732193) q[0];
sx q[0];
rz(0.98261181) q[0];
rz(-1.8334552) q[1];
sx q[1];
rz(-1.9255226) q[1];
sx q[1];
rz(1.4391724) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2932201) q[0];
sx q[0];
rz(-1.9008027) q[0];
sx q[0];
rz(2.7820935) q[0];
rz(-pi) q[1];
rz(1.4365044) q[2];
sx q[2];
rz(-0.52596131) q[2];
sx q[2];
rz(0.44912042) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8285259) q[1];
sx q[1];
rz(-1.2432352) q[1];
sx q[1];
rz(2.2799973) q[1];
rz(-0.26994522) q[3];
sx q[3];
rz(-2.649123) q[3];
sx q[3];
rz(0.47583285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.16284379) q[2];
sx q[2];
rz(-0.89357251) q[2];
sx q[2];
rz(-0.54235512) q[2];
rz(-2.1584623) q[3];
sx q[3];
rz(-1.1636846) q[3];
sx q[3];
rz(0.94016176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
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
rz(0.56383175) q[0];
sx q[0];
rz(-1.7153808) q[0];
sx q[0];
rz(-0.78052178) q[0];
rz(-1.1753987) q[1];
sx q[1];
rz(-2.5927717) q[1];
sx q[1];
rz(-1.0343879) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7808409) q[0];
sx q[0];
rz(-2.7324317) q[0];
sx q[0];
rz(-0.43171127) q[0];
rz(3.0597866) q[2];
sx q[2];
rz(-1.4282287) q[2];
sx q[2];
rz(-0.35065251) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.24496291) q[1];
sx q[1];
rz(-1.9707929) q[1];
sx q[1];
rz(-1.4463615) q[1];
x q[2];
rz(0.42203776) q[3];
sx q[3];
rz(-2.5521818) q[3];
sx q[3];
rz(-0.60822076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9913651) q[2];
sx q[2];
rz(-1.7346202) q[2];
sx q[2];
rz(3.0312209) q[2];
rz(-1.8433833) q[3];
sx q[3];
rz(-1.7896264) q[3];
sx q[3];
rz(2.8492294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2503535) q[0];
sx q[0];
rz(-1.013101) q[0];
sx q[0];
rz(1.891834) q[0];
rz(-0.091014422) q[1];
sx q[1];
rz(-2.4200771) q[1];
sx q[1];
rz(-0.7116085) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34214333) q[0];
sx q[0];
rz(-2.3490688) q[0];
sx q[0];
rz(1.5613129) q[0];
rz(-0.79475148) q[2];
sx q[2];
rz(-2.134179) q[2];
sx q[2];
rz(1.3842441) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2818799) q[1];
sx q[1];
rz(-1.2285474) q[1];
sx q[1];
rz(-0.48515353) q[1];
x q[2];
rz(-2.3859714) q[3];
sx q[3];
rz(-1.3174302) q[3];
sx q[3];
rz(2.8112335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8083501) q[2];
sx q[2];
rz(-1.658172) q[2];
sx q[2];
rz(1.9688155) q[2];
rz(-1.6138389) q[3];
sx q[3];
rz(-2.0634191) q[3];
sx q[3];
rz(3.0465904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.461819) q[0];
sx q[0];
rz(-0.12400308) q[0];
sx q[0];
rz(-1.5754196) q[0];
rz(-2.7836986) q[1];
sx q[1];
rz(-2.0707668) q[1];
sx q[1];
rz(0.90248743) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47115147) q[0];
sx q[0];
rz(-2.4629399) q[0];
sx q[0];
rz(-1.9891692) q[0];
x q[1];
rz(-3.112215) q[2];
sx q[2];
rz(-2.4973719) q[2];
sx q[2];
rz(-0.51960301) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6186699) q[1];
sx q[1];
rz(-0.82397193) q[1];
sx q[1];
rz(2.0342779) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4290513) q[3];
sx q[3];
rz(-0.55484164) q[3];
sx q[3];
rz(2.6920163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0365399) q[2];
sx q[2];
rz(-0.5883216) q[2];
sx q[2];
rz(1.6416637) q[2];
rz(-0.52122742) q[3];
sx q[3];
rz(-0.14385496) q[3];
sx q[3];
rz(-0.95411333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70984107) q[0];
sx q[0];
rz(-2.3713645) q[0];
sx q[0];
rz(1.4159528) q[0];
rz(0.00090986666) q[1];
sx q[1];
rz(-1.4716499) q[1];
sx q[1];
rz(1.6843527) q[1];
rz(-1.4191237) q[2];
sx q[2];
rz(-0.57874741) q[2];
sx q[2];
rz(-0.49606965) q[2];
rz(-2.2542027) q[3];
sx q[3];
rz(-0.94684623) q[3];
sx q[3];
rz(2.0317247) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
