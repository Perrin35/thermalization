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
rz(-1.2008774) q[0];
sx q[0];
rz(-0.90056363) q[0];
sx q[0];
rz(-0.23298921) q[0];
rz(1.6147344) q[1];
sx q[1];
rz(-1.072071) q[1];
sx q[1];
rz(-1.1195247) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0454691) q[0];
sx q[0];
rz(-1.5579559) q[0];
sx q[0];
rz(1.8052285) q[0];
rz(-0.018466516) q[2];
sx q[2];
rz(-2.5800843) q[2];
sx q[2];
rz(-0.36040053) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.865316) q[1];
sx q[1];
rz(-0.36082339) q[1];
sx q[1];
rz(-1.3555384) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60910881) q[3];
sx q[3];
rz(-1.2909596) q[3];
sx q[3];
rz(0.5718872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.69922525) q[2];
sx q[2];
rz(-1.0590326) q[2];
sx q[2];
rz(-0.11288682) q[2];
rz(-2.8804307) q[3];
sx q[3];
rz(-1.7966725) q[3];
sx q[3];
rz(1.8908709) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3498822) q[0];
sx q[0];
rz(-3.0630906) q[0];
sx q[0];
rz(-3.0872524) q[0];
rz(-2.9229274) q[1];
sx q[1];
rz(-1.668914) q[1];
sx q[1];
rz(-0.36453077) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.88663) q[0];
sx q[0];
rz(-0.97575649) q[0];
sx q[0];
rz(2.182385) q[0];
x q[1];
rz(-2.0979375) q[2];
sx q[2];
rz(-2.4424565) q[2];
sx q[2];
rz(1.8086036) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0987072) q[1];
sx q[1];
rz(-1.646893) q[1];
sx q[1];
rz(2.5185381) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2317936) q[3];
sx q[3];
rz(-1.7425797) q[3];
sx q[3];
rz(1.671227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9008122) q[2];
sx q[2];
rz(-1.6996926) q[2];
sx q[2];
rz(2.0330632) q[2];
rz(-0.19567868) q[3];
sx q[3];
rz(-1.207573) q[3];
sx q[3];
rz(-1.6746563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7473258) q[0];
sx q[0];
rz(-2.1013923) q[0];
sx q[0];
rz(2.105383) q[0];
rz(0.58468435) q[1];
sx q[1];
rz(-1.4671289) q[1];
sx q[1];
rz(-2.5856957) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5824048) q[0];
sx q[0];
rz(-2.3055861) q[0];
sx q[0];
rz(-0.57484267) q[0];
rz(1.6481208) q[2];
sx q[2];
rz(-1.5333042) q[2];
sx q[2];
rz(-2.897597) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7207009) q[1];
sx q[1];
rz(-1.7016422) q[1];
sx q[1];
rz(2.2761005) q[1];
rz(1.8238277) q[3];
sx q[3];
rz(-2.0076499) q[3];
sx q[3];
rz(1.7010613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3802152) q[2];
sx q[2];
rz(-0.20647241) q[2];
sx q[2];
rz(-1.192344) q[2];
rz(-0.10382593) q[3];
sx q[3];
rz(-1.1622279) q[3];
sx q[3];
rz(-1.1473568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.035606774) q[0];
sx q[0];
rz(-0.62653956) q[0];
sx q[0];
rz(-1.4242127) q[0];
rz(-0.72319889) q[1];
sx q[1];
rz(-2.9056748) q[1];
sx q[1];
rz(2.9211488) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0940277) q[0];
sx q[0];
rz(-1.936543) q[0];
sx q[0];
rz(-2.0630552) q[0];
rz(-1.0102369) q[2];
sx q[2];
rz(-1.3187485) q[2];
sx q[2];
rz(0.28365669) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4749149) q[1];
sx q[1];
rz(-1.4581469) q[1];
sx q[1];
rz(-1.2126318) q[1];
rz(-pi) q[2];
rz(-1.1573767) q[3];
sx q[3];
rz(-1.9154906) q[3];
sx q[3];
rz(-2.7935296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7580938) q[2];
sx q[2];
rz(-1.8269962) q[2];
sx q[2];
rz(-0.51551762) q[2];
rz(0.085144194) q[3];
sx q[3];
rz(-0.65535039) q[3];
sx q[3];
rz(-2.3084124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67007095) q[0];
sx q[0];
rz(-0.18297289) q[0];
sx q[0];
rz(2.4772189) q[0];
rz(-0.5270671) q[1];
sx q[1];
rz(-1.138569) q[1];
sx q[1];
rz(-1.9648431) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4445317) q[0];
sx q[0];
rz(-0.089655487) q[0];
sx q[0];
rz(1.2517002) q[0];
rz(-pi) q[1];
rz(2.1161304) q[2];
sx q[2];
rz(-2.3787913) q[2];
sx q[2];
rz(-0.57833407) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.47065319) q[1];
sx q[1];
rz(-1.8090418) q[1];
sx q[1];
rz(1.5508503) q[1];
rz(-pi) q[2];
rz(-0.17950157) q[3];
sx q[3];
rz(-2.13509) q[3];
sx q[3];
rz(1.5416984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1145757) q[2];
sx q[2];
rz(-0.21275529) q[2];
sx q[2];
rz(-2.2034755) q[2];
rz(2.0959057) q[3];
sx q[3];
rz(-2.9595879) q[3];
sx q[3];
rz(2.3312881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.806458) q[0];
sx q[0];
rz(-1.198575) q[0];
sx q[0];
rz(-1.8916116) q[0];
rz(-2.9373923) q[1];
sx q[1];
rz(-2.4649492) q[1];
sx q[1];
rz(-0.85320371) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9850904) q[0];
sx q[0];
rz(-1.2285875) q[0];
sx q[0];
rz(1.2669105) q[0];
x q[1];
rz(2.1825024) q[2];
sx q[2];
rz(-2.0995127) q[2];
sx q[2];
rz(-2.6564646) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0029162) q[1];
sx q[1];
rz(-0.87453784) q[1];
sx q[1];
rz(0.55714861) q[1];
rz(1.2291992) q[3];
sx q[3];
rz(-1.4506243) q[3];
sx q[3];
rz(-0.12334331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0975254) q[2];
sx q[2];
rz(-1.3419469) q[2];
sx q[2];
rz(1.1996783) q[2];
rz(1.0058282) q[3];
sx q[3];
rz(-0.89322105) q[3];
sx q[3];
rz(0.83542663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6658039) q[0];
sx q[0];
rz(-2.8732193) q[0];
sx q[0];
rz(2.1589808) q[0];
rz(1.8334552) q[1];
sx q[1];
rz(-1.2160701) q[1];
sx q[1];
rz(-1.7024202) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4339695) q[0];
sx q[0];
rz(-2.6585007) q[0];
sx q[0];
rz(0.77204319) q[0];
rz(1.0487532) q[2];
sx q[2];
rz(-1.5035275) q[2];
sx q[2];
rz(-1.0053588) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3130668) q[1];
sx q[1];
rz(-1.2432352) q[1];
sx q[1];
rz(2.2799973) q[1];
rz(2.8716474) q[3];
sx q[3];
rz(-2.649123) q[3];
sx q[3];
rz(-2.6657598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.16284379) q[2];
sx q[2];
rz(-0.89357251) q[2];
sx q[2];
rz(-0.54235512) q[2];
rz(0.98313037) q[3];
sx q[3];
rz(-1.1636846) q[3];
sx q[3];
rz(0.94016176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56383175) q[0];
sx q[0];
rz(-1.7153808) q[0];
sx q[0];
rz(-0.78052178) q[0];
rz(-1.966194) q[1];
sx q[1];
rz(-0.54882097) q[1];
sx q[1];
rz(2.1072047) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3607517) q[0];
sx q[0];
rz(-0.409161) q[0];
sx q[0];
rz(0.43171127) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.081806094) q[2];
sx q[2];
rz(-1.4282287) q[2];
sx q[2];
rz(-0.35065251) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.065830439) q[1];
sx q[1];
rz(-0.41790634) q[1];
sx q[1];
rz(2.8560545) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3034399) q[3];
sx q[3];
rz(-2.1026097) q[3];
sx q[3];
rz(0.11296266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15022755) q[2];
sx q[2];
rz(-1.7346202) q[2];
sx q[2];
rz(0.11037174) q[2];
rz(-1.2982093) q[3];
sx q[3];
rz(-1.7896264) q[3];
sx q[3];
rz(-2.8492294) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2503535) q[0];
sx q[0];
rz(-1.013101) q[0];
sx q[0];
rz(-1.2497586) q[0];
rz(3.0505782) q[1];
sx q[1];
rz(-2.4200771) q[1];
sx q[1];
rz(-0.7116085) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7994493) q[0];
sx q[0];
rz(-2.3490688) q[0];
sx q[0];
rz(-1.5802797) q[0];
rz(-pi) q[1];
rz(2.3468412) q[2];
sx q[2];
rz(-1.0074136) q[2];
sx q[2];
rz(1.7573485) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8636057) q[1];
sx q[1];
rz(-2.5558439) q[1];
sx q[1];
rz(-0.65237712) q[1];
rz(-pi) q[2];
rz(-0.75562128) q[3];
sx q[3];
rz(-1.3174302) q[3];
sx q[3];
rz(-2.8112335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8083501) q[2];
sx q[2];
rz(-1.4834206) q[2];
sx q[2];
rz(-1.1727772) q[2];
rz(1.6138389) q[3];
sx q[3];
rz(-2.0634191) q[3];
sx q[3];
rz(0.095002256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67977366) q[0];
sx q[0];
rz(-0.12400308) q[0];
sx q[0];
rz(1.5754196) q[0];
rz(-2.7836986) q[1];
sx q[1];
rz(-1.0708258) q[1];
sx q[1];
rz(2.2391052) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6704412) q[0];
sx q[0];
rz(-2.4629399) q[0];
sx q[0];
rz(1.1524234) q[0];
rz(-pi) q[1];
rz(0.029377653) q[2];
sx q[2];
rz(-0.64422078) q[2];
sx q[2];
rz(-2.6219896) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6186699) q[1];
sx q[1];
rz(-0.82397193) q[1];
sx q[1];
rz(1.1073148) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4290513) q[3];
sx q[3];
rz(-0.55484164) q[3];
sx q[3];
rz(-2.6920163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0365399) q[2];
sx q[2];
rz(-0.5883216) q[2];
sx q[2];
rz(1.6416637) q[2];
rz(0.52122742) q[3];
sx q[3];
rz(-0.14385496) q[3];
sx q[3];
rz(0.95411333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.4317516) q[0];
sx q[0];
rz(-2.3713645) q[0];
sx q[0];
rz(1.4159528) q[0];
rz(3.1406828) q[1];
sx q[1];
rz(-1.6699427) q[1];
sx q[1];
rz(-1.4572399) q[1];
rz(-2.1442689) q[2];
sx q[2];
rz(-1.6535342) q[2];
sx q[2];
rz(0.94746305) q[2];
rz(2.2542027) q[3];
sx q[3];
rz(-2.1947464) q[3];
sx q[3];
rz(-1.1098679) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
