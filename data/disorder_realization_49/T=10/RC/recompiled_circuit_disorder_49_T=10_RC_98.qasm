OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5306659) q[0];
sx q[0];
rz(4.0806169) q[0];
sx q[0];
rz(9.4299849) q[0];
rz(1.7967254) q[1];
sx q[1];
rz(4.2978573) q[1];
sx q[1];
rz(8.2351091) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1253163) q[0];
sx q[0];
rz(-1.2287041) q[0];
sx q[0];
rz(0.54434158) q[0];
x q[1];
rz(-1.7132684) q[2];
sx q[2];
rz(-1.7490897) q[2];
sx q[2];
rz(-1.1536191) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.78566879) q[1];
sx q[1];
rz(-2.7969116) q[1];
sx q[1];
rz(1.1204526) q[1];
rz(-pi) q[2];
rz(-1.0359997) q[3];
sx q[3];
rz(-2.4403604) q[3];
sx q[3];
rz(-2.0555156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4253915) q[2];
sx q[2];
rz(-2.3925245) q[2];
sx q[2];
rz(2.544196) q[2];
rz(-1.3655837) q[3];
sx q[3];
rz(-1.3436915) q[3];
sx q[3];
rz(1.2805773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7146724) q[0];
sx q[0];
rz(-0.5517813) q[0];
sx q[0];
rz(2.8080217) q[0];
rz(-1.0936273) q[1];
sx q[1];
rz(-0.90197864) q[1];
sx q[1];
rz(0.11322583) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2213759) q[0];
sx q[0];
rz(-2.7144103) q[0];
sx q[0];
rz(-1.4772052) q[0];
rz(-pi) q[1];
rz(1.5905154) q[2];
sx q[2];
rz(-1.644051) q[2];
sx q[2];
rz(-0.8578701) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3791703) q[1];
sx q[1];
rz(-1.1602853) q[1];
sx q[1];
rz(1.1344086) q[1];
rz(2.9037335) q[3];
sx q[3];
rz(-0.97994643) q[3];
sx q[3];
rz(-1.295134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1229646) q[2];
sx q[2];
rz(-2.661442) q[2];
sx q[2];
rz(0.22228995) q[2];
rz(0.22953454) q[3];
sx q[3];
rz(-2.4042606) q[3];
sx q[3];
rz(0.078331746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6753321) q[0];
sx q[0];
rz(-1.570913) q[0];
sx q[0];
rz(0.88071841) q[0];
rz(2.0942028) q[1];
sx q[1];
rz(-1.9347582) q[1];
sx q[1];
rz(-3.0139794) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1193082) q[0];
sx q[0];
rz(-0.89953178) q[0];
sx q[0];
rz(-2.9191201) q[0];
x q[1];
rz(2.3348546) q[2];
sx q[2];
rz(-0.85959896) q[2];
sx q[2];
rz(-2.9850609) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.7319273) q[1];
sx q[1];
rz(-1.8095784) q[1];
sx q[1];
rz(-0.084652918) q[1];
rz(1.2283986) q[3];
sx q[3];
rz(-0.6558334) q[3];
sx q[3];
rz(1.1046023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3601274) q[2];
sx q[2];
rz(-2.0121622) q[2];
sx q[2];
rz(-2.8313417) q[2];
rz(-2.7919853) q[3];
sx q[3];
rz(-0.6844371) q[3];
sx q[3];
rz(-1.413697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4629102) q[0];
sx q[0];
rz(-1.8197729) q[0];
sx q[0];
rz(2.983685) q[0];
rz(0.36610106) q[1];
sx q[1];
rz(-0.78916517) q[1];
sx q[1];
rz(2.7691832) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87293738) q[0];
sx q[0];
rz(-0.71632179) q[0];
sx q[0];
rz(1.7563845) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1521058) q[2];
sx q[2];
rz(-1.2876236) q[2];
sx q[2];
rz(0.95209939) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.47632521) q[1];
sx q[1];
rz(-1.421126) q[1];
sx q[1];
rz(1.9515576) q[1];
rz(-0.72752556) q[3];
sx q[3];
rz(-1.1480867) q[3];
sx q[3];
rz(-0.62733516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7238414) q[2];
sx q[2];
rz(-2.4390742) q[2];
sx q[2];
rz(-0.27238971) q[2];
rz(-2.2327936) q[3];
sx q[3];
rz(-0.89509982) q[3];
sx q[3];
rz(-2.6866384) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2588147) q[0];
sx q[0];
rz(-0.60084501) q[0];
sx q[0];
rz(-2.0102665) q[0];
rz(2.5436026) q[1];
sx q[1];
rz(-2.2441041) q[1];
sx q[1];
rz(-2.4647443) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8880496) q[0];
sx q[0];
rz(-1.6497668) q[0];
sx q[0];
rz(0.44102863) q[0];
rz(-1.3555688) q[2];
sx q[2];
rz(-0.85126801) q[2];
sx q[2];
rz(-2.1446251) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.685806) q[1];
sx q[1];
rz(-1.5965441) q[1];
sx q[1];
rz(-0.2548494) q[1];
x q[2];
rz(-1.5464209) q[3];
sx q[3];
rz(-1.9642324) q[3];
sx q[3];
rz(-0.76373053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0040434917) q[2];
sx q[2];
rz(-1.153959) q[2];
sx q[2];
rz(1.1452902) q[2];
rz(-2.9925313) q[3];
sx q[3];
rz(-0.59864932) q[3];
sx q[3];
rz(-2.1242274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0155708) q[0];
sx q[0];
rz(-0.64783043) q[0];
sx q[0];
rz(1.6824678) q[0];
rz(-0.74516621) q[1];
sx q[1];
rz(-0.35062733) q[1];
sx q[1];
rz(0.93313342) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4695078) q[0];
sx q[0];
rz(-0.55969319) q[0];
sx q[0];
rz(-0.99225386) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2007347) q[2];
sx q[2];
rz(-1.0002631) q[2];
sx q[2];
rz(-2.7119315) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7416523) q[1];
sx q[1];
rz(-1.5730904) q[1];
sx q[1];
rz(1.9386577) q[1];
rz(-0.91306367) q[3];
sx q[3];
rz(-1.7115895) q[3];
sx q[3];
rz(1.6933683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3559945) q[2];
sx q[2];
rz(-0.22452536) q[2];
sx q[2];
rz(2.8549426) q[2];
rz(0.24946985) q[3];
sx q[3];
rz(-1.1688066) q[3];
sx q[3];
rz(2.6732895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0934802) q[0];
sx q[0];
rz(-1.1756287) q[0];
sx q[0];
rz(1.4666784) q[0];
rz(1.02007) q[1];
sx q[1];
rz(-1.6721098) q[1];
sx q[1];
rz(2.1405623) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2350378) q[0];
sx q[0];
rz(-2.8585003) q[0];
sx q[0];
rz(-1.7049768) q[0];
rz(-pi) q[1];
x q[1];
rz(2.311003) q[2];
sx q[2];
rz(-2.6579654) q[2];
sx q[2];
rz(2.4110576) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.96879362) q[1];
sx q[1];
rz(-1.7263004) q[1];
sx q[1];
rz(-1.1771727) q[1];
rz(-0.6506728) q[3];
sx q[3];
rz(-1.4015897) q[3];
sx q[3];
rz(0.31864877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.004185685) q[2];
sx q[2];
rz(-1.2265394) q[2];
sx q[2];
rz(-3.0380847) q[2];
rz(2.5497656) q[3];
sx q[3];
rz(-1.3261565) q[3];
sx q[3];
rz(-0.83759585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8822534) q[0];
sx q[0];
rz(-1.0162202) q[0];
sx q[0];
rz(2.5296339) q[0];
rz(1.3423963) q[1];
sx q[1];
rz(-1.1609158) q[1];
sx q[1];
rz(2.7517095) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0355426) q[0];
sx q[0];
rz(-0.70362008) q[0];
sx q[0];
rz(-0.57530595) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50750081) q[2];
sx q[2];
rz(-1.2050873) q[2];
sx q[2];
rz(0.65069228) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4854359) q[1];
sx q[1];
rz(-2.0011144) q[1];
sx q[1];
rz(0.3635316) q[1];
rz(-pi) q[2];
x q[2];
rz(2.249986) q[3];
sx q[3];
rz(-1.7779751) q[3];
sx q[3];
rz(-0.87219119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.3020246) q[2];
sx q[2];
rz(-2.6111157) q[2];
sx q[2];
rz(0.71425444) q[2];
rz(2.2438625) q[3];
sx q[3];
rz(-2.0102746) q[3];
sx q[3];
rz(-0.57202655) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4488688) q[0];
sx q[0];
rz(-0.435193) q[0];
sx q[0];
rz(-2.3642448) q[0];
rz(2.2946987) q[1];
sx q[1];
rz(-2.6234026) q[1];
sx q[1];
rz(-2.8651967) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9696635) q[0];
sx q[0];
rz(-1.164325) q[0];
sx q[0];
rz(-1.3967692) q[0];
rz(0.22819613) q[2];
sx q[2];
rz(-1.7595292) q[2];
sx q[2];
rz(-2.2224094) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8917577) q[1];
sx q[1];
rz(-1.7913622) q[1];
sx q[1];
rz(-2.9794429) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3474465) q[3];
sx q[3];
rz(-0.63784079) q[3];
sx q[3];
rz(2.4362302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.85236621) q[2];
sx q[2];
rz(-1.5091395) q[2];
sx q[2];
rz(-0.12750553) q[2];
rz(-0.036711983) q[3];
sx q[3];
rz(-0.82009411) q[3];
sx q[3];
rz(-1.9071074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7413095) q[0];
sx q[0];
rz(-2.648634) q[0];
sx q[0];
rz(-0.35274831) q[0];
rz(-2.5648975) q[1];
sx q[1];
rz(-0.99027094) q[1];
sx q[1];
rz(-1.030285) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.83537) q[0];
sx q[0];
rz(-1.3754002) q[0];
sx q[0];
rz(-0.72619254) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7416517) q[2];
sx q[2];
rz(-0.72657864) q[2];
sx q[2];
rz(-2.4253997) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.13223091) q[1];
sx q[1];
rz(-1.2538326) q[1];
sx q[1];
rz(-2.6793381) q[1];
rz(-pi) q[2];
rz(-0.61334765) q[3];
sx q[3];
rz(-1.8457335) q[3];
sx q[3];
rz(-3.1230694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0742005) q[2];
sx q[2];
rz(-1.3211162) q[2];
sx q[2];
rz(-2.8004048) q[2];
rz(-0.7406922) q[3];
sx q[3];
rz(-2.1518028) q[3];
sx q[3];
rz(-0.091879524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13070233) q[0];
sx q[0];
rz(-1.1477092) q[0];
sx q[0];
rz(-0.64684091) q[0];
rz(-2.4717992) q[1];
sx q[1];
rz(-0.61129807) q[1];
sx q[1];
rz(1.5652464) q[1];
rz(-0.40217051) q[2];
sx q[2];
rz(-2.675534) q[2];
sx q[2];
rz(0.13327394) q[2];
rz(2.4767247) q[3];
sx q[3];
rz(-0.80857279) q[3];
sx q[3];
rz(1.6894658) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
