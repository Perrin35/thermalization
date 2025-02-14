OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.2494994) q[0];
sx q[0];
rz(-1.9568994) q[0];
sx q[0];
rz(1.1353593) q[0];
rz(-0.42179498) q[1];
sx q[1];
rz(-2.3611958) q[1];
sx q[1];
rz(2.5193522) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.529015) q[0];
sx q[0];
rz(-1.4473995) q[0];
sx q[0];
rz(-1.7181267) q[0];
rz(1.6517992) q[2];
sx q[2];
rz(-2.6930444) q[2];
sx q[2];
rz(2.0392921) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.63401081) q[1];
sx q[1];
rz(-0.9764428) q[1];
sx q[1];
rz(-1.5567012) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1417326) q[3];
sx q[3];
rz(-1.7594764) q[3];
sx q[3];
rz(-2.7267371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3368095) q[2];
sx q[2];
rz(-1.4562891) q[2];
sx q[2];
rz(0.08610227) q[2];
rz(-0.51186776) q[3];
sx q[3];
rz(-2.0972926) q[3];
sx q[3];
rz(-0.69387236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43857273) q[0];
sx q[0];
rz(-0.95412552) q[0];
sx q[0];
rz(2.2254206) q[0];
rz(2.84962) q[1];
sx q[1];
rz(-2.5501854) q[1];
sx q[1];
rz(2.3672262) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60304927) q[0];
sx q[0];
rz(-1.9218245) q[0];
sx q[0];
rz(-2.8013264) q[0];
x q[1];
rz(1.2254481) q[2];
sx q[2];
rz(-1.5651476) q[2];
sx q[2];
rz(-2.408509) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9271597) q[1];
sx q[1];
rz(-1.4478874) q[1];
sx q[1];
rz(-1.6729808) q[1];
x q[2];
rz(-2.1470469) q[3];
sx q[3];
rz(-0.835383) q[3];
sx q[3];
rz(1.6540838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3463717) q[2];
sx q[2];
rz(-2.5796311) q[2];
sx q[2];
rz(2.5566761) q[2];
rz(-1.4784721) q[3];
sx q[3];
rz(-1.1774747) q[3];
sx q[3];
rz(-1.4896721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4855578) q[0];
sx q[0];
rz(-1.0233044) q[0];
sx q[0];
rz(2.5265332) q[0];
rz(-2.8133605) q[1];
sx q[1];
rz(-1.0098207) q[1];
sx q[1];
rz(-1.0438017) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39343483) q[0];
sx q[0];
rz(-1.5183477) q[0];
sx q[0];
rz(1.5758874) q[0];
rz(-pi) q[1];
rz(-2.7132665) q[2];
sx q[2];
rz(-1.5751683) q[2];
sx q[2];
rz(-2.0010889) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.97745132) q[1];
sx q[1];
rz(-0.79557997) q[1];
sx q[1];
rz(-2.1358016) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39643328) q[3];
sx q[3];
rz(-1.9045489) q[3];
sx q[3];
rz(2.9891356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7528167) q[2];
sx q[2];
rz(-2.7478168) q[2];
sx q[2];
rz(2.8326995) q[2];
rz(0.4872407) q[3];
sx q[3];
rz(-1.7083218) q[3];
sx q[3];
rz(-2.2545831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.3798645) q[0];
sx q[0];
rz(-2.3966615) q[0];
sx q[0];
rz(0.47759011) q[0];
rz(1.8567122) q[1];
sx q[1];
rz(-1.4297337) q[1];
sx q[1];
rz(-2.6755948) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14529252) q[0];
sx q[0];
rz(-1.0178653) q[0];
sx q[0];
rz(-1.6607222) q[0];
x q[1];
rz(2.6902507) q[2];
sx q[2];
rz(-1.1548496) q[2];
sx q[2];
rz(-2.0030264) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1124777) q[1];
sx q[1];
rz(-2.3211109) q[1];
sx q[1];
rz(-0.047767834) q[1];
rz(-pi) q[2];
rz(0.68434207) q[3];
sx q[3];
rz(-2.0668233) q[3];
sx q[3];
rz(3.0471714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8183225) q[2];
sx q[2];
rz(-0.44760901) q[2];
sx q[2];
rz(-1.6935879) q[2];
rz(-1.7285534) q[3];
sx q[3];
rz(-1.5331242) q[3];
sx q[3];
rz(-1.1477227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2959761) q[0];
sx q[0];
rz(-2.4281261) q[0];
sx q[0];
rz(3.0699068) q[0];
rz(1.4777615) q[1];
sx q[1];
rz(-1.7774589) q[1];
sx q[1];
rz(-0.62612265) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5324931) q[0];
sx q[0];
rz(-0.83864826) q[0];
sx q[0];
rz(-0.12153836) q[0];
x q[1];
rz(-1.6501872) q[2];
sx q[2];
rz(-0.65020665) q[2];
sx q[2];
rz(-2.2967867) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71264919) q[1];
sx q[1];
rz(-0.82368851) q[1];
sx q[1];
rz(2.9045543) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6780217) q[3];
sx q[3];
rz(-2.7509973) q[3];
sx q[3];
rz(2.5724346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6757297) q[2];
sx q[2];
rz(-2.6079123) q[2];
sx q[2];
rz(-1.6443058) q[2];
rz(2.7346482) q[3];
sx q[3];
rz(-2.1502083) q[3];
sx q[3];
rz(-0.86010325) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1740455) q[0];
sx q[0];
rz(-1.8438735) q[0];
sx q[0];
rz(-0.26200727) q[0];
rz(1.920248) q[1];
sx q[1];
rz(-1.5120993) q[1];
sx q[1];
rz(1.2006203) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1977494) q[0];
sx q[0];
rz(-2.245994) q[0];
sx q[0];
rz(-2.4989216) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1706743) q[2];
sx q[2];
rz(-0.80764233) q[2];
sx q[2];
rz(1.7852448) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8194325) q[1];
sx q[1];
rz(-1.3272078) q[1];
sx q[1];
rz(-2.0889682) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2030199) q[3];
sx q[3];
rz(-0.098417131) q[3];
sx q[3];
rz(-0.84151387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0210586) q[2];
sx q[2];
rz(-2.6594682) q[2];
sx q[2];
rz(0.47312197) q[2];
rz(-2.1038697) q[3];
sx q[3];
rz(-2.7389052) q[3];
sx q[3];
rz(1.13824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6556019) q[0];
sx q[0];
rz(-0.36360535) q[0];
sx q[0];
rz(-2.3714491) q[0];
rz(-0.4153525) q[1];
sx q[1];
rz(-1.628592) q[1];
sx q[1];
rz(1.279668) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3949385) q[0];
sx q[0];
rz(-1.5956559) q[0];
sx q[0];
rz(-0.8381054) q[0];
rz(-2.0023022) q[2];
sx q[2];
rz(-1.6864634) q[2];
sx q[2];
rz(-2.3097599) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.34977594) q[1];
sx q[1];
rz(-2.8703051) q[1];
sx q[1];
rz(-0.94601814) q[1];
rz(-pi) q[2];
rz(0.47904738) q[3];
sx q[3];
rz(-1.5622514) q[3];
sx q[3];
rz(2.6984877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8659849) q[2];
sx q[2];
rz(-1.393968) q[2];
sx q[2];
rz(2.7371791) q[2];
rz(1.154493) q[3];
sx q[3];
rz(-1.5279852) q[3];
sx q[3];
rz(1.9563458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.6605717) q[0];
sx q[0];
rz(-2.4206929) q[0];
sx q[0];
rz(-0.50203669) q[0];
rz(1.0188811) q[1];
sx q[1];
rz(-1.6682245) q[1];
sx q[1];
rz(-0.92207164) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2699091) q[0];
sx q[0];
rz(-0.26342621) q[0];
sx q[0];
rz(-2.7676537) q[0];
x q[1];
rz(-1.8174632) q[2];
sx q[2];
rz(-0.96117678) q[2];
sx q[2];
rz(-2.7423046) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0623297) q[1];
sx q[1];
rz(-0.070186071) q[1];
sx q[1];
rz(0.90026469) q[1];
rz(-pi) q[2];
rz(-1.5516547) q[3];
sx q[3];
rz(-0.12775207) q[3];
sx q[3];
rz(0.16523339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8036538) q[2];
sx q[2];
rz(-1.8622082) q[2];
sx q[2];
rz(0.12019084) q[2];
rz(0.0029314824) q[3];
sx q[3];
rz(-0.35792297) q[3];
sx q[3];
rz(-3.1357989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040076479) q[0];
sx q[0];
rz(-2.4637971) q[0];
sx q[0];
rz(-1.5090322) q[0];
rz(-2.312233) q[1];
sx q[1];
rz(-0.50193915) q[1];
sx q[1];
rz(-2.0109743) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7735159) q[0];
sx q[0];
rz(-0.46958554) q[0];
sx q[0];
rz(2.3115524) q[0];
rz(-pi) q[1];
rz(-2.4254029) q[2];
sx q[2];
rz(-1.9702868) q[2];
sx q[2];
rz(-1.3376118) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77525951) q[1];
sx q[1];
rz(-2.383552) q[1];
sx q[1];
rz(-2.232928) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9958903) q[3];
sx q[3];
rz(-0.55820528) q[3];
sx q[3];
rz(-0.74287117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.067554615) q[2];
sx q[2];
rz(-2.270547) q[2];
sx q[2];
rz(1.3646431) q[2];
rz(0.015297628) q[3];
sx q[3];
rz(-2.3099895) q[3];
sx q[3];
rz(-1.3219705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(2.7010654) q[0];
sx q[0];
rz(-2.4897713) q[0];
sx q[0];
rz(-0.1524674) q[0];
rz(-1.7561779) q[1];
sx q[1];
rz(-2.0962174) q[1];
sx q[1];
rz(-2.8130031) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91157527) q[0];
sx q[0];
rz(-1.5571638) q[0];
sx q[0];
rz(1.5799119) q[0];
rz(-pi) q[1];
rz(-2.3797437) q[2];
sx q[2];
rz(-1.3947316) q[2];
sx q[2];
rz(-3.0886285) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.32544245) q[1];
sx q[1];
rz(-2.901536) q[1];
sx q[1];
rz(-1.6961967) q[1];
x q[2];
rz(-1.2166546) q[3];
sx q[3];
rz(-2.1855857) q[3];
sx q[3];
rz(1.5686017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.80455989) q[2];
sx q[2];
rz(-2.0064662) q[2];
sx q[2];
rz(-2.0648773) q[2];
rz(-0.95450258) q[3];
sx q[3];
rz(-1.4408305) q[3];
sx q[3];
rz(-0.29630989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3267219) q[0];
sx q[0];
rz(-2.1537415) q[0];
sx q[0];
rz(1.4065773) q[0];
rz(-1.3181435) q[1];
sx q[1];
rz(-1.6271918) q[1];
sx q[1];
rz(0.79961332) q[1];
rz(-1.5158657) q[2];
sx q[2];
rz(-1.5274897) q[2];
sx q[2];
rz(-2.4804583) q[2];
rz(-1.7794505) q[3];
sx q[3];
rz(-0.91051523) q[3];
sx q[3];
rz(-0.12403535) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
