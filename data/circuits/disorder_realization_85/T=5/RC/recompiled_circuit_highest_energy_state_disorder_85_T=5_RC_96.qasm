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
rz(2.3020701) q[0];
sx q[0];
rz(-1.5538202) q[0];
sx q[0];
rz(2.177218) q[0];
rz(0.40274611) q[1];
sx q[1];
rz(-1.6778814) q[1];
sx q[1];
rz(-0.95199624) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6168686) q[0];
sx q[0];
rz(-1.6646806) q[0];
sx q[0];
rz(2.062172) q[0];
rz(-pi) q[1];
rz(-0.098564221) q[2];
sx q[2];
rz(-1.1505039) q[2];
sx q[2];
rz(1.9914876) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1237028) q[1];
sx q[1];
rz(-0.54819878) q[1];
sx q[1];
rz(-2.3831559) q[1];
x q[2];
rz(-2.3856569) q[3];
sx q[3];
rz(-1.9256251) q[3];
sx q[3];
rz(-0.95099041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6160148) q[2];
sx q[2];
rz(-1.792045) q[2];
sx q[2];
rz(-0.40526059) q[2];
rz(1.1303834) q[3];
sx q[3];
rz(-0.78961343) q[3];
sx q[3];
rz(-1.9270886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90527117) q[0];
sx q[0];
rz(-3.1089678) q[0];
sx q[0];
rz(0.76500934) q[0];
rz(-2.8969823) q[1];
sx q[1];
rz(-1.2449934) q[1];
sx q[1];
rz(0.6388706) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9332704) q[0];
sx q[0];
rz(-0.92678419) q[0];
sx q[0];
rz(-1.6025376) q[0];
x q[1];
rz(-2.6354796) q[2];
sx q[2];
rz(-2.4547336) q[2];
sx q[2];
rz(-3.0619301) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8117026) q[1];
sx q[1];
rz(-2.5479758) q[1];
sx q[1];
rz(2.4827186) q[1];
rz(-pi) q[2];
rz(1.2876533) q[3];
sx q[3];
rz(-2.3842709) q[3];
sx q[3];
rz(0.18268302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2834393) q[2];
sx q[2];
rz(-2.6946113) q[2];
sx q[2];
rz(-2.8972674) q[2];
rz(-2.0325932) q[3];
sx q[3];
rz(-1.2764443) q[3];
sx q[3];
rz(0.23787704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.37794381) q[0];
sx q[0];
rz(-2.0574594) q[0];
sx q[0];
rz(-1.9042683) q[0];
rz(1.7297277) q[1];
sx q[1];
rz(-2.6928163) q[1];
sx q[1];
rz(-1.3317187) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.004382172) q[0];
sx q[0];
rz(-1.4829652) q[0];
sx q[0];
rz(-2.3669404) q[0];
rz(-pi) q[1];
rz(1.2197415) q[2];
sx q[2];
rz(-1.3655541) q[2];
sx q[2];
rz(-0.83265162) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0983372) q[1];
sx q[1];
rz(-1.5374633) q[1];
sx q[1];
rz(-0.055524932) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94966763) q[3];
sx q[3];
rz(-1.0499718) q[3];
sx q[3];
rz(-2.4609712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9776913) q[2];
sx q[2];
rz(-2.0745664) q[2];
sx q[2];
rz(2.4816371) q[2];
rz(-1.6948949) q[3];
sx q[3];
rz(-1.7132672) q[3];
sx q[3];
rz(2.0383294) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25809836) q[0];
sx q[0];
rz(-0.75089184) q[0];
sx q[0];
rz(1.1335491) q[0];
rz(-2.6253888) q[1];
sx q[1];
rz(-1.8323003) q[1];
sx q[1];
rz(2.5925327) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7507197) q[0];
sx q[0];
rz(-1.735512) q[0];
sx q[0];
rz(0.14589845) q[0];
rz(-pi) q[1];
x q[1];
rz(1.334515) q[2];
sx q[2];
rz(-0.60271128) q[2];
sx q[2];
rz(0.38379764) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8827277) q[1];
sx q[1];
rz(-2.8587864) q[1];
sx q[1];
rz(-2.905719) q[1];
rz(-pi) q[2];
rz(-0.28259773) q[3];
sx q[3];
rz(-1.188084) q[3];
sx q[3];
rz(-3.0846918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0893726) q[2];
sx q[2];
rz(-1.0129656) q[2];
sx q[2];
rz(0.29903665) q[2];
rz(2.5236169) q[3];
sx q[3];
rz(-2.0682014) q[3];
sx q[3];
rz(-1.3581902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1842136) q[0];
sx q[0];
rz(-3.1035677) q[0];
sx q[0];
rz(-2.6755565) q[0];
rz(0.86712962) q[1];
sx q[1];
rz(-2.6301818) q[1];
sx q[1];
rz(2.3066511) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3059495) q[0];
sx q[0];
rz(-1.5993274) q[0];
sx q[0];
rz(-3.0898407) q[0];
rz(-1.4915799) q[2];
sx q[2];
rz(-1.4822373) q[2];
sx q[2];
rz(1.4350325) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1620491) q[1];
sx q[1];
rz(-2.0503373) q[1];
sx q[1];
rz(-0.060258106) q[1];
rz(-pi) q[2];
rz(-0.72839177) q[3];
sx q[3];
rz(-2.0660704) q[3];
sx q[3];
rz(-0.033897696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.20778188) q[2];
sx q[2];
rz(-2.2489579) q[2];
sx q[2];
rz(3.1100698) q[2];
rz(2.3837714) q[3];
sx q[3];
rz(-2.0407245) q[3];
sx q[3];
rz(2.5789564) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9071478) q[0];
sx q[0];
rz(-1.7666768) q[0];
sx q[0];
rz(0.031540792) q[0];
rz(-0.08983287) q[1];
sx q[1];
rz(-2.641771) q[1];
sx q[1];
rz(2.8057742) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9851097) q[0];
sx q[0];
rz(-1.8976189) q[0];
sx q[0];
rz(2.7669897) q[0];
rz(-pi) q[1];
rz(-1.01891) q[2];
sx q[2];
rz(-1.2794211) q[2];
sx q[2];
rz(-0.62770972) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5988859) q[1];
sx q[1];
rz(-0.60208354) q[1];
sx q[1];
rz(2.00804) q[1];
rz(-pi) q[2];
rz(-1.1169711) q[3];
sx q[3];
rz(-1.3012851) q[3];
sx q[3];
rz(-2.5170642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1699367) q[2];
sx q[2];
rz(-1.1145096) q[2];
sx q[2];
rz(1.8709095) q[2];
rz(3.0830248) q[3];
sx q[3];
rz(-1.4437557) q[3];
sx q[3];
rz(-0.34846714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.612065) q[0];
sx q[0];
rz(-2.6030774) q[0];
sx q[0];
rz(-0.2659604) q[0];
rz(0.82414857) q[1];
sx q[1];
rz(-1.8310841) q[1];
sx q[1];
rz(1.5305653) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96722051) q[0];
sx q[0];
rz(-1.8215113) q[0];
sx q[0];
rz(-2.870737) q[0];
rz(0.058892756) q[2];
sx q[2];
rz(-2.4097171) q[2];
sx q[2];
rz(2.1443387) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2648523) q[1];
sx q[1];
rz(-1.0050259) q[1];
sx q[1];
rz(-1.8840204) q[1];
rz(-pi) q[2];
rz(1.2208185) q[3];
sx q[3];
rz(-0.82982291) q[3];
sx q[3];
rz(-2.8184824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.019235762) q[2];
sx q[2];
rz(-1.1160206) q[2];
sx q[2];
rz(-0.96987152) q[2];
rz(0.14361778) q[3];
sx q[3];
rz(-0.93846455) q[3];
sx q[3];
rz(2.5041049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73760167) q[0];
sx q[0];
rz(-2.8856475) q[0];
sx q[0];
rz(3.0373489) q[0];
rz(0.741611) q[1];
sx q[1];
rz(-2.9545018) q[1];
sx q[1];
rz(-2.7033477) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61442626) q[0];
sx q[0];
rz(-2.2402555) q[0];
sx q[0];
rz(1.7015083) q[0];
x q[1];
rz(2.1684219) q[2];
sx q[2];
rz(-1.3750374) q[2];
sx q[2];
rz(2.6225066) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.91854535) q[1];
sx q[1];
rz(-2.1761736) q[1];
sx q[1];
rz(3.1001841) q[1];
rz(-pi) q[2];
rz(2.8986071) q[3];
sx q[3];
rz(-2.5930517) q[3];
sx q[3];
rz(1.0787971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7354342) q[2];
sx q[2];
rz(-0.27190748) q[2];
sx q[2];
rz(-0.3212277) q[2];
rz(-0.99212232) q[3];
sx q[3];
rz(-2.0480053) q[3];
sx q[3];
rz(1.4120302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2883437) q[0];
sx q[0];
rz(-3.0217331) q[0];
sx q[0];
rz(-2.993592) q[0];
rz(-2.0026228) q[1];
sx q[1];
rz(-0.6593467) q[1];
sx q[1];
rz(-2.151087) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73511998) q[0];
sx q[0];
rz(-0.48790259) q[0];
sx q[0];
rz(-0.16012971) q[0];
x q[1];
rz(-2.4100346) q[2];
sx q[2];
rz(-1.2148982) q[2];
sx q[2];
rz(1.4479835) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3618092) q[1];
sx q[1];
rz(-1.6020007) q[1];
sx q[1];
rz(3.1364016) q[1];
rz(-pi) q[2];
rz(-2.2049974) q[3];
sx q[3];
rz(-0.59118012) q[3];
sx q[3];
rz(1.180296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1691576) q[2];
sx q[2];
rz(-0.7908228) q[2];
sx q[2];
rz(0.54231826) q[2];
rz(-1.4646685) q[3];
sx q[3];
rz(-1.2277536) q[3];
sx q[3];
rz(0.98384583) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2391613) q[0];
sx q[0];
rz(-0.79817525) q[0];
sx q[0];
rz(0.27957988) q[0];
rz(2.3565049) q[1];
sx q[1];
rz(-2.3468192) q[1];
sx q[1];
rz(-0.40207544) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3668985) q[0];
sx q[0];
rz(-2.0089564) q[0];
sx q[0];
rz(-1.9209981) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3763972) q[2];
sx q[2];
rz(-2.1957261) q[2];
sx q[2];
rz(-0.13002061) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2759034) q[1];
sx q[1];
rz(-0.69178693) q[1];
sx q[1];
rz(-1.8744286) q[1];
rz(-2.8567186) q[3];
sx q[3];
rz(-1.9683629) q[3];
sx q[3];
rz(-1.3314825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1606007) q[2];
sx q[2];
rz(-2.4805785) q[2];
sx q[2];
rz(0.92657363) q[2];
rz(2.5911234) q[3];
sx q[3];
rz(-0.87369839) q[3];
sx q[3];
rz(1.2800823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67113039) q[0];
sx q[0];
rz(-1.067938) q[0];
sx q[0];
rz(-2.9965042) q[0];
rz(3.0628224) q[1];
sx q[1];
rz(-1.4046184) q[1];
sx q[1];
rz(1.2974993) q[1];
rz(2.3433122) q[2];
sx q[2];
rz(-1.5816806) q[2];
sx q[2];
rz(0.47368373) q[2];
rz(-2.0918431) q[3];
sx q[3];
rz(-1.4561984) q[3];
sx q[3];
rz(0.39312058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
