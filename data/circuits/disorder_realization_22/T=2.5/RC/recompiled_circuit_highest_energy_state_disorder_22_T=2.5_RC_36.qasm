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
rz(0.15777388) q[0];
sx q[0];
rz(-1.1717492) q[0];
sx q[0];
rz(1.2494614) q[0];
rz(-0.14021048) q[1];
sx q[1];
rz(-1.6970716) q[1];
sx q[1];
rz(0.061847774) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01067804) q[0];
sx q[0];
rz(-2.317551) q[0];
sx q[0];
rz(0.86573059) q[0];
rz(-pi) q[1];
rz(-1.5711478) q[2];
sx q[2];
rz(-1.6612189) q[2];
sx q[2];
rz(-3.0302559) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0272408) q[1];
sx q[1];
rz(-1.6682965) q[1];
sx q[1];
rz(0.84836324) q[1];
rz(2.2592945) q[3];
sx q[3];
rz(-1.2475444) q[3];
sx q[3];
rz(-2.0417449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7141815) q[2];
sx q[2];
rz(-2.2990871) q[2];
sx q[2];
rz(-2.8851435) q[2];
rz(-0.17476684) q[3];
sx q[3];
rz(-1.9758965) q[3];
sx q[3];
rz(-0.83785653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3278219) q[0];
sx q[0];
rz(-2.796266) q[0];
sx q[0];
rz(-0.12369618) q[0];
rz(-2.9028614) q[1];
sx q[1];
rz(-2.3614466) q[1];
sx q[1];
rz(3.0366268) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86228131) q[0];
sx q[0];
rz(-2.263592) q[0];
sx q[0];
rz(0.19578085) q[0];
rz(2.6506422) q[2];
sx q[2];
rz(-0.46451312) q[2];
sx q[2];
rz(2.2238942) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8387977) q[1];
sx q[1];
rz(-2.0640336) q[1];
sx q[1];
rz(-2.8431945) q[1];
x q[2];
rz(1.9267004) q[3];
sx q[3];
rz(-0.89733636) q[3];
sx q[3];
rz(1.2237024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3671942) q[2];
sx q[2];
rz(-1.4906733) q[2];
sx q[2];
rz(1.4614089) q[2];
rz(-1.8558308) q[3];
sx q[3];
rz(-1.6971842) q[3];
sx q[3];
rz(0.27209601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9937781) q[0];
sx q[0];
rz(-2.9312134) q[0];
sx q[0];
rz(2.126597) q[0];
rz(0.89730942) q[1];
sx q[1];
rz(-1.6474479) q[1];
sx q[1];
rz(-2.4651333) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7583026) q[0];
sx q[0];
rz(-1.5967909) q[0];
sx q[0];
rz(0.95375188) q[0];
rz(-pi) q[1];
rz(0.93134201) q[2];
sx q[2];
rz(-2.4081875) q[2];
sx q[2];
rz(-2.1426107) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9046717) q[1];
sx q[1];
rz(-2.7075012) q[1];
sx q[1];
rz(-1.0755811) q[1];
rz(-1.5900988) q[3];
sx q[3];
rz(-0.75452828) q[3];
sx q[3];
rz(-1.4135264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.93495381) q[2];
sx q[2];
rz(-2.166344) q[2];
sx q[2];
rz(-2.3812531) q[2];
rz(-1.9350516) q[3];
sx q[3];
rz(-0.81086603) q[3];
sx q[3];
rz(0.84347239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43561414) q[0];
sx q[0];
rz(-2.5137081) q[0];
sx q[0];
rz(-1.5950369) q[0];
rz(2.4226923) q[1];
sx q[1];
rz(-1.3201821) q[1];
sx q[1];
rz(0.23153201) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2255428) q[0];
sx q[0];
rz(-2.0839491) q[0];
sx q[0];
rz(1.0731359) q[0];
rz(2.3489352) q[2];
sx q[2];
rz(-1.5735627) q[2];
sx q[2];
rz(-2.1994928) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0147103) q[1];
sx q[1];
rz(-2.0538804) q[1];
sx q[1];
rz(1.4654069) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60482408) q[3];
sx q[3];
rz(-1.2422274) q[3];
sx q[3];
rz(-0.25024763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4270758) q[2];
sx q[2];
rz(-2.7526553) q[2];
sx q[2];
rz(0.60849774) q[2];
rz(-2.9454339) q[3];
sx q[3];
rz(-1.4213296) q[3];
sx q[3];
rz(-0.99854809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8683559) q[0];
sx q[0];
rz(-2.4899857) q[0];
sx q[0];
rz(-2.3625145) q[0];
rz(-0.92998663) q[1];
sx q[1];
rz(-1.168074) q[1];
sx q[1];
rz(-1.1735865) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0239651) q[0];
sx q[0];
rz(-1.5668586) q[0];
sx q[0];
rz(1.2327475) q[0];
x q[1];
rz(0.9554602) q[2];
sx q[2];
rz(-1.6809856) q[2];
sx q[2];
rz(-0.9166536) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.47347676) q[1];
sx q[1];
rz(-1.2537696) q[1];
sx q[1];
rz(0.067639785) q[1];
x q[2];
rz(-2.7435674) q[3];
sx q[3];
rz(-2.6809691) q[3];
sx q[3];
rz(1.4945791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.32213) q[2];
sx q[2];
rz(-0.16699114) q[2];
sx q[2];
rz(-2.4665311) q[2];
rz(0.86554646) q[3];
sx q[3];
rz(-1.5138488) q[3];
sx q[3];
rz(-1.9338098) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.753767) q[0];
sx q[0];
rz(-0.017711552) q[0];
sx q[0];
rz(-0.87131635) q[0];
rz(2.9246092) q[1];
sx q[1];
rz(-1.5996409) q[1];
sx q[1];
rz(1.3380231) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3335452) q[0];
sx q[0];
rz(-1.4879972) q[0];
sx q[0];
rz(2.4890578) q[0];
x q[1];
rz(-0.84514825) q[2];
sx q[2];
rz(-1.4652227) q[2];
sx q[2];
rz(2.0279864) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.034778193) q[1];
sx q[1];
rz(-2.9902774) q[1];
sx q[1];
rz(-1.2447912) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35049482) q[3];
sx q[3];
rz(-0.51096254) q[3];
sx q[3];
rz(3.0929589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1285105) q[2];
sx q[2];
rz(-1.4980114) q[2];
sx q[2];
rz(-0.80005542) q[2];
rz(0.99177805) q[3];
sx q[3];
rz(-2.1938775) q[3];
sx q[3];
rz(0.94314027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8868788) q[0];
sx q[0];
rz(-2.7006221) q[0];
sx q[0];
rz(-0.15175858) q[0];
rz(1.7249379) q[1];
sx q[1];
rz(-0.85985008) q[1];
sx q[1];
rz(-0.83622611) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4164616) q[0];
sx q[0];
rz(-0.90354474) q[0];
sx q[0];
rz(-2.6362772) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.012846666) q[2];
sx q[2];
rz(-0.57672665) q[2];
sx q[2];
rz(1.4801413) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9735896) q[1];
sx q[1];
rz(-2.8781946) q[1];
sx q[1];
rz(1.2348639) q[1];
x q[2];
rz(2.8854675) q[3];
sx q[3];
rz(-2.6626592) q[3];
sx q[3];
rz(1.1161982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6553216) q[2];
sx q[2];
rz(-0.063491193) q[2];
sx q[2];
rz(-3.0625694) q[2];
rz(0.71845734) q[3];
sx q[3];
rz(-1.6301165) q[3];
sx q[3];
rz(-2.2169936) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.786161) q[0];
sx q[0];
rz(-2.5434255) q[0];
sx q[0];
rz(-0.86791903) q[0];
rz(-0.68880853) q[1];
sx q[1];
rz(-0.30615607) q[1];
sx q[1];
rz(-0.93592962) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42396069) q[0];
sx q[0];
rz(-2.0905295) q[0];
sx q[0];
rz(-1.7802159) q[0];
x q[1];
rz(-0.34078074) q[2];
sx q[2];
rz(-1.9592957) q[2];
sx q[2];
rz(-2.4430371) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1327428) q[1];
sx q[1];
rz(-0.32948438) q[1];
sx q[1];
rz(-1.4047755) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3981188) q[3];
sx q[3];
rz(-1.5718565) q[3];
sx q[3];
rz(0.27976945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9800637) q[2];
sx q[2];
rz(-0.65902013) q[2];
sx q[2];
rz(2.8847983) q[2];
rz(-2.1234546) q[3];
sx q[3];
rz(-1.1192106) q[3];
sx q[3];
rz(-2.170678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7629906) q[0];
sx q[0];
rz(-2.584223) q[0];
sx q[0];
rz(-0.85365224) q[0];
rz(1.8866106) q[1];
sx q[1];
rz(-1.1458784) q[1];
sx q[1];
rz(2.8010211) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84682631) q[0];
sx q[0];
rz(-1.7309523) q[0];
sx q[0];
rz(-1.7677884) q[0];
rz(0.30605028) q[2];
sx q[2];
rz(-0.69717625) q[2];
sx q[2];
rz(2.9070284) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7813959) q[1];
sx q[1];
rz(-2.2193647) q[1];
sx q[1];
rz(0.99332033) q[1];
x q[2];
rz(-2.0331618) q[3];
sx q[3];
rz(-1.8183072) q[3];
sx q[3];
rz(2.4570297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.42369947) q[2];
sx q[2];
rz(-1.7607949) q[2];
sx q[2];
rz(-0.37810668) q[2];
rz(-2.9041491) q[3];
sx q[3];
rz(-2.0707371) q[3];
sx q[3];
rz(2.8606991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.045192748) q[0];
sx q[0];
rz(-1.0843596) q[0];
sx q[0];
rz(0.64724809) q[0];
rz(1.1680394) q[1];
sx q[1];
rz(-2.2868575) q[1];
sx q[1];
rz(0.38356575) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7591568) q[0];
sx q[0];
rz(-1.54017) q[0];
sx q[0];
rz(-1.5923862) q[0];
rz(0.53159376) q[2];
sx q[2];
rz(-1.2832912) q[2];
sx q[2];
rz(0.93898857) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3473222) q[1];
sx q[1];
rz(-1.3065265) q[1];
sx q[1];
rz(1.3572973) q[1];
rz(-pi) q[2];
rz(-2.0862574) q[3];
sx q[3];
rz(-1.664242) q[3];
sx q[3];
rz(-2.6938113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.26455227) q[2];
sx q[2];
rz(-2.2525747) q[2];
sx q[2];
rz(-1.5569347) q[2];
rz(1.5133739) q[3];
sx q[3];
rz(-1.3686562) q[3];
sx q[3];
rz(-0.42678601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7814816) q[0];
sx q[0];
rz(-1.7088912) q[0];
sx q[0];
rz(0.29722469) q[0];
rz(-1.9602641) q[1];
sx q[1];
rz(-1.344463) q[1];
sx q[1];
rz(-0.32483473) q[1];
rz(-1.1245341) q[2];
sx q[2];
rz(-2.6206803) q[2];
sx q[2];
rz(-0.16823106) q[2];
rz(-0.060130974) q[3];
sx q[3];
rz(-2.6562666) q[3];
sx q[3];
rz(-0.2851214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
