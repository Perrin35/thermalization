OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2774529) q[0];
sx q[0];
rz(-1.5885408) q[0];
sx q[0];
rz(1.5074402) q[0];
rz(1.5965257) q[1];
sx q[1];
rz(-0.59626055) q[1];
sx q[1];
rz(-2.526386) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98696729) q[0];
sx q[0];
rz(-1.1095424) q[0];
sx q[0];
rz(2.2485562) q[0];
rz(-pi) q[1];
rz(2.6745559) q[2];
sx q[2];
rz(-0.39614284) q[2];
sx q[2];
rz(0.064183891) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.93278904) q[1];
sx q[1];
rz(-2.1425769) q[1];
sx q[1];
rz(-0.51704452) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7513566) q[3];
sx q[3];
rz(-1.2689586) q[3];
sx q[3];
rz(-0.80871049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.4404099) q[2];
sx q[2];
rz(-1.6117233) q[2];
sx q[2];
rz(0.33828503) q[2];
rz(-1.7017378) q[3];
sx q[3];
rz(-2.2262636) q[3];
sx q[3];
rz(-2.2556944) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97025362) q[0];
sx q[0];
rz(-0.71115029) q[0];
sx q[0];
rz(3.1112444) q[0];
rz(0.066210315) q[1];
sx q[1];
rz(-0.98774424) q[1];
sx q[1];
rz(-1.5240086) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.671316) q[0];
sx q[0];
rz(-2.94194) q[0];
sx q[0];
rz(1.5623564) q[0];
x q[1];
rz(-1.5288058) q[2];
sx q[2];
rz(-2.6868372) q[2];
sx q[2];
rz(0.08200478) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.46088947) q[1];
sx q[1];
rz(-1.5017121) q[1];
sx q[1];
rz(0.99980385) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3832983) q[3];
sx q[3];
rz(-1.4783995) q[3];
sx q[3];
rz(0.20860162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7559738) q[2];
sx q[2];
rz(-0.9884584) q[2];
sx q[2];
rz(1.1478109) q[2];
rz(1.3267481) q[3];
sx q[3];
rz(-1.3245405) q[3];
sx q[3];
rz(-0.23708788) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0148934) q[0];
sx q[0];
rz(-0.48148695) q[0];
sx q[0];
rz(-2.8258064) q[0];
rz(-2.2029927) q[1];
sx q[1];
rz(-1.4626075) q[1];
sx q[1];
rz(2.8895203) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2825851) q[0];
sx q[0];
rz(-2.1827645) q[0];
sx q[0];
rz(0.99143272) q[0];
rz(-pi) q[1];
rz(2.4263072) q[2];
sx q[2];
rz(-0.89865696) q[2];
sx q[2];
rz(0.31179024) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7604916) q[1];
sx q[1];
rz(-2.4931506) q[1];
sx q[1];
rz(1.8849424) q[1];
rz(-pi) q[2];
rz(-0.73022233) q[3];
sx q[3];
rz(-0.50958868) q[3];
sx q[3];
rz(-2.2024221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0198274) q[2];
sx q[2];
rz(-0.44744197) q[2];
sx q[2];
rz(-3.1075409) q[2];
rz(0.017459067) q[3];
sx q[3];
rz(-1.7826467) q[3];
sx q[3];
rz(-1.0954558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62717342) q[0];
sx q[0];
rz(-1.592941) q[0];
sx q[0];
rz(-1.6148286) q[0];
rz(-1.0871672) q[1];
sx q[1];
rz(-0.68030578) q[1];
sx q[1];
rz(0.70708752) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2061545) q[0];
sx q[0];
rz(-1.1797138) q[0];
sx q[0];
rz(-2.6576256) q[0];
x q[1];
rz(2.884042) q[2];
sx q[2];
rz(-0.45986816) q[2];
sx q[2];
rz(-2.3515153) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6882119) q[1];
sx q[1];
rz(-2.4773295) q[1];
sx q[1];
rz(1.3761671) q[1];
rz(-pi) q[2];
rz(-2.1726923) q[3];
sx q[3];
rz(-2.7728191) q[3];
sx q[3];
rz(0.2975279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2924071) q[2];
sx q[2];
rz(-1.1881928) q[2];
sx q[2];
rz(-0.46009955) q[2];
rz(-1.7442616) q[3];
sx q[3];
rz(-1.5528691) q[3];
sx q[3];
rz(0.24266711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8206772) q[0];
sx q[0];
rz(-2.9292332) q[0];
sx q[0];
rz(-1.3943577) q[0];
rz(-1.0955411) q[1];
sx q[1];
rz(-1.54116) q[1];
sx q[1];
rz(-2.8869693) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1063227) q[0];
sx q[0];
rz(-1.584504) q[0];
sx q[0];
rz(-1.6499707) q[0];
rz(1.5484372) q[2];
sx q[2];
rz(-2.5639113) q[2];
sx q[2];
rz(2.3235842) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.188272) q[1];
sx q[1];
rz(-2.4362262) q[1];
sx q[1];
rz(2.7642247) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7084064) q[3];
sx q[3];
rz(-2.232589) q[3];
sx q[3];
rz(1.8720686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.022481) q[2];
sx q[2];
rz(-2.9412061) q[2];
sx q[2];
rz(1.7648034) q[2];
rz(1.4962176) q[3];
sx q[3];
rz(-1.508537) q[3];
sx q[3];
rz(-2.0549324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6901533) q[0];
sx q[0];
rz(-0.7398766) q[0];
sx q[0];
rz(2.8421463) q[0];
rz(-1.0401789) q[1];
sx q[1];
rz(-1.4458011) q[1];
sx q[1];
rz(-0.20656955) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5305938) q[0];
sx q[0];
rz(-1.8078513) q[0];
sx q[0];
rz(0.80942746) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.39482306) q[2];
sx q[2];
rz(-2.0645112) q[2];
sx q[2];
rz(2.0815108) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5094604) q[1];
sx q[1];
rz(-1.3941947) q[1];
sx q[1];
rz(-0.86061865) q[1];
rz(0.30287403) q[3];
sx q[3];
rz(-2.2040963) q[3];
sx q[3];
rz(2.5686222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2999337) q[2];
sx q[2];
rz(-1.4175697) q[2];
sx q[2];
rz(0.28277961) q[2];
rz(-0.81280604) q[3];
sx q[3];
rz(-0.41906425) q[3];
sx q[3];
rz(-2.9747484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.92641002) q[0];
sx q[0];
rz(-2.0880501) q[0];
sx q[0];
rz(-2.7600631) q[0];
rz(2.5577257) q[1];
sx q[1];
rz(-0.54324141) q[1];
sx q[1];
rz(1.8136224) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5663985) q[0];
sx q[0];
rz(-1.1325206) q[0];
sx q[0];
rz(3.0153494) q[0];
rz(-pi) q[1];
rz(2.1778657) q[2];
sx q[2];
rz(-2.4237195) q[2];
sx q[2];
rz(1.8075862) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.59203397) q[1];
sx q[1];
rz(-2.4394819) q[1];
sx q[1];
rz(1.3105427) q[1];
rz(-1.9973203) q[3];
sx q[3];
rz(-1.0970308) q[3];
sx q[3];
rz(1.9667448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5232089) q[2];
sx q[2];
rz(-2.0760459) q[2];
sx q[2];
rz(-1.3605114) q[2];
rz(1.7112188) q[3];
sx q[3];
rz(-1.0083219) q[3];
sx q[3];
rz(-0.84806228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9946063) q[0];
sx q[0];
rz(-1.1598347) q[0];
sx q[0];
rz(-2.9597136) q[0];
rz(0.47422844) q[1];
sx q[1];
rz(-2.1209746) q[1];
sx q[1];
rz(-2.1906733) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5247105) q[0];
sx q[0];
rz(-1.8572154) q[0];
sx q[0];
rz(-0.25047238) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3046706) q[2];
sx q[2];
rz(-1.5106491) q[2];
sx q[2];
rz(-1.6595449) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.37590313) q[1];
sx q[1];
rz(-1.6357058) q[1];
sx q[1];
rz(-2.9325571) q[1];
x q[2];
rz(1.9740281) q[3];
sx q[3];
rz(-1.9208761) q[3];
sx q[3];
rz(-2.0389327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9178847) q[2];
sx q[2];
rz(-0.48030883) q[2];
sx q[2];
rz(0.075909464) q[2];
rz(0.54801303) q[3];
sx q[3];
rz(-1.8173822) q[3];
sx q[3];
rz(-2.5089335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6178745) q[0];
sx q[0];
rz(-2.0848367) q[0];
sx q[0];
rz(1.7653718) q[0];
rz(2.7245522) q[1];
sx q[1];
rz(-1.4191671) q[1];
sx q[1];
rz(2.4818647) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10044554) q[0];
sx q[0];
rz(-2.1491096) q[0];
sx q[0];
rz(-2.5006177) q[0];
rz(-2.9218036) q[2];
sx q[2];
rz(-0.79384365) q[2];
sx q[2];
rz(-1.1500051) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1253423) q[1];
sx q[1];
rz(-0.94031912) q[1];
sx q[1];
rz(-2.7793105) q[1];
rz(-pi) q[2];
rz(-3.0845853) q[3];
sx q[3];
rz(-2.1087286) q[3];
sx q[3];
rz(1.2008592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62347162) q[2];
sx q[2];
rz(-2.3770964) q[2];
sx q[2];
rz(-1.0260322) q[2];
rz(2.9927411) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(0.025645105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24348564) q[0];
sx q[0];
rz(-2.4004816) q[0];
sx q[0];
rz(1.3056668) q[0];
rz(1.1765515) q[1];
sx q[1];
rz(-1.2780317) q[1];
sx q[1];
rz(2.1059039) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4982088) q[0];
sx q[0];
rz(-1.1632825) q[0];
sx q[0];
rz(-1.9252752) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87558526) q[2];
sx q[2];
rz(-1.1165809) q[2];
sx q[2];
rz(0.21739612) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.021691572) q[1];
sx q[1];
rz(-2.7530115) q[1];
sx q[1];
rz(-2.4619224) q[1];
x q[2];
rz(1.3003179) q[3];
sx q[3];
rz(-0.83993739) q[3];
sx q[3];
rz(-0.65281103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.62853652) q[2];
sx q[2];
rz(-2.3302902) q[2];
sx q[2];
rz(1.9899842) q[2];
rz(2.628905) q[3];
sx q[3];
rz(-2.0435464) q[3];
sx q[3];
rz(-0.36469665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7619027) q[0];
sx q[0];
rz(-1.7871478) q[0];
sx q[0];
rz(-2.3085069) q[0];
rz(1.5079386) q[1];
sx q[1];
rz(-0.58273756) q[1];
sx q[1];
rz(2.6599463) q[1];
rz(-0.41082906) q[2];
sx q[2];
rz(-1.5317691) q[2];
sx q[2];
rz(-0.017824235) q[2];
rz(-3.1291943) q[3];
sx q[3];
rz(-2.5231902) q[3];
sx q[3];
rz(-1.4004422) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
