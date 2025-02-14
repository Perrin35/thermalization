OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.088012785) q[0];
sx q[0];
rz(3.454257) q[0];
sx q[0];
rz(9.5587048) q[0];
rz(3.1349831) q[1];
sx q[1];
rz(-1.5899038) q[1];
sx q[1];
rz(2.7166727) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44533396) q[0];
sx q[0];
rz(-0.15275341) q[0];
sx q[0];
rz(-0.85102083) q[0];
rz(-pi) q[1];
rz(-1.0966461) q[2];
sx q[2];
rz(-0.72410781) q[2];
sx q[2];
rz(0.5257789) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.9953138) q[1];
sx q[1];
rz(-1.9094719) q[1];
sx q[1];
rz(2.2168753) q[1];
x q[2];
rz(-2.4272301) q[3];
sx q[3];
rz(-2.2788725) q[3];
sx q[3];
rz(1.8797415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6801844) q[2];
sx q[2];
rz(-0.37698656) q[2];
sx q[2];
rz(0.11412966) q[2];
rz(-2.7528609) q[3];
sx q[3];
rz(-0.26624334) q[3];
sx q[3];
rz(1.8069327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-0.0085793063) q[0];
sx q[0];
rz(-0.40224922) q[0];
sx q[0];
rz(3.114793) q[0];
rz(1.9642448) q[1];
sx q[1];
rz(-0.64759308) q[1];
sx q[1];
rz(3.1039544) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5863754) q[0];
sx q[0];
rz(-0.77459413) q[0];
sx q[0];
rz(3.1223749) q[0];
rz(1.3635719) q[2];
sx q[2];
rz(-1.3517016) q[2];
sx q[2];
rz(-1.6082282) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.37477949) q[1];
sx q[1];
rz(-1.0008924) q[1];
sx q[1];
rz(0.5809231) q[1];
x q[2];
rz(0.82017558) q[3];
sx q[3];
rz(-1.7595152) q[3];
sx q[3];
rz(2.3403326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2521952) q[2];
sx q[2];
rz(-2.1571721) q[2];
sx q[2];
rz(-0.21749116) q[2];
rz(-2.7018231) q[3];
sx q[3];
rz(-1.7871126) q[3];
sx q[3];
rz(0.10194889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(0.52461034) q[0];
sx q[0];
rz(-3.1341902) q[0];
sx q[0];
rz(0.55968416) q[0];
rz(0.89645487) q[1];
sx q[1];
rz(-0.47210109) q[1];
sx q[1];
rz(-0.40759531) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2893163) q[0];
sx q[0];
rz(-1.455515) q[0];
sx q[0];
rz(0.98457054) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0733068) q[2];
sx q[2];
rz(-0.90203055) q[2];
sx q[2];
rz(-2.1055773) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7771908) q[1];
sx q[1];
rz(-2.6620416) q[1];
sx q[1];
rz(-0.24141356) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0225917) q[3];
sx q[3];
rz(-1.8361409) q[3];
sx q[3];
rz(3.1309553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.018230351) q[2];
sx q[2];
rz(-1.9841649) q[2];
sx q[2];
rz(2.1165712) q[2];
rz(1.7068663) q[3];
sx q[3];
rz(-2.6985109) q[3];
sx q[3];
rz(-2.4923435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(0.0593112) q[0];
sx q[0];
rz(-2.870443) q[0];
sx q[0];
rz(-2.6442288) q[0];
rz(1.8479895) q[1];
sx q[1];
rz(-1.0062287) q[1];
sx q[1];
rz(1.2327548) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13834223) q[0];
sx q[0];
rz(-1.3673393) q[0];
sx q[0];
rz(0.5749216) q[0];
rz(-pi) q[1];
rz(1.1511939) q[2];
sx q[2];
rz(-0.91165724) q[2];
sx q[2];
rz(1.4105547) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2805505) q[1];
sx q[1];
rz(-0.32233626) q[1];
sx q[1];
rz(-0.68962421) q[1];
rz(-pi) q[2];
x q[2];
rz(2.74128) q[3];
sx q[3];
rz(-1.5535206) q[3];
sx q[3];
rz(-3.0028385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5645912) q[2];
sx q[2];
rz(-0.61508721) q[2];
sx q[2];
rz(0.1423398) q[2];
rz(-1.1764935) q[3];
sx q[3];
rz(-2.1409972) q[3];
sx q[3];
rz(0.4308027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7560526) q[0];
sx q[0];
rz(-0.98243326) q[0];
sx q[0];
rz(-0.2734215) q[0];
rz(2.7418819) q[1];
sx q[1];
rz(-0.81396657) q[1];
sx q[1];
rz(-2.8846557) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98824153) q[0];
sx q[0];
rz(-2.0978161) q[0];
sx q[0];
rz(-0.37470438) q[0];
x q[1];
rz(0.24946282) q[2];
sx q[2];
rz(-1.9469065) q[2];
sx q[2];
rz(0.94406908) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7796799) q[1];
sx q[1];
rz(-3.0579902) q[1];
sx q[1];
rz(2.1474775) q[1];
x q[2];
rz(2.8572548) q[3];
sx q[3];
rz(-0.91129843) q[3];
sx q[3];
rz(-2.4267765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5895245) q[2];
sx q[2];
rz(-2.7858211) q[2];
sx q[2];
rz(-0.75999981) q[2];
rz(-2.149557) q[3];
sx q[3];
rz(-1.9325247) q[3];
sx q[3];
rz(1.1442643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8511667) q[0];
sx q[0];
rz(-2.5100584) q[0];
sx q[0];
rz(-2.5354711) q[0];
rz(2.0947314) q[1];
sx q[1];
rz(-1.650834) q[1];
sx q[1];
rz(-0.077233888) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.58409) q[0];
sx q[0];
rz(-1.5833921) q[0];
sx q[0];
rz(-1.6112616) q[0];
x q[1];
rz(2.8229643) q[2];
sx q[2];
rz(-1.1201522) q[2];
sx q[2];
rz(-0.81540996) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8171766) q[1];
sx q[1];
rz(-2.5253173) q[1];
sx q[1];
rz(-2.3072427) q[1];
x q[2];
rz(0.24155946) q[3];
sx q[3];
rz(-0.67725855) q[3];
sx q[3];
rz(0.23158555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.91654009) q[2];
sx q[2];
rz(-0.70987916) q[2];
sx q[2];
rz(-0.27099657) q[2];
rz(-2.5535876) q[3];
sx q[3];
rz(-2.3376412) q[3];
sx q[3];
rz(2.7325381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3410909) q[0];
sx q[0];
rz(-1.5301457) q[0];
sx q[0];
rz(-2.8588168) q[0];
rz(2.458789) q[1];
sx q[1];
rz(-0.86137259) q[1];
sx q[1];
rz(0.85321325) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5763007) q[0];
sx q[0];
rz(-1.9423663) q[0];
sx q[0];
rz(-2.316733) q[0];
x q[1];
rz(0.81099895) q[2];
sx q[2];
rz(-1.2110707) q[2];
sx q[2];
rz(-0.0051509858) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4567199) q[1];
sx q[1];
rz(-1.8032852) q[1];
sx q[1];
rz(-2.9732735) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0045347) q[3];
sx q[3];
rz(-0.2314724) q[3];
sx q[3];
rz(1.2616518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0736488) q[2];
sx q[2];
rz(-0.51764071) q[2];
sx q[2];
rz(-2.850387) q[2];
rz(-1.9368885) q[3];
sx q[3];
rz(-0.57345814) q[3];
sx q[3];
rz(-1.5709491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7912927) q[0];
sx q[0];
rz(-1.6181823) q[0];
sx q[0];
rz(0.6268025) q[0];
rz(1.0621915) q[1];
sx q[1];
rz(-0.79963446) q[1];
sx q[1];
rz(-2.3041384) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0499247) q[0];
sx q[0];
rz(-0.44418884) q[0];
sx q[0];
rz(1.8416406) q[0];
x q[1];
rz(-2.8026514) q[2];
sx q[2];
rz(-1.7941448) q[2];
sx q[2];
rz(-1.3517584) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7211322) q[1];
sx q[1];
rz(-2.0803806) q[1];
sx q[1];
rz(2.9886118) q[1];
x q[2];
rz(-1.4222048) q[3];
sx q[3];
rz(-1.0576384) q[3];
sx q[3];
rz(-0.10529127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8219139) q[2];
sx q[2];
rz(-1.2205114) q[2];
sx q[2];
rz(2.1004045) q[2];
rz(-0.90211165) q[3];
sx q[3];
rz(-1.9769042) q[3];
sx q[3];
rz(-2.5743217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.64004552) q[0];
sx q[0];
rz(-0.29842672) q[0];
sx q[0];
rz(3.0409467) q[0];
rz(1.3238662) q[1];
sx q[1];
rz(-2.3034425) q[1];
sx q[1];
rz(-2.5460338) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57306615) q[0];
sx q[0];
rz(-1.0254481) q[0];
sx q[0];
rz(-2.5667356) q[0];
rz(-2.9352292) q[2];
sx q[2];
rz(-2.2225755) q[2];
sx q[2];
rz(0.18557063) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7337949) q[1];
sx q[1];
rz(-0.76939728) q[1];
sx q[1];
rz(0.35843884) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15138126) q[3];
sx q[3];
rz(-2.3321242) q[3];
sx q[3];
rz(2.4224506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9969534) q[2];
sx q[2];
rz(-2.3914631) q[2];
sx q[2];
rz(-1.7993125) q[2];
rz(0.73623776) q[3];
sx q[3];
rz(-2.8056371) q[3];
sx q[3];
rz(3.0430702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1243684) q[0];
sx q[0];
rz(-2.4898744) q[0];
sx q[0];
rz(0.29618725) q[0];
rz(1.724297) q[1];
sx q[1];
rz(-2.2139151) q[1];
sx q[1];
rz(-0.16253026) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4341287) q[0];
sx q[0];
rz(-1.5445827) q[0];
sx q[0];
rz(-3.1349584) q[0];
rz(-pi) q[1];
rz(-2.7835566) q[2];
sx q[2];
rz(-1.0834143) q[2];
sx q[2];
rz(2.0646937) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1028984) q[1];
sx q[1];
rz(-1.4331281) q[1];
sx q[1];
rz(-0.18935151) q[1];
rz(-pi) q[2];
rz(-0.53402114) q[3];
sx q[3];
rz(-2.0724618) q[3];
sx q[3];
rz(2.925774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.51809239) q[2];
sx q[2];
rz(-1.9237498) q[2];
sx q[2];
rz(0.45563844) q[2];
rz(-1.0738922) q[3];
sx q[3];
rz(-0.62993252) q[3];
sx q[3];
rz(-0.58551252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0940654) q[0];
sx q[0];
rz(-1.3205262) q[0];
sx q[0];
rz(2.5393215) q[0];
rz(0.31996721) q[1];
sx q[1];
rz(-1.1155557) q[1];
sx q[1];
rz(-1.3394578) q[1];
rz(-0.336774) q[2];
sx q[2];
rz(-1.5143186) q[2];
sx q[2];
rz(1.0568525) q[2];
rz(1.6781758) q[3];
sx q[3];
rz(-0.79859514) q[3];
sx q[3];
rz(2.7664281) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
