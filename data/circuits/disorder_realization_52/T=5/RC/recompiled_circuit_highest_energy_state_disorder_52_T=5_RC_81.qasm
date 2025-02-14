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
rz(1.7088543) q[0];
sx q[0];
rz(-2.813485) q[0];
sx q[0];
rz(2.3018667) q[0];
rz(1.4108763) q[1];
sx q[1];
rz(-1.5411935) q[1];
sx q[1];
rz(-2.3720001) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1043848) q[0];
sx q[0];
rz(-1.6055371) q[0];
sx q[0];
rz(1.6503235) q[0];
rz(-pi) q[1];
rz(1.5985477) q[2];
sx q[2];
rz(-2.395438) q[2];
sx q[2];
rz(2.1389824) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.855866) q[1];
sx q[1];
rz(-1.2190423) q[1];
sx q[1];
rz(-3.1409164) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64483492) q[3];
sx q[3];
rz(-1.1434237) q[3];
sx q[3];
rz(-2.4348034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8771693) q[2];
sx q[2];
rz(-1.2494272) q[2];
sx q[2];
rz(-1.1733615) q[2];
rz(2.7133283) q[3];
sx q[3];
rz(-2.5731125) q[3];
sx q[3];
rz(2.4885524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0083171) q[0];
sx q[0];
rz(-2.7011217) q[0];
sx q[0];
rz(1.0622729) q[0];
rz(1.5509037) q[1];
sx q[1];
rz(-2.5804602) q[1];
sx q[1];
rz(-1.0947469) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0355201) q[0];
sx q[0];
rz(-2.8815334) q[0];
sx q[0];
rz(1.3706657) q[0];
rz(2.6954912) q[2];
sx q[2];
rz(-0.94233905) q[2];
sx q[2];
rz(2.7035509) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.96910406) q[1];
sx q[1];
rz(-1.8575694) q[1];
sx q[1];
rz(2.7412834) q[1];
x q[2];
rz(2.9238104) q[3];
sx q[3];
rz(-2.3541303) q[3];
sx q[3];
rz(1.2125208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4636479) q[2];
sx q[2];
rz(-2.8539168) q[2];
sx q[2];
rz(-2.2410683) q[2];
rz(1.0362222) q[3];
sx q[3];
rz(-0.13768727) q[3];
sx q[3];
rz(3.1305967) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66505945) q[0];
sx q[0];
rz(-2.0158975) q[0];
sx q[0];
rz(-1.6735459) q[0];
rz(-2.2131069) q[1];
sx q[1];
rz(-2.1378345) q[1];
sx q[1];
rz(1.4220062) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4741459) q[0];
sx q[0];
rz(-2.8774539) q[0];
sx q[0];
rz(-0.35281065) q[0];
x q[1];
rz(-2.456383) q[2];
sx q[2];
rz(-2.1985719) q[2];
sx q[2];
rz(2.1653099) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3513438) q[1];
sx q[1];
rz(-0.66748744) q[1];
sx q[1];
rz(2.5363467) q[1];
rz(0.4246131) q[3];
sx q[3];
rz(-1.8311005) q[3];
sx q[3];
rz(-1.2411959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15098393) q[2];
sx q[2];
rz(-1.432212) q[2];
sx q[2];
rz(2.3217311) q[2];
rz(-0.12623434) q[3];
sx q[3];
rz(-1.9641967) q[3];
sx q[3];
rz(-1.646515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6821297) q[0];
sx q[0];
rz(-1.7403025) q[0];
sx q[0];
rz(-1.6023741) q[0];
rz(2.050121) q[1];
sx q[1];
rz(-2.3075054) q[1];
sx q[1];
rz(1.9713255) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0905492) q[0];
sx q[0];
rz(-1.1950332) q[0];
sx q[0];
rz(-2.7997478) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2762047) q[2];
sx q[2];
rz(-2.3995993) q[2];
sx q[2];
rz(-0.8823673) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2728426) q[1];
sx q[1];
rz(-0.12050546) q[1];
sx q[1];
rz(-2.4122448) q[1];
rz(-1.3552731) q[3];
sx q[3];
rz(-1.2053262) q[3];
sx q[3];
rz(0.48367031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9200661) q[2];
sx q[2];
rz(-0.92851323) q[2];
sx q[2];
rz(-0.083219223) q[2];
rz(-0.65256882) q[3];
sx q[3];
rz(-0.021952732) q[3];
sx q[3];
rz(0.86161247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6457152) q[0];
sx q[0];
rz(-2.8577514) q[0];
sx q[0];
rz(-0.19677095) q[0];
rz(1.1605877) q[1];
sx q[1];
rz(-2.6996758) q[1];
sx q[1];
rz(1.7534076) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8179479) q[0];
sx q[0];
rz(-2.9716688) q[0];
sx q[0];
rz(-0.0090881149) q[0];
x q[1];
rz(1.722253) q[2];
sx q[2];
rz(-1.3706319) q[2];
sx q[2];
rz(1.4520926) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7418688) q[1];
sx q[1];
rz(-2.464964) q[1];
sx q[1];
rz(-1.2938188) q[1];
x q[2];
rz(-2.7775473) q[3];
sx q[3];
rz(-0.87556404) q[3];
sx q[3];
rz(-0.95330584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5768726) q[2];
sx q[2];
rz(-1.2083961) q[2];
sx q[2];
rz(1.8604856) q[2];
rz(-2.7992904) q[3];
sx q[3];
rz(-3.0908995) q[3];
sx q[3];
rz(-2.2127693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36011919) q[0];
sx q[0];
rz(-2.0442648) q[0];
sx q[0];
rz(-2.1088364) q[0];
rz(-2.4646941) q[1];
sx q[1];
rz(-2.758226) q[1];
sx q[1];
rz(-2.3597609) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6581377) q[0];
sx q[0];
rz(-1.0080604) q[0];
sx q[0];
rz(-1.9938009) q[0];
rz(-pi) q[1];
rz(1.5196477) q[2];
sx q[2];
rz(-1.6031853) q[2];
sx q[2];
rz(0.50078228) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5128153) q[1];
sx q[1];
rz(-2.1330962) q[1];
sx q[1];
rz(-1.2059216) q[1];
rz(-pi) q[2];
rz(0.12198066) q[3];
sx q[3];
rz(-2.8013419) q[3];
sx q[3];
rz(1.5911136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3013762) q[2];
sx q[2];
rz(-0.8679114) q[2];
sx q[2];
rz(-2.2678579) q[2];
rz(0.78911632) q[3];
sx q[3];
rz(-1.5849761) q[3];
sx q[3];
rz(-2.6553787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18272045) q[0];
sx q[0];
rz(-0.046165753) q[0];
sx q[0];
rz(0.090855457) q[0];
rz(-0.60984045) q[1];
sx q[1];
rz(-1.5900541) q[1];
sx q[1];
rz(-0.59744936) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0598568) q[0];
sx q[0];
rz(-1.6993521) q[0];
sx q[0];
rz(-0.14692326) q[0];
rz(-pi) q[1];
rz(0.84709527) q[2];
sx q[2];
rz(-2.6604746) q[2];
sx q[2];
rz(-1.4916949) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.3143464) q[1];
sx q[1];
rz(-2.0318077) q[1];
sx q[1];
rz(2.6516312) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1037056) q[3];
sx q[3];
rz(-1.9901313) q[3];
sx q[3];
rz(-1.208272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7647543) q[2];
sx q[2];
rz(-0.54543269) q[2];
sx q[2];
rz(3.0010014) q[2];
rz(1.8600672) q[3];
sx q[3];
rz(-1.8257273) q[3];
sx q[3];
rz(-2.6144821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.918688) q[0];
sx q[0];
rz(-1.0162901) q[0];
sx q[0];
rz(-1.6131529) q[0];
rz(0.81360045) q[1];
sx q[1];
rz(-3.0366812) q[1];
sx q[1];
rz(-2.334107) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4190237) q[0];
sx q[0];
rz(-1.3815855) q[0];
sx q[0];
rz(-1.7721227) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1296223) q[2];
sx q[2];
rz(-2.0106914) q[2];
sx q[2];
rz(2.9578046) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.85881104) q[1];
sx q[1];
rz(-2.3859897) q[1];
sx q[1];
rz(1.0048466) q[1];
rz(-pi) q[2];
rz(2.1106312) q[3];
sx q[3];
rz(-1.9412517) q[3];
sx q[3];
rz(-2.3453804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7878824) q[2];
sx q[2];
rz(-1.7870125) q[2];
sx q[2];
rz(-0.14459571) q[2];
rz(-1.1041798) q[3];
sx q[3];
rz(-2.927533) q[3];
sx q[3];
rz(2.1000699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5398194) q[0];
sx q[0];
rz(-1.1310534) q[0];
sx q[0];
rz(0.55651504) q[0];
rz(-2.3717608) q[1];
sx q[1];
rz(-0.12431215) q[1];
sx q[1];
rz(2.6803023) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71660626) q[0];
sx q[0];
rz(-1.5231774) q[0];
sx q[0];
rz(-0.99999962) q[0];
x q[1];
rz(2.2940852) q[2];
sx q[2];
rz(-0.88107785) q[2];
sx q[2];
rz(-2.412852) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.76595054) q[1];
sx q[1];
rz(-1.8643171) q[1];
sx q[1];
rz(-0.75325812) q[1];
x q[2];
rz(1.1848524) q[3];
sx q[3];
rz(-2.1084774) q[3];
sx q[3];
rz(-0.73962921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7028659) q[2];
sx q[2];
rz(-2.8361969) q[2];
sx q[2];
rz(2.3766282) q[2];
rz(-1.2139828) q[3];
sx q[3];
rz(-0.58911222) q[3];
sx q[3];
rz(0.37863076) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42534378) q[0];
sx q[0];
rz(-3.0916164) q[0];
sx q[0];
rz(-0.36323994) q[0];
rz(-1.4092457) q[1];
sx q[1];
rz(-2.1556518) q[1];
sx q[1];
rz(0.22663103) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5674681) q[0];
sx q[0];
rz(-0.13524817) q[0];
sx q[0];
rz(3.1194206) q[0];
rz(2.2025467) q[2];
sx q[2];
rz(-2.2151196) q[2];
sx q[2];
rz(-2.7959888) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5635101) q[1];
sx q[1];
rz(-0.77111608) q[1];
sx q[1];
rz(-1.095849) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5814621) q[3];
sx q[3];
rz(-1.9086325) q[3];
sx q[3];
rz(-0.84381553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.54214415) q[2];
sx q[2];
rz(-0.52079529) q[2];
sx q[2];
rz(-0.64860541) q[2];
rz(3.0829698) q[3];
sx q[3];
rz(-2.5128745) q[3];
sx q[3];
rz(2.4962943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36515737) q[0];
sx q[0];
rz(-1.9226274) q[0];
sx q[0];
rz(2.6489039) q[0];
rz(1.0538712) q[1];
sx q[1];
rz(-2.5851879) q[1];
sx q[1];
rz(-2.7256706) q[1];
rz(-3.1171534) q[2];
sx q[2];
rz(-1.8069488) q[2];
sx q[2];
rz(-1.9843742) q[2];
rz(-2.8063227) q[3];
sx q[3];
rz(-2.2380968) q[3];
sx q[3];
rz(2.5401881) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
