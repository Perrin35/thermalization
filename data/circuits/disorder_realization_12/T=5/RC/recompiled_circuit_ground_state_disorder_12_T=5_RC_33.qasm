OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0491068) q[0];
sx q[0];
rz(-1.8422814) q[0];
sx q[0];
rz(0.045825034) q[0];
rz(1.2996281) q[1];
sx q[1];
rz(-1.8027432) q[1];
sx q[1];
rz(-1.8884698) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0631615) q[0];
sx q[0];
rz(-2.6303311) q[0];
sx q[0];
rz(-1.0529165) q[0];
x q[1];
rz(-2.9304977) q[2];
sx q[2];
rz(-2.935754) q[2];
sx q[2];
rz(1.7072276) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.10132566) q[1];
sx q[1];
rz(-0.76124746) q[1];
sx q[1];
rz(-2.9684307) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6982444) q[3];
sx q[3];
rz(-1.5747254) q[3];
sx q[3];
rz(1.9779825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8200298) q[2];
sx q[2];
rz(-1.7367312) q[2];
sx q[2];
rz(-0.31537867) q[2];
rz(-3.0553715) q[3];
sx q[3];
rz(-0.092970522) q[3];
sx q[3];
rz(2.2975547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.4013937) q[0];
sx q[0];
rz(-1.711015) q[0];
sx q[0];
rz(-2.8061818) q[0];
rz(2.7566578) q[1];
sx q[1];
rz(-2.3454869) q[1];
sx q[1];
rz(1.8523432) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9612564) q[0];
sx q[0];
rz(-1.633005) q[0];
sx q[0];
rz(2.7520685) q[0];
rz(-pi) q[1];
rz(2.2122967) q[2];
sx q[2];
rz(-1.4946862) q[2];
sx q[2];
rz(-2.2141505) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1114167) q[1];
sx q[1];
rz(-1.1892288) q[1];
sx q[1];
rz(1.8757344) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13522526) q[3];
sx q[3];
rz(-1.3298508) q[3];
sx q[3];
rz(-1.2775482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5424767) q[2];
sx q[2];
rz(-1.3834388) q[2];
sx q[2];
rz(1.6790338) q[2];
rz(-2.2595432) q[3];
sx q[3];
rz(-2.4817011) q[3];
sx q[3];
rz(1.5006458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3682692) q[0];
sx q[0];
rz(-2.6550846) q[0];
sx q[0];
rz(2.4988556) q[0];
rz(0.87359387) q[1];
sx q[1];
rz(-1.3106376) q[1];
sx q[1];
rz(2.1547623) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5499134) q[0];
sx q[0];
rz(-2.4351623) q[0];
sx q[0];
rz(1.3384992) q[0];
x q[1];
rz(3.0796649) q[2];
sx q[2];
rz(-1.0871097) q[2];
sx q[2];
rz(2.7360423) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3915186) q[1];
sx q[1];
rz(-1.0315391) q[1];
sx q[1];
rz(-2.2184994) q[1];
x q[2];
rz(1.945472) q[3];
sx q[3];
rz(-1.5244532) q[3];
sx q[3];
rz(-2.1425193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.70283908) q[2];
sx q[2];
rz(-3.1162016) q[2];
sx q[2];
rz(-1.0566443) q[2];
rz(1.3422525) q[3];
sx q[3];
rz(-1.6868351) q[3];
sx q[3];
rz(2.2630283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1037647) q[0];
sx q[0];
rz(-0.28452888) q[0];
sx q[0];
rz(1.1430662) q[0];
rz(2.4567538) q[1];
sx q[1];
rz(-1.7999444) q[1];
sx q[1];
rz(-0.24982223) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9651325) q[0];
sx q[0];
rz(-1.1145381) q[0];
sx q[0];
rz(0.81911032) q[0];
rz(-pi) q[1];
rz(-1.5748137) q[2];
sx q[2];
rz(-2.145014) q[2];
sx q[2];
rz(1.9030273) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6982242) q[1];
sx q[1];
rz(-2.0495546) q[1];
sx q[1];
rz(0.67826346) q[1];
x q[2];
rz(1.5771082) q[3];
sx q[3];
rz(-2.8996116) q[3];
sx q[3];
rz(-2.210287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7217241) q[2];
sx q[2];
rz(-1.3519752) q[2];
sx q[2];
rz(-1.295759) q[2];
rz(-1.7660247) q[3];
sx q[3];
rz(-1.4294759) q[3];
sx q[3];
rz(-0.099055722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7527723) q[0];
sx q[0];
rz(-1.3714014) q[0];
sx q[0];
rz(2.2950628) q[0];
rz(1.1835774) q[1];
sx q[1];
rz(-1.1776244) q[1];
sx q[1];
rz(-1.902098) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0060109) q[0];
sx q[0];
rz(-0.2811037) q[0];
sx q[0];
rz(1.1233888) q[0];
rz(0.34244142) q[2];
sx q[2];
rz(-1.7424118) q[2];
sx q[2];
rz(-1.1521074) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.848865) q[1];
sx q[1];
rz(-1.5588817) q[1];
sx q[1];
rz(1.5148276) q[1];
rz(0.97364135) q[3];
sx q[3];
rz(-1.6203364) q[3];
sx q[3];
rz(1.3508391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6401297) q[2];
sx q[2];
rz(-1.632246) q[2];
sx q[2];
rz(-0.35422361) q[2];
rz(0.7192449) q[3];
sx q[3];
rz(-0.57643276) q[3];
sx q[3];
rz(-0.80938068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3767913) q[0];
sx q[0];
rz(-0.23623315) q[0];
sx q[0];
rz(-2.7948622) q[0];
rz(2.4193343) q[1];
sx q[1];
rz(-1.6112695) q[1];
sx q[1];
rz(0.17682704) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.090031) q[0];
sx q[0];
rz(-1.430184) q[0];
sx q[0];
rz(-2.0363801) q[0];
x q[1];
rz(2.0697303) q[2];
sx q[2];
rz(-1.4889354) q[2];
sx q[2];
rz(2.6112542) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2894779) q[1];
sx q[1];
rz(-1.3657709) q[1];
sx q[1];
rz(1.6774898) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1336055) q[3];
sx q[3];
rz(-1.4848108) q[3];
sx q[3];
rz(-0.028349625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.56883812) q[2];
sx q[2];
rz(-0.73358959) q[2];
sx q[2];
rz(-1.0432165) q[2];
rz(-0.15626945) q[3];
sx q[3];
rz(-1.0633435) q[3];
sx q[3];
rz(0.094303057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8909376) q[0];
sx q[0];
rz(-1.0303048) q[0];
sx q[0];
rz(1.2283196) q[0];
rz(2.7078775) q[1];
sx q[1];
rz(-2.3408196) q[1];
sx q[1];
rz(2.0965651) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7997847) q[0];
sx q[0];
rz(-1.5130454) q[0];
sx q[0];
rz(1.7736797) q[0];
rz(-pi) q[1];
rz(-0.82257338) q[2];
sx q[2];
rz(-1.2116836) q[2];
sx q[2];
rz(2.1445779) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1365388) q[1];
sx q[1];
rz(-0.95816678) q[1];
sx q[1];
rz(3.0023252) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0923474) q[3];
sx q[3];
rz(-2.2766067) q[3];
sx q[3];
rz(1.7569831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9016483) q[2];
sx q[2];
rz(-1.0239235) q[2];
sx q[2];
rz(0.18590064) q[2];
rz(-1.2402041) q[3];
sx q[3];
rz(-1.5408206) q[3];
sx q[3];
rz(2.127229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0094725322) q[0];
sx q[0];
rz(-1.8931696) q[0];
sx q[0];
rz(-0.15952071) q[0];
rz(2.3347143) q[1];
sx q[1];
rz(-1.2346377) q[1];
sx q[1];
rz(0.0085011403) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5497564) q[0];
sx q[0];
rz(-1.4167794) q[0];
sx q[0];
rz(-1.9817673) q[0];
rz(-pi) q[1];
x q[1];
rz(0.10430704) q[2];
sx q[2];
rz(-2.3416714) q[2];
sx q[2];
rz(-1.4664354) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.22903012) q[1];
sx q[1];
rz(-2.8797315) q[1];
sx q[1];
rz(-0.53066855) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2448266) q[3];
sx q[3];
rz(-2.634832) q[3];
sx q[3];
rz(2.4017577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.27935091) q[2];
sx q[2];
rz(-1.778435) q[2];
sx q[2];
rz(2.9478493) q[2];
rz(1.6019609) q[3];
sx q[3];
rz(-1.1774457) q[3];
sx q[3];
rz(2.0120473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2368471) q[0];
sx q[0];
rz(-1.756825) q[0];
sx q[0];
rz(-2.4054476) q[0];
rz(-0.88788095) q[1];
sx q[1];
rz(-0.45279756) q[1];
sx q[1];
rz(0.050051659) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0607292) q[0];
sx q[0];
rz(-1.9115149) q[0];
sx q[0];
rz(1.7572035) q[0];
rz(-pi) q[1];
rz(2.8269269) q[2];
sx q[2];
rz(-1.6648714) q[2];
sx q[2];
rz(-1.5621746) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1117275) q[1];
sx q[1];
rz(-1.7741412) q[1];
sx q[1];
rz(-0.3158658) q[1];
x q[2];
rz(-2.5727651) q[3];
sx q[3];
rz(-2.1095368) q[3];
sx q[3];
rz(-2.8728268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8686409) q[2];
sx q[2];
rz(-1.7014039) q[2];
sx q[2];
rz(3.1325373) q[2];
rz(-2.791259) q[3];
sx q[3];
rz(-2.83559) q[3];
sx q[3];
rz(2.8370324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41651273) q[0];
sx q[0];
rz(-1.059499) q[0];
sx q[0];
rz(2.1685725) q[0];
rz(-1.5618207) q[1];
sx q[1];
rz(-0.90619722) q[1];
sx q[1];
rz(-0.80950338) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0334629) q[0];
sx q[0];
rz(-1.9282883) q[0];
sx q[0];
rz(-2.3603159) q[0];
x q[1];
rz(-1.3415178) q[2];
sx q[2];
rz(-1.3827168) q[2];
sx q[2];
rz(-0.089872472) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.36724801) q[1];
sx q[1];
rz(-1.4092236) q[1];
sx q[1];
rz(2.0366686) q[1];
rz(-pi) q[2];
rz(2.6736835) q[3];
sx q[3];
rz(-2.10026) q[3];
sx q[3];
rz(2.5264945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1805264) q[2];
sx q[2];
rz(-2.3193391) q[2];
sx q[2];
rz(-0.52214885) q[2];
rz(-0.3565878) q[3];
sx q[3];
rz(-0.57727376) q[3];
sx q[3];
rz(3.0780415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6575573) q[0];
sx q[0];
rz(-2.7900896) q[0];
sx q[0];
rz(-1.8650613) q[0];
rz(2.7784078) q[1];
sx q[1];
rz(-2.6138432) q[1];
sx q[1];
rz(1.6539727) q[1];
rz(-0.54375081) q[2];
sx q[2];
rz(-0.85493543) q[2];
sx q[2];
rz(-2.2731352) q[2];
rz(1.894886) q[3];
sx q[3];
rz(-1.0151328) q[3];
sx q[3];
rz(-0.48965164) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
