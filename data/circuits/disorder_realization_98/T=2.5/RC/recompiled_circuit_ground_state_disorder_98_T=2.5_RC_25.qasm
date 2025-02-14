OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.15220517) q[0];
sx q[0];
rz(-0.21259354) q[0];
sx q[0];
rz(0.73417443) q[0];
rz(-1.9250159) q[1];
sx q[1];
rz(-0.34595481) q[1];
sx q[1];
rz(3.1168361) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5261794) q[0];
sx q[0];
rz(-0.075881474) q[0];
sx q[0];
rz(-0.60657255) q[0];
rz(-1.8046298) q[2];
sx q[2];
rz(-1.6811835) q[2];
sx q[2];
rz(1.259492) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1851968) q[1];
sx q[1];
rz(-1.6231771) q[1];
sx q[1];
rz(2.8422794) q[1];
rz(-1.7832028) q[3];
sx q[3];
rz(-1.6254566) q[3];
sx q[3];
rz(2.2509991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.88432246) q[2];
sx q[2];
rz(-1.4929644) q[2];
sx q[2];
rz(-2.8327668) q[2];
rz(2.3158) q[3];
sx q[3];
rz(-2.1335996) q[3];
sx q[3];
rz(-2.3776313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9708213) q[0];
sx q[0];
rz(-1.229137) q[0];
sx q[0];
rz(2.1139076) q[0];
rz(-0.81614256) q[1];
sx q[1];
rz(-2.0232537) q[1];
sx q[1];
rz(2.3291086) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6861758) q[0];
sx q[0];
rz(-2.7551965) q[0];
sx q[0];
rz(-2.0040345) q[0];
rz(1.585) q[2];
sx q[2];
rz(-1.1974632) q[2];
sx q[2];
rz(0.65711428) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5066487) q[1];
sx q[1];
rz(-2.6104984) q[1];
sx q[1];
rz(-1.2960766) q[1];
rz(-pi) q[2];
rz(-2.3458293) q[3];
sx q[3];
rz(-0.47945346) q[3];
sx q[3];
rz(-0.68414068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.57427788) q[2];
sx q[2];
rz(-1.652176) q[2];
sx q[2];
rz(-0.91892773) q[2];
rz(-0.53141665) q[3];
sx q[3];
rz(-1.0454949) q[3];
sx q[3];
rz(-0.04118583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4384005) q[0];
sx q[0];
rz(-2.0719318) q[0];
sx q[0];
rz(1.137314) q[0];
rz(1.7510341) q[1];
sx q[1];
rz(-0.5568234) q[1];
sx q[1];
rz(-2.4086319) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82405381) q[0];
sx q[0];
rz(-1.8630872) q[0];
sx q[0];
rz(-0.43265105) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5733092) q[2];
sx q[2];
rz(-2.097762) q[2];
sx q[2];
rz(0.87430853) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3324686) q[1];
sx q[1];
rz(-0.66624236) q[1];
sx q[1];
rz(-2.3205643) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0065782733) q[3];
sx q[3];
rz(-2.0195144) q[3];
sx q[3];
rz(-2.6980163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1375894) q[2];
sx q[2];
rz(-1.4121476) q[2];
sx q[2];
rz(-2.1027193) q[2];
rz(-1.002257) q[3];
sx q[3];
rz(-1.301731) q[3];
sx q[3];
rz(0.29016289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20659474) q[0];
sx q[0];
rz(-2.0552141) q[0];
sx q[0];
rz(-0.90718734) q[0];
rz(0.15779237) q[1];
sx q[1];
rz(-2.1267499) q[1];
sx q[1];
rz(-1.8603604) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7665721) q[0];
sx q[0];
rz(-2.2849884) q[0];
sx q[0];
rz(-0.71582224) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1195391) q[2];
sx q[2];
rz(-0.76821487) q[2];
sx q[2];
rz(0.31203416) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.43739) q[1];
sx q[1];
rz(-1.1466195) q[1];
sx q[1];
rz(-0.68628879) q[1];
x q[2];
rz(-2.2415555) q[3];
sx q[3];
rz(-2.0931912) q[3];
sx q[3];
rz(2.1417422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24718757) q[2];
sx q[2];
rz(-2.2053568) q[2];
sx q[2];
rz(1.8518764) q[2];
rz(1.3442518) q[3];
sx q[3];
rz(-1.684609) q[3];
sx q[3];
rz(-2.4414506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6808788) q[0];
sx q[0];
rz(-2.6432156) q[0];
sx q[0];
rz(1.0182678) q[0];
rz(-1.1030039) q[1];
sx q[1];
rz(-0.7119199) q[1];
sx q[1];
rz(-1.3404554) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8333574) q[0];
sx q[0];
rz(-1.0923704) q[0];
sx q[0];
rz(-1.1638948) q[0];
rz(-0.37968882) q[2];
sx q[2];
rz(-1.4799397) q[2];
sx q[2];
rz(2.1374709) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3892438) q[1];
sx q[1];
rz(-1.443814) q[1];
sx q[1];
rz(-1.4402706) q[1];
rz(-1.0304673) q[3];
sx q[3];
rz(-1.6202462) q[3];
sx q[3];
rz(-1.7211308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5420142) q[2];
sx q[2];
rz(-0.88853637) q[2];
sx q[2];
rz(-2.1395903) q[2];
rz(0.69139785) q[3];
sx q[3];
rz(-1.1804429) q[3];
sx q[3];
rz(-1.5207759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6028676) q[0];
sx q[0];
rz(-1.069101) q[0];
sx q[0];
rz(0.62527239) q[0];
rz(-1.9484776) q[1];
sx q[1];
rz(-0.68980491) q[1];
sx q[1];
rz(2.5742721) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1330018) q[0];
sx q[0];
rz(-2.3390305) q[0];
sx q[0];
rz(0.24548291) q[0];
rz(1.5789323) q[2];
sx q[2];
rz(-2.3971791) q[2];
sx q[2];
rz(-1.2071963) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7116136) q[1];
sx q[1];
rz(-0.95250722) q[1];
sx q[1];
rz(-3.1097163) q[1];
rz(2.3472957) q[3];
sx q[3];
rz(-2.8927589) q[3];
sx q[3];
rz(2.2373696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.83711964) q[2];
sx q[2];
rz(-1.2834872) q[2];
sx q[2];
rz(1.5817969) q[2];
rz(0.61257735) q[3];
sx q[3];
rz(-1.160683) q[3];
sx q[3];
rz(-0.15577236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0953858) q[0];
sx q[0];
rz(-1.950773) q[0];
sx q[0];
rz(2.0147391) q[0];
rz(1.2696179) q[1];
sx q[1];
rz(-1.9536628) q[1];
sx q[1];
rz(-2.6205305) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57018897) q[0];
sx q[0];
rz(-0.91556433) q[0];
sx q[0];
rz(2.6960424) q[0];
rz(-pi) q[1];
rz(-0.032835788) q[2];
sx q[2];
rz(-1.2419789) q[2];
sx q[2];
rz(-1.1177899) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4849629) q[1];
sx q[1];
rz(-0.99492517) q[1];
sx q[1];
rz(0.63159512) q[1];
rz(-pi) q[2];
rz(0.43989681) q[3];
sx q[3];
rz(-2.5554113) q[3];
sx q[3];
rz(-0.90914721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1178939) q[2];
sx q[2];
rz(-1.3641027) q[2];
sx q[2];
rz(1.2269616) q[2];
rz(1.4620694) q[3];
sx q[3];
rz(-1.5173802) q[3];
sx q[3];
rz(-0.44155651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11314497) q[0];
sx q[0];
rz(-2.1118836) q[0];
sx q[0];
rz(-0.5471158) q[0];
rz(3.0534577) q[1];
sx q[1];
rz(-1.3476177) q[1];
sx q[1];
rz(-0.42253447) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28332253) q[0];
sx q[0];
rz(-1.8517324) q[0];
sx q[0];
rz(3.1188117) q[0];
rz(-1.2014548) q[2];
sx q[2];
rz(-1.0779194) q[2];
sx q[2];
rz(0.77134672) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9692058) q[1];
sx q[1];
rz(-1.5017461) q[1];
sx q[1];
rz(0.056180908) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2219031) q[3];
sx q[3];
rz(-2.9125179) q[3];
sx q[3];
rz(-3.0034415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2253458) q[2];
sx q[2];
rz(-1.5044745) q[2];
sx q[2];
rz(2.8543191) q[2];
rz(0.54287994) q[3];
sx q[3];
rz(-2.3714239) q[3];
sx q[3];
rz(3.0834901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3995689) q[0];
sx q[0];
rz(-2.4574807) q[0];
sx q[0];
rz(2.3642819) q[0];
rz(-0.9264535) q[1];
sx q[1];
rz(-1.9994206) q[1];
sx q[1];
rz(-1.3949589) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3789916) q[0];
sx q[0];
rz(-1.9585573) q[0];
sx q[0];
rz(-0.34404018) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5424273) q[2];
sx q[2];
rz(-1.3627421) q[2];
sx q[2];
rz(-0.64917246) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3923126) q[1];
sx q[1];
rz(-0.5764851) q[1];
sx q[1];
rz(2.2063401) q[1];
rz(-0.80628245) q[3];
sx q[3];
rz(-0.7976992) q[3];
sx q[3];
rz(-0.21569852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.42561439) q[2];
sx q[2];
rz(-1.7691111) q[2];
sx q[2];
rz(-2.051029) q[2];
rz(-2.1770832) q[3];
sx q[3];
rz(-1.0114074) q[3];
sx q[3];
rz(-1.9264268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5665117) q[0];
sx q[0];
rz(-2.8031741) q[0];
sx q[0];
rz(1.6315208) q[0];
rz(0.034320023) q[1];
sx q[1];
rz(-1.3395373) q[1];
sx q[1];
rz(2.7095749) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1572947) q[0];
sx q[0];
rz(-1.0327067) q[0];
sx q[0];
rz(0.068416083) q[0];
rz(-pi) q[1];
rz(-3.0598214) q[2];
sx q[2];
rz(-1.9207947) q[2];
sx q[2];
rz(3.0584665) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2434531) q[1];
sx q[1];
rz(-1.6695078) q[1];
sx q[1];
rz(2.9961186) q[1];
rz(-1.9625147) q[3];
sx q[3];
rz(-0.61810247) q[3];
sx q[3];
rz(2.6276692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4711275) q[2];
sx q[2];
rz(-1.9894783) q[2];
sx q[2];
rz(2.067789) q[2];
rz(-0.18887575) q[3];
sx q[3];
rz(-2.5329068) q[3];
sx q[3];
rz(2.7591738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043924532) q[0];
sx q[0];
rz(-2.8207939) q[0];
sx q[0];
rz(2.4378142) q[0];
rz(1.7882998) q[1];
sx q[1];
rz(-2.4663993) q[1];
sx q[1];
rz(-0.83723062) q[1];
rz(-2.4953669) q[2];
sx q[2];
rz(-2.2868509) q[2];
sx q[2];
rz(1.7522191) q[2];
rz(0.45936361) q[3];
sx q[3];
rz(-1.50926) q[3];
sx q[3];
rz(2.0987233) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
