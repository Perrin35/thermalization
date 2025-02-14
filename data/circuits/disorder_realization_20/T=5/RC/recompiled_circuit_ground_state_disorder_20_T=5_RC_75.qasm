OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5372758) q[0];
sx q[0];
rz(-0.24157) q[0];
sx q[0];
rz(0.33302745) q[0];
rz(-1.1355407) q[1];
sx q[1];
rz(-2.314664) q[1];
sx q[1];
rz(0.64396042) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49227958) q[0];
sx q[0];
rz(-0.67848533) q[0];
sx q[0];
rz(-1.1301103) q[0];
x q[1];
rz(1.4139373) q[2];
sx q[2];
rz(-1.867395) q[2];
sx q[2];
rz(3.0576884) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.8713069) q[1];
sx q[1];
rz(-1.3340923) q[1];
sx q[1];
rz(-2.7530297) q[1];
x q[2];
rz(2.7936739) q[3];
sx q[3];
rz(-1.9609465) q[3];
sx q[3];
rz(1.8371234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.48646271) q[2];
sx q[2];
rz(-2.6772406) q[2];
sx q[2];
rz(0.74074024) q[2];
rz(1.5517392) q[3];
sx q[3];
rz(-0.71151763) q[3];
sx q[3];
rz(0.5927425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7084259) q[0];
sx q[0];
rz(-1.1335224) q[0];
sx q[0];
rz(0.11696996) q[0];
rz(-2.6372657) q[1];
sx q[1];
rz(-1.3386936) q[1];
sx q[1];
rz(-2.8574944) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4000716) q[0];
sx q[0];
rz(-1.0062381) q[0];
sx q[0];
rz(0.098640504) q[0];
x q[1];
rz(-2.1524736) q[2];
sx q[2];
rz(-2.2788413) q[2];
sx q[2];
rz(3.0281554) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5631905) q[1];
sx q[1];
rz(-0.12688533) q[1];
sx q[1];
rz(1.2940426) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25155622) q[3];
sx q[3];
rz(-1.8402035) q[3];
sx q[3];
rz(-0.081307383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32610193) q[2];
sx q[2];
rz(-2.2559866) q[2];
sx q[2];
rz(0.76078129) q[2];
rz(2.271999) q[3];
sx q[3];
rz(-1.875501) q[3];
sx q[3];
rz(2.6884955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78075439) q[0];
sx q[0];
rz(-2.0212845) q[0];
sx q[0];
rz(-0.25217062) q[0];
rz(1.699327) q[1];
sx q[1];
rz(-0.95528722) q[1];
sx q[1];
rz(2.7853277) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8460444) q[0];
sx q[0];
rz(-1.6303807) q[0];
sx q[0];
rz(-0.0061959717) q[0];
rz(-0.54889955) q[2];
sx q[2];
rz(-1.2147012) q[2];
sx q[2];
rz(0.20129542) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1504476) q[1];
sx q[1];
rz(-2.4887062) q[1];
sx q[1];
rz(-1.3351403) q[1];
rz(3.1396418) q[3];
sx q[3];
rz(-1.3144799) q[3];
sx q[3];
rz(-1.0904097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1189271) q[2];
sx q[2];
rz(-1.7513195) q[2];
sx q[2];
rz(1.6290132) q[2];
rz(3.1318943) q[3];
sx q[3];
rz(-1.2303979) q[3];
sx q[3];
rz(2.5201918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6464722) q[0];
sx q[0];
rz(-2.5992114) q[0];
sx q[0];
rz(-0.25076184) q[0];
rz(-1.3465025) q[1];
sx q[1];
rz(-0.52918068) q[1];
sx q[1];
rz(-1.4432602) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8478717) q[0];
sx q[0];
rz(-1.3547055) q[0];
sx q[0];
rz(-0.8190852) q[0];
rz(-pi) q[1];
rz(0.94130959) q[2];
sx q[2];
rz(-2.513859) q[2];
sx q[2];
rz(2.8986487) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5233744) q[1];
sx q[1];
rz(-0.95444767) q[1];
sx q[1];
rz(-0.90202721) q[1];
rz(-pi) q[2];
rz(-0.26224995) q[3];
sx q[3];
rz(-1.3104386) q[3];
sx q[3];
rz(-2.8182056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5859588) q[2];
sx q[2];
rz(-1.1052174) q[2];
sx q[2];
rz(0.24331681) q[2];
rz(-2.5569052) q[3];
sx q[3];
rz(-0.57716113) q[3];
sx q[3];
rz(0.028133597) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1481767) q[0];
sx q[0];
rz(-0.36574829) q[0];
sx q[0];
rz(-2.962033) q[0];
rz(1.9163632) q[1];
sx q[1];
rz(-1.7370677) q[1];
sx q[1];
rz(-0.45825759) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4438547) q[0];
sx q[0];
rz(-0.30618822) q[0];
sx q[0];
rz(-0.26060391) q[0];
rz(1.304428) q[2];
sx q[2];
rz(-2.2812216) q[2];
sx q[2];
rz(-1.923234) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0441664) q[1];
sx q[1];
rz(-1.4961745) q[1];
sx q[1];
rz(1.4358372) q[1];
rz(2.2716801) q[3];
sx q[3];
rz(-1.2906089) q[3];
sx q[3];
rz(-0.82893574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3852343) q[2];
sx q[2];
rz(-2.7172654) q[2];
sx q[2];
rz(-2.2198086) q[2];
rz(-1.8900169) q[3];
sx q[3];
rz(-1.7332964) q[3];
sx q[3];
rz(0.061669953) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.09403041) q[0];
sx q[0];
rz(-2.229409) q[0];
sx q[0];
rz(-2.9929274) q[0];
rz(2.8246763) q[1];
sx q[1];
rz(-0.47873679) q[1];
sx q[1];
rz(2.5247578) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6233404) q[0];
sx q[0];
rz(-3.0419311) q[0];
sx q[0];
rz(-0.90604337) q[0];
rz(-pi) q[1];
rz(2.7814034) q[2];
sx q[2];
rz(-1.9860886) q[2];
sx q[2];
rz(-0.54999888) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.65766636) q[1];
sx q[1];
rz(-1.0063419) q[1];
sx q[1];
rz(0.29740833) q[1];
x q[2];
rz(-2.0944164) q[3];
sx q[3];
rz(-1.8111808) q[3];
sx q[3];
rz(-2.2464858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.8646912) q[2];
sx q[2];
rz(-1.7128877) q[2];
sx q[2];
rz(0.0083262715) q[2];
rz(0.62638038) q[3];
sx q[3];
rz(-0.71599394) q[3];
sx q[3];
rz(3.1410419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3061227) q[0];
sx q[0];
rz(-1.1518814) q[0];
sx q[0];
rz(-0.48208153) q[0];
rz(-3.112402) q[1];
sx q[1];
rz(-1.8455285) q[1];
sx q[1];
rz(-0.76404244) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69224834) q[0];
sx q[0];
rz(-1.5034165) q[0];
sx q[0];
rz(1.199556) q[0];
x q[1];
rz(-0.025310658) q[2];
sx q[2];
rz(-0.9111852) q[2];
sx q[2];
rz(0.74756223) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.75023821) q[1];
sx q[1];
rz(-1.7403894) q[1];
sx q[1];
rz(2.1743348) q[1];
rz(-1.5519616) q[3];
sx q[3];
rz(-0.40098396) q[3];
sx q[3];
rz(-2.8249521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.46334106) q[2];
sx q[2];
rz(-1.3193139) q[2];
sx q[2];
rz(2.6574262) q[2];
rz(-0.68228996) q[3];
sx q[3];
rz(-0.65410084) q[3];
sx q[3];
rz(-3.0068523) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.057137) q[0];
sx q[0];
rz(-2.3637922) q[0];
sx q[0];
rz(-1.8498259) q[0];
rz(2.6765587) q[1];
sx q[1];
rz(-2.6227622) q[1];
sx q[1];
rz(0.30716392) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43873337) q[0];
sx q[0];
rz(-2.4552058) q[0];
sx q[0];
rz(0.86875654) q[0];
rz(-2.1910153) q[2];
sx q[2];
rz(-0.75846106) q[2];
sx q[2];
rz(1.72067) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0548361) q[1];
sx q[1];
rz(-2.1385178) q[1];
sx q[1];
rz(-0.18594976) q[1];
rz(-pi) q[2];
x q[2];
rz(0.028178111) q[3];
sx q[3];
rz(-2.5163076) q[3];
sx q[3];
rz(-2.6373088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.18053599) q[2];
sx q[2];
rz(-2.7013216) q[2];
sx q[2];
rz(-2.4214936) q[2];
rz(-2.2911206) q[3];
sx q[3];
rz(-1.1668147) q[3];
sx q[3];
rz(-1.4115964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20486031) q[0];
sx q[0];
rz(-2.9731049) q[0];
sx q[0];
rz(-2.4334461) q[0];
rz(1.6234966) q[1];
sx q[1];
rz(-2.203439) q[1];
sx q[1];
rz(0.35619563) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4771381) q[0];
sx q[0];
rz(-1.150048) q[0];
sx q[0];
rz(-0.53945213) q[0];
rz(-0.61788606) q[2];
sx q[2];
rz(-1.5466585) q[2];
sx q[2];
rz(1.9113505) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1005786) q[1];
sx q[1];
rz(-1.5407526) q[1];
sx q[1];
rz(0.64893367) q[1];
rz(-pi) q[2];
rz(-2.735504) q[3];
sx q[3];
rz(-1.1392987) q[3];
sx q[3];
rz(-0.0270947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0492964) q[2];
sx q[2];
rz(-2.3342817) q[2];
sx q[2];
rz(2.2558007) q[2];
rz(-2.2057335) q[3];
sx q[3];
rz(-1.1805308) q[3];
sx q[3];
rz(0.22127557) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43779272) q[0];
sx q[0];
rz(-0.047310345) q[0];
sx q[0];
rz(-1.6920775) q[0];
rz(-1.0225147) q[1];
sx q[1];
rz(-2.6578564) q[1];
sx q[1];
rz(-1.449301) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.463844) q[0];
sx q[0];
rz(-2.7263256) q[0];
sx q[0];
rz(-1.6941316) q[0];
x q[1];
rz(-2.3751276) q[2];
sx q[2];
rz(-2.577707) q[2];
sx q[2];
rz(1.6414798) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3811581) q[1];
sx q[1];
rz(-1.0236003) q[1];
sx q[1];
rz(-2.8924499) q[1];
rz(-2.0725771) q[3];
sx q[3];
rz(-0.3007362) q[3];
sx q[3];
rz(-1.5141443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.25541043) q[2];
sx q[2];
rz(-1.1899199) q[2];
sx q[2];
rz(-2.4667242) q[2];
rz(-0.97149649) q[3];
sx q[3];
rz(-0.70236218) q[3];
sx q[3];
rz(-1.6500047) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7991199) q[0];
sx q[0];
rz(-1.5457038) q[0];
sx q[0];
rz(-0.85734838) q[0];
rz(2.9091861) q[1];
sx q[1];
rz(-1.0127761) q[1];
sx q[1];
rz(1.3565328) q[1];
rz(1.0398374) q[2];
sx q[2];
rz(-1.2086443) q[2];
sx q[2];
rz(2.3966387) q[2];
rz(-0.078217004) q[3];
sx q[3];
rz(-0.84368869) q[3];
sx q[3];
rz(0.61144184) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
