OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6668532) q[0];
sx q[0];
rz(3.9711877) q[0];
sx q[0];
rz(9.2708099) q[0];
rz(0.83377588) q[1];
sx q[1];
rz(-2.1492465) q[1];
sx q[1];
rz(-0.33831236) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96903893) q[0];
sx q[0];
rz(-1.8882897) q[0];
sx q[0];
rz(-2.88455) q[0];
rz(-pi) q[1];
rz(2.7726735) q[2];
sx q[2];
rz(-0.92637617) q[2];
sx q[2];
rz(1.8298139) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.675128) q[1];
sx q[1];
rz(-0.29106859) q[1];
sx q[1];
rz(2.9558099) q[1];
rz(-pi) q[2];
rz(0.015720856) q[3];
sx q[3];
rz(-1.058488) q[3];
sx q[3];
rz(0.44954625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.14264318) q[2];
sx q[2];
rz(-2.8012186) q[2];
sx q[2];
rz(-1.9677229) q[2];
rz(0.075803444) q[3];
sx q[3];
rz(-1.1444164) q[3];
sx q[3];
rz(0.092806667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3409815) q[0];
sx q[0];
rz(-2.0759463) q[0];
sx q[0];
rz(0.064963438) q[0];
rz(2.5669572) q[1];
sx q[1];
rz(-2.7119633) q[1];
sx q[1];
rz(-1.8992791) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7466465) q[0];
sx q[0];
rz(-1.442369) q[0];
sx q[0];
rz(0.24982474) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93810268) q[2];
sx q[2];
rz(-1.4422851) q[2];
sx q[2];
rz(3.1046257) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4878792) q[1];
sx q[1];
rz(-2.5839845) q[1];
sx q[1];
rz(-2.8370268) q[1];
rz(-pi) q[2];
rz(2.4317125) q[3];
sx q[3];
rz(-1.9544365) q[3];
sx q[3];
rz(0.82304728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3339281) q[2];
sx q[2];
rz(-1.0753205) q[2];
sx q[2];
rz(-2.5578965) q[2];
rz(2.5675473) q[3];
sx q[3];
rz(-2.0160926) q[3];
sx q[3];
rz(3.0103502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42049256) q[0];
sx q[0];
rz(-0.89389602) q[0];
sx q[0];
rz(-2.4131391) q[0];
rz(1.4942253) q[1];
sx q[1];
rz(-2.7431226) q[1];
sx q[1];
rz(2.1247991) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66660488) q[0];
sx q[0];
rz(-0.089086108) q[0];
sx q[0];
rz(-2.7276917) q[0];
rz(-pi) q[1];
x q[1];
rz(2.617308) q[2];
sx q[2];
rz(-2.0031843) q[2];
sx q[2];
rz(2.1081032) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.96117561) q[1];
sx q[1];
rz(-0.75913402) q[1];
sx q[1];
rz(1.3707861) q[1];
rz(-pi) q[2];
rz(-3.0924762) q[3];
sx q[3];
rz(-1.8521063) q[3];
sx q[3];
rz(-2.1360872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3399405) q[2];
sx q[2];
rz(-1.6093971) q[2];
sx q[2];
rz(1.0495079) q[2];
rz(-2.5028051) q[3];
sx q[3];
rz(-0.63126606) q[3];
sx q[3];
rz(-1.1857741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3291572) q[0];
sx q[0];
rz(-1.8742467) q[0];
sx q[0];
rz(1.4720434) q[0];
rz(0.73515785) q[1];
sx q[1];
rz(-2.3627294) q[1];
sx q[1];
rz(-2.8947815) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4287764) q[0];
sx q[0];
rz(-1.719559) q[0];
sx q[0];
rz(-0.87417283) q[0];
rz(-1.3601801) q[2];
sx q[2];
rz(-2.0401376) q[2];
sx q[2];
rz(2.5460555) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2553195) q[1];
sx q[1];
rz(-1.6233994) q[1];
sx q[1];
rz(-2.7707151) q[1];
rz(-pi) q[2];
rz(0.95440063) q[3];
sx q[3];
rz(-1.8064926) q[3];
sx q[3];
rz(-2.1637722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6639158) q[2];
sx q[2];
rz(-2.0229979) q[2];
sx q[2];
rz(1.5412615) q[2];
rz(0.70704308) q[3];
sx q[3];
rz(-1.0975081) q[3];
sx q[3];
rz(2.2533806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8108869) q[0];
sx q[0];
rz(-2.4208477) q[0];
sx q[0];
rz(-1.3274308) q[0];
rz(1.56303) q[1];
sx q[1];
rz(-2.6674318) q[1];
sx q[1];
rz(2.8932103) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43863338) q[0];
sx q[0];
rz(-2.6002433) q[0];
sx q[0];
rz(-0.066141733) q[0];
rz(2.4400473) q[2];
sx q[2];
rz(-0.24214673) q[2];
sx q[2];
rz(2.0602351) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0866962) q[1];
sx q[1];
rz(-0.61453648) q[1];
sx q[1];
rz(-0.30026786) q[1];
rz(-pi) q[2];
rz(-1.8099269) q[3];
sx q[3];
rz(-2.2056747) q[3];
sx q[3];
rz(-2.2374416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7053232) q[2];
sx q[2];
rz(-2.1537809) q[2];
sx q[2];
rz(-2.3441337) q[2];
rz(0.38875368) q[3];
sx q[3];
rz(-0.60384408) q[3];
sx q[3];
rz(-0.50271547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0080863) q[0];
sx q[0];
rz(-0.07645034) q[0];
sx q[0];
rz(-1.7957934) q[0];
rz(-2.0603518) q[1];
sx q[1];
rz(-1.9045647) q[1];
sx q[1];
rz(3.016901) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42449441) q[0];
sx q[0];
rz(-1.6569123) q[0];
sx q[0];
rz(-2.9647102) q[0];
x q[1];
rz(-0.40839809) q[2];
sx q[2];
rz(-0.64986594) q[2];
sx q[2];
rz(-1.0330531) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60285073) q[1];
sx q[1];
rz(-0.83487836) q[1];
sx q[1];
rz(-1.49453) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3966339) q[3];
sx q[3];
rz(-2.8248441) q[3];
sx q[3];
rz(1.6572286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6283915) q[2];
sx q[2];
rz(-1.1867563) q[2];
sx q[2];
rz(-0.61895269) q[2];
rz(-1.0533054) q[3];
sx q[3];
rz(-0.17799938) q[3];
sx q[3];
rz(-2.2119904) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58182794) q[0];
sx q[0];
rz(-1.3904089) q[0];
sx q[0];
rz(-1.0429617) q[0];
rz(-2.6783121) q[1];
sx q[1];
rz(-2.0279341) q[1];
sx q[1];
rz(-2.0708864) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28744222) q[0];
sx q[0];
rz(-1.5713912) q[0];
sx q[0];
rz(-0.48405148) q[0];
rz(-pi) q[1];
rz(-1.7037017) q[2];
sx q[2];
rz(-1.9153567) q[2];
sx q[2];
rz(0.29495707) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0740944) q[1];
sx q[1];
rz(-2.3405582) q[1];
sx q[1];
rz(-1.7098411) q[1];
rz(2.8093852) q[3];
sx q[3];
rz(-1.5899961) q[3];
sx q[3];
rz(1.7393877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.36859194) q[2];
sx q[2];
rz(-0.7395145) q[2];
sx q[2];
rz(0.32361844) q[2];
rz(2.1598024) q[3];
sx q[3];
rz(-2.2798645) q[3];
sx q[3];
rz(-1.0872844) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6535646) q[0];
sx q[0];
rz(-1.7249148) q[0];
sx q[0];
rz(-1.2063684) q[0];
rz(1.9288829) q[1];
sx q[1];
rz(-2.2884463) q[1];
sx q[1];
rz(0.94747296) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065814106) q[0];
sx q[0];
rz(-1.4697207) q[0];
sx q[0];
rz(2.1094735) q[0];
rz(-pi) q[1];
rz(0.1768441) q[2];
sx q[2];
rz(-1.4143922) q[2];
sx q[2];
rz(-0.96207372) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.35093388) q[1];
sx q[1];
rz(-2.5194063) q[1];
sx q[1];
rz(1.15637) q[1];
x q[2];
rz(-1.9318337) q[3];
sx q[3];
rz(-2.3449538) q[3];
sx q[3];
rz(-2.3494997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.56132135) q[2];
sx q[2];
rz(-1.3803955) q[2];
sx q[2];
rz(0.46978152) q[2];
rz(-1.8404768) q[3];
sx q[3];
rz(-1.7062635) q[3];
sx q[3];
rz(-0.27967134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25205055) q[0];
sx q[0];
rz(-0.38953504) q[0];
sx q[0];
rz(1.3289733) q[0];
rz(-0.7912311) q[1];
sx q[1];
rz(-0.33114094) q[1];
sx q[1];
rz(-2.9387617) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.348939) q[0];
sx q[0];
rz(-1.5968423) q[0];
sx q[0];
rz(0.039020122) q[0];
rz(-pi) q[1];
rz(1.3615666) q[2];
sx q[2];
rz(-1.7893357) q[2];
sx q[2];
rz(-0.92313672) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.030414) q[1];
sx q[1];
rz(-1.5215538) q[1];
sx q[1];
rz(-1.4700252) q[1];
rz(-pi) q[2];
rz(-1.4167452) q[3];
sx q[3];
rz(-2.027958) q[3];
sx q[3];
rz(2.9726213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7245076) q[2];
sx q[2];
rz(-2.8736726) q[2];
sx q[2];
rz(0.99651304) q[2];
rz(-2.7881682) q[3];
sx q[3];
rz(-2.3963908) q[3];
sx q[3];
rz(-2.3341808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.080169454) q[0];
sx q[0];
rz(-2.3286979) q[0];
sx q[0];
rz(0.18173519) q[0];
rz(0.043047992) q[1];
sx q[1];
rz(-0.64518607) q[1];
sx q[1];
rz(2.8607686) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1998394) q[0];
sx q[0];
rz(-2.0240677) q[0];
sx q[0];
rz(-1.0928632) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7634723) q[2];
sx q[2];
rz(-1.780605) q[2];
sx q[2];
rz(3.0378621) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9210789) q[1];
sx q[1];
rz(-1.5346569) q[1];
sx q[1];
rz(1.6284579) q[1];
rz(-pi) q[2];
rz(1.9032352) q[3];
sx q[3];
rz(-1.7566163) q[3];
sx q[3];
rz(2.9591054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8250371) q[2];
sx q[2];
rz(-1.2544158) q[2];
sx q[2];
rz(-2.5184856) q[2];
rz(-1.0021707) q[3];
sx q[3];
rz(-1.802417) q[3];
sx q[3];
rz(-0.56308693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4939209) q[0];
sx q[0];
rz(-1.5734084) q[0];
sx q[0];
rz(-1.5403803) q[0];
rz(-2.2676246) q[1];
sx q[1];
rz(-1.0653492) q[1];
sx q[1];
rz(-3.031562) q[1];
rz(0.80550823) q[2];
sx q[2];
rz(-2.0279573) q[2];
sx q[2];
rz(1.3062994) q[2];
rz(-0.35691805) q[3];
sx q[3];
rz(-1.1834984) q[3];
sx q[3];
rz(-0.055565861) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
