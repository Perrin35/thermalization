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
rz(-1.902154) q[0];
sx q[0];
rz(-1.3286989) q[0];
sx q[0];
rz(-2.9297096) q[0];
rz(1.7243241) q[1];
sx q[1];
rz(-0.53242004) q[1];
sx q[1];
rz(2.7636757) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48403063) q[0];
sx q[0];
rz(-0.0095417984) q[0];
sx q[0];
rz(1.6173167) q[0];
x q[1];
rz(1.2075519) q[2];
sx q[2];
rz(-0.9245199) q[2];
sx q[2];
rz(-0.56528795) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4204882) q[1];
sx q[1];
rz(-1.2961565) q[1];
sx q[1];
rz(2.4824449) q[1];
x q[2];
rz(-2.3232949) q[3];
sx q[3];
rz(-2.1967271) q[3];
sx q[3];
rz(-1.4656386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8523031) q[2];
sx q[2];
rz(-1.0785582) q[2];
sx q[2];
rz(-3.0618073) q[2];
rz(0.96528178) q[3];
sx q[3];
rz(-1.691317) q[3];
sx q[3];
rz(0.7307581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47736436) q[0];
sx q[0];
rz(-1.0034765) q[0];
sx q[0];
rz(0.088951237) q[0];
rz(1.3720007) q[1];
sx q[1];
rz(-1.7141432) q[1];
sx q[1];
rz(-2.037183) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4084642) q[0];
sx q[0];
rz(-2.1845572) q[0];
sx q[0];
rz(-1.852714) q[0];
x q[1];
rz(2.0599819) q[2];
sx q[2];
rz(-2.2913165) q[2];
sx q[2];
rz(-0.94948506) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4648393) q[1];
sx q[1];
rz(-1.1285121) q[1];
sx q[1];
rz(0.16525903) q[1];
rz(-pi) q[2];
rz(-0.16544754) q[3];
sx q[3];
rz(-0.54537994) q[3];
sx q[3];
rz(-2.489413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8770807) q[2];
sx q[2];
rz(-2.0161714) q[2];
sx q[2];
rz(-1.4667) q[2];
rz(3.1125715) q[3];
sx q[3];
rz(-1.0163739) q[3];
sx q[3];
rz(-2.4513054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90917176) q[0];
sx q[0];
rz(-0.066815289) q[0];
sx q[0];
rz(0.27798852) q[0];
rz(1.4311283) q[1];
sx q[1];
rz(-2.2093096) q[1];
sx q[1];
rz(0.40036449) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5059197) q[0];
sx q[0];
rz(-0.95209661) q[0];
sx q[0];
rz(-1.0158422) q[0];
rz(-pi) q[1];
rz(1.3157) q[2];
sx q[2];
rz(-0.9005024) q[2];
sx q[2];
rz(-2.3885661) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6975721) q[1];
sx q[1];
rz(-0.50067798) q[1];
sx q[1];
rz(0.33659192) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2648029) q[3];
sx q[3];
rz(-1.8974432) q[3];
sx q[3];
rz(1.5375002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6843188) q[2];
sx q[2];
rz(-1.532734) q[2];
sx q[2];
rz(-2.686783) q[2];
rz(1.8999892) q[3];
sx q[3];
rz(-2.1865215) q[3];
sx q[3];
rz(-2.4212867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.95551816) q[0];
sx q[0];
rz(-1.3294514) q[0];
sx q[0];
rz(0.00051001471) q[0];
rz(-0.60091248) q[1];
sx q[1];
rz(-2.3020703) q[1];
sx q[1];
rz(-0.14437637) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2287284) q[0];
sx q[0];
rz(-1.0236386) q[0];
sx q[0];
rz(0.24258258) q[0];
x q[1];
rz(-1.4207441) q[2];
sx q[2];
rz(-1.320082) q[2];
sx q[2];
rz(1.5992407) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2541009) q[1];
sx q[1];
rz(-0.93587592) q[1];
sx q[1];
rz(-1.2870406) q[1];
x q[2];
rz(2.9668429) q[3];
sx q[3];
rz(-0.55395444) q[3];
sx q[3];
rz(-0.21077158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.942975) q[2];
sx q[2];
rz(-1.696442) q[2];
sx q[2];
rz(-2.1790738) q[2];
rz(1.6019542) q[3];
sx q[3];
rz(-1.4055777) q[3];
sx q[3];
rz(0.34077728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2365504) q[0];
sx q[0];
rz(-1.4555229) q[0];
sx q[0];
rz(-1.1849674) q[0];
rz(2.9153337) q[1];
sx q[1];
rz(-2.2626651) q[1];
sx q[1];
rz(1.0386946) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13030355) q[0];
sx q[0];
rz(-1.4892007) q[0];
sx q[0];
rz(2.5528583) q[0];
rz(-pi) q[1];
rz(-0.67579999) q[2];
sx q[2];
rz(-0.86188176) q[2];
sx q[2];
rz(1.0647286) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1490399) q[1];
sx q[1];
rz(-2.2254308) q[1];
sx q[1];
rz(0.66019411) q[1];
rz(-pi) q[2];
x q[2];
rz(0.065537621) q[3];
sx q[3];
rz(-2.0100694) q[3];
sx q[3];
rz(0.36079839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7096536) q[2];
sx q[2];
rz(-0.69810549) q[2];
sx q[2];
rz(-0.31614885) q[2];
rz(-1.6759253) q[3];
sx q[3];
rz(-1.1301872) q[3];
sx q[3];
rz(-1.7620311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6370711) q[0];
sx q[0];
rz(-0.9032473) q[0];
sx q[0];
rz(2.0380518) q[0];
rz(-2.6612813) q[1];
sx q[1];
rz(-0.66245285) q[1];
sx q[1];
rz(0.59741098) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99973122) q[0];
sx q[0];
rz(-1.7787361) q[0];
sx q[0];
rz(-2.9491762) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74540794) q[2];
sx q[2];
rz(-1.7495724) q[2];
sx q[2];
rz(2.488236) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.5525517) q[1];
sx q[1];
rz(-1.3197761) q[1];
sx q[1];
rz(1.4000721) q[1];
rz(1.1998981) q[3];
sx q[3];
rz(-1.4850332) q[3];
sx q[3];
rz(0.24468064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.20208134) q[2];
sx q[2];
rz(-1.6263522) q[2];
sx q[2];
rz(2.0651979) q[2];
rz(1.9988029) q[3];
sx q[3];
rz(-0.79311526) q[3];
sx q[3];
rz(0.49016652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6550605) q[0];
sx q[0];
rz(-1.158411) q[0];
sx q[0];
rz(3.1031188) q[0];
rz(-0.064727457) q[1];
sx q[1];
rz(-1.3836626) q[1];
sx q[1];
rz(0.23385349) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7215639) q[0];
sx q[0];
rz(-1.7917243) q[0];
sx q[0];
rz(1.3807644) q[0];
rz(-0.31189708) q[2];
sx q[2];
rz(-1.9210749) q[2];
sx q[2];
rz(-0.52679449) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.95215248) q[1];
sx q[1];
rz(-2.2169211) q[1];
sx q[1];
rz(1.6829856) q[1];
rz(-1.6285628) q[3];
sx q[3];
rz(-1.4411297) q[3];
sx q[3];
rz(-0.71021375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6528299) q[2];
sx q[2];
rz(-1.5583928) q[2];
sx q[2];
rz(-2.5433507) q[2];
rz(-0.11387842) q[3];
sx q[3];
rz(-1.7555534) q[3];
sx q[3];
rz(2.2940476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4050196) q[0];
sx q[0];
rz(-0.84252715) q[0];
sx q[0];
rz(-2.8952428) q[0];
rz(1.8440638) q[1];
sx q[1];
rz(-1.1187226) q[1];
sx q[1];
rz(-2.6447703) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56820541) q[0];
sx q[0];
rz(-2.3620689) q[0];
sx q[0];
rz(0.60978344) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70513983) q[2];
sx q[2];
rz(-1.459483) q[2];
sx q[2];
rz(1.4347347) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0223479) q[1];
sx q[1];
rz(-1.8059219) q[1];
sx q[1];
rz(1.0757331) q[1];
rz(-0.15487352) q[3];
sx q[3];
rz(-1.0836156) q[3];
sx q[3];
rz(0.39184141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42314998) q[2];
sx q[2];
rz(-2.6079874) q[2];
sx q[2];
rz(-1.9971087) q[2];
rz(-1.1540958) q[3];
sx q[3];
rz(-1.4242947) q[3];
sx q[3];
rz(-0.19788876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6941187) q[0];
sx q[0];
rz(-0.23451528) q[0];
sx q[0];
rz(3.0754454) q[0];
rz(-1.1478395) q[1];
sx q[1];
rz(-1.3775974) q[1];
sx q[1];
rz(0.60417169) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2691374) q[0];
sx q[0];
rz(-2.7949484) q[0];
sx q[0];
rz(-0.64328648) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2653052) q[2];
sx q[2];
rz(-0.54232208) q[2];
sx q[2];
rz(3.0408183) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79163247) q[1];
sx q[1];
rz(-0.59795982) q[1];
sx q[1];
rz(-0.34611846) q[1];
rz(-pi) q[2];
rz(-0.6571397) q[3];
sx q[3];
rz(-2.1517188) q[3];
sx q[3];
rz(0.37347735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8784647) q[2];
sx q[2];
rz(-2.2915514) q[2];
sx q[2];
rz(-2.7086332) q[2];
rz(1.912502) q[3];
sx q[3];
rz(-1.9357598) q[3];
sx q[3];
rz(1.8142726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94806725) q[0];
sx q[0];
rz(-2.0663517) q[0];
sx q[0];
rz(1.8684335) q[0];
rz(2.676447) q[1];
sx q[1];
rz(-1.7756614) q[1];
sx q[1];
rz(0.21496162) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8037655) q[0];
sx q[0];
rz(-1.318232) q[0];
sx q[0];
rz(0.65852491) q[0];
rz(-2.6097492) q[2];
sx q[2];
rz(-1.1924679) q[2];
sx q[2];
rz(-1.391562) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.82483208) q[1];
sx q[1];
rz(-0.90710282) q[1];
sx q[1];
rz(0.65726991) q[1];
x q[2];
rz(-0.56193476) q[3];
sx q[3];
rz(-1.7896381) q[3];
sx q[3];
rz(0.17499017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8421858) q[2];
sx q[2];
rz(-2.3278548) q[2];
sx q[2];
rz(2.1827533) q[2];
rz(-1.1622608) q[3];
sx q[3];
rz(-1.4872888) q[3];
sx q[3];
rz(1.6963262) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6017629) q[0];
sx q[0];
rz(-1.5400664) q[0];
sx q[0];
rz(-1.6590317) q[0];
rz(-2.4330347) q[1];
sx q[1];
rz(-2.9500912) q[1];
sx q[1];
rz(-0.74831829) q[1];
rz(1.2461927) q[2];
sx q[2];
rz(-2.1094252) q[2];
sx q[2];
rz(-0.88627041) q[2];
rz(0.21227588) q[3];
sx q[3];
rz(-0.73321453) q[3];
sx q[3];
rz(-1.1238255) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
