OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.38874415) q[0];
sx q[0];
rz(3.677877) q[0];
sx q[0];
rz(10.372547) q[0];
rz(-1.3287969) q[1];
sx q[1];
rz(-1.8741908) q[1];
sx q[1];
rz(1.0277494) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0004814) q[0];
sx q[0];
rz(-1.1950462) q[0];
sx q[0];
rz(-1.5462589) q[0];
x q[1];
rz(-0.22123863) q[2];
sx q[2];
rz(-1.0916296) q[2];
sx q[2];
rz(-2.0146807) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7386908) q[1];
sx q[1];
rz(-1.9539023) q[1];
sx q[1];
rz(-2.4163567) q[1];
rz(2.5563452) q[3];
sx q[3];
rz(-2.4260776) q[3];
sx q[3];
rz(-1.9714718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0119005) q[2];
sx q[2];
rz(-1.4346069) q[2];
sx q[2];
rz(-2.091308) q[2];
rz(-1.1132647) q[3];
sx q[3];
rz(-2.2498825) q[3];
sx q[3];
rz(-3.0734857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.072409078) q[0];
sx q[0];
rz(-1.8962815) q[0];
sx q[0];
rz(-2.8438399) q[0];
rz(-0.61966664) q[1];
sx q[1];
rz(-1.0071808) q[1];
sx q[1];
rz(-1.108095) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71582687) q[0];
sx q[0];
rz(-2.1799488) q[0];
sx q[0];
rz(-2.8858375) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96252243) q[2];
sx q[2];
rz(-1.2108742) q[2];
sx q[2];
rz(-2.4134709) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.019921692) q[1];
sx q[1];
rz(-0.53598511) q[1];
sx q[1];
rz(-2.3977445) q[1];
x q[2];
rz(0.70257367) q[3];
sx q[3];
rz(-1.3630023) q[3];
sx q[3];
rz(2.7489565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.46253282) q[2];
sx q[2];
rz(-2.1647537) q[2];
sx q[2];
rz(-0.97529808) q[2];
rz(-0.9179999) q[3];
sx q[3];
rz(-1.8564329) q[3];
sx q[3];
rz(0.29618922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.179203) q[0];
sx q[0];
rz(-2.2900892) q[0];
sx q[0];
rz(2.5986824) q[0];
rz(2.2593598) q[1];
sx q[1];
rz(-2.0062607) q[1];
sx q[1];
rz(0.96484819) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1436413) q[0];
sx q[0];
rz(-0.43295857) q[0];
sx q[0];
rz(2.4233682) q[0];
rz(-pi) q[1];
rz(-1.6134878) q[2];
sx q[2];
rz(-1.8674769) q[2];
sx q[2];
rz(1.3670849) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.94273401) q[1];
sx q[1];
rz(-3.004289) q[1];
sx q[1];
rz(1.0820461) q[1];
rz(0.19823719) q[3];
sx q[3];
rz(-1.1612411) q[3];
sx q[3];
rz(-0.73892456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0470011) q[2];
sx q[2];
rz(-0.61085218) q[2];
sx q[2];
rz(-1.1331406) q[2];
rz(2.9099693) q[3];
sx q[3];
rz(-1.8685721) q[3];
sx q[3];
rz(-0.75727063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8621181) q[0];
sx q[0];
rz(-0.010443895) q[0];
sx q[0];
rz(1.7650771) q[0];
rz(0.51849413) q[1];
sx q[1];
rz(-1.8771749) q[1];
sx q[1];
rz(-2.8994697) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74894416) q[0];
sx q[0];
rz(-2.2376275) q[0];
sx q[0];
rz(1.3632266) q[0];
rz(1.5577199) q[2];
sx q[2];
rz(-1.1411975) q[2];
sx q[2];
rz(-1.0193046) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0603234) q[1];
sx q[1];
rz(-0.4821061) q[1];
sx q[1];
rz(-0.37270765) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4297156) q[3];
sx q[3];
rz(-1.9851079) q[3];
sx q[3];
rz(-2.8431818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.018192856) q[2];
sx q[2];
rz(-2.172956) q[2];
sx q[2];
rz(0.094853178) q[2];
rz(1.799396) q[3];
sx q[3];
rz(-1.3972524) q[3];
sx q[3];
rz(-2.7024787) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3577394) q[0];
sx q[0];
rz(-1.3107603) q[0];
sx q[0];
rz(-0.36079303) q[0];
rz(1.7533253) q[1];
sx q[1];
rz(-1.3307945) q[1];
sx q[1];
rz(-2.0070019) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9103694) q[0];
sx q[0];
rz(-2.0632422) q[0];
sx q[0];
rz(-0.72427303) q[0];
rz(-pi) q[1];
rz(-1.1410494) q[2];
sx q[2];
rz(-1.8998713) q[2];
sx q[2];
rz(-2.3171901) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6280917) q[1];
sx q[1];
rz(-2.6673632) q[1];
sx q[1];
rz(-0.34631108) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4868823) q[3];
sx q[3];
rz(-2.0062371) q[3];
sx q[3];
rz(0.30121379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0052884) q[2];
sx q[2];
rz(-1.4207999) q[2];
sx q[2];
rz(-0.57265442) q[2];
rz(-2.2128361) q[3];
sx q[3];
rz(-2.6199665) q[3];
sx q[3];
rz(-1.1675534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8144433) q[0];
sx q[0];
rz(-2.0454018) q[0];
sx q[0];
rz(1.8922528) q[0];
rz(1.2231187) q[1];
sx q[1];
rz(-1.616281) q[1];
sx q[1];
rz(-1.1522326) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.457068) q[0];
sx q[0];
rz(-3.0355434) q[0];
sx q[0];
rz(-1.4507136) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4695775) q[2];
sx q[2];
rz(-0.61215559) q[2];
sx q[2];
rz(1.9415346) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7631543) q[1];
sx q[1];
rz(-1.3707146) q[1];
sx q[1];
rz(1.4640019) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4310962) q[3];
sx q[3];
rz(-2.3264255) q[3];
sx q[3];
rz(1.1991771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6727009) q[2];
sx q[2];
rz(-1.2142618) q[2];
sx q[2];
rz(-0.99096283) q[2];
rz(0.64783603) q[3];
sx q[3];
rz(-2.1760553) q[3];
sx q[3];
rz(-0.34753862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25794849) q[0];
sx q[0];
rz(-0.22709665) q[0];
sx q[0];
rz(-0.062285475) q[0];
rz(-2.9557872) q[1];
sx q[1];
rz(-1.4567016) q[1];
sx q[1];
rz(-2.7468162) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61790723) q[0];
sx q[0];
rz(-2.6751408) q[0];
sx q[0];
rz(2.594069) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4418418) q[2];
sx q[2];
rz(-0.47669461) q[2];
sx q[2];
rz(-1.9112019) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2784087) q[1];
sx q[1];
rz(-1.3235958) q[1];
sx q[1];
rz(-2.7017038) q[1];
rz(-pi) q[2];
rz(2.9069101) q[3];
sx q[3];
rz(-2.9838786) q[3];
sx q[3];
rz(-1.371067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.64849598) q[2];
sx q[2];
rz(-2.811537) q[2];
sx q[2];
rz(0.27302343) q[2];
rz(-1.3027044) q[3];
sx q[3];
rz(-1.8282993) q[3];
sx q[3];
rz(0.31204143) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39772314) q[0];
sx q[0];
rz(-1.1345154) q[0];
sx q[0];
rz(1.8564818) q[0];
rz(1.6400736) q[1];
sx q[1];
rz(-1.7506426) q[1];
sx q[1];
rz(-1.8008908) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2368603) q[0];
sx q[0];
rz(-1.245984) q[0];
sx q[0];
rz(1.9843285) q[0];
x q[1];
rz(0.69581823) q[2];
sx q[2];
rz(-1.6098445) q[2];
sx q[2];
rz(0.85862904) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.50780523) q[1];
sx q[1];
rz(-1.2101189) q[1];
sx q[1];
rz(2.7651869) q[1];
x q[2];
rz(2.0128485) q[3];
sx q[3];
rz(-1.8574517) q[3];
sx q[3];
rz(-3.0451881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0354707) q[2];
sx q[2];
rz(-2.3351228) q[2];
sx q[2];
rz(1.0236053) q[2];
rz(-2.9566531) q[3];
sx q[3];
rz(-2.7513294) q[3];
sx q[3];
rz(-2.8997054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0446562) q[0];
sx q[0];
rz(-0.99656314) q[0];
sx q[0];
rz(-1.5203083) q[0];
rz(2.8114491) q[1];
sx q[1];
rz(-1.9338927) q[1];
sx q[1];
rz(-0.83713371) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1820113) q[0];
sx q[0];
rz(-0.89996979) q[0];
sx q[0];
rz(-1.3120033) q[0];
rz(-pi) q[1];
rz(-1.8673973) q[2];
sx q[2];
rz(-2.1269848) q[2];
sx q[2];
rz(-2.0402758) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.21048927) q[1];
sx q[1];
rz(-1.466202) q[1];
sx q[1];
rz(1.4180257) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0666396) q[3];
sx q[3];
rz(-1.3131724) q[3];
sx q[3];
rz(-2.6356634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5403486) q[2];
sx q[2];
rz(-2.0566172) q[2];
sx q[2];
rz(-2.9619651) q[2];
rz(-0.99572292) q[3];
sx q[3];
rz(-1.2487753) q[3];
sx q[3];
rz(-1.8306336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3778465) q[0];
sx q[0];
rz(-2.7959931) q[0];
sx q[0];
rz(2.0843704) q[0];
rz(3.0341042) q[1];
sx q[1];
rz(-1.2534393) q[1];
sx q[1];
rz(2.1616139) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36469034) q[0];
sx q[0];
rz(-1.5199465) q[0];
sx q[0];
rz(-1.6965673) q[0];
rz(-pi) q[1];
rz(-2.6238742) q[2];
sx q[2];
rz(-1.1283518) q[2];
sx q[2];
rz(-0.25416086) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.500463) q[1];
sx q[1];
rz(-0.49912057) q[1];
sx q[1];
rz(-1.0541381) q[1];
rz(-1.2471334) q[3];
sx q[3];
rz(-1.666288) q[3];
sx q[3];
rz(-3.132706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8367299) q[2];
sx q[2];
rz(-1.9367846) q[2];
sx q[2];
rz(-2.1137962) q[2];
rz(-1.7547539) q[3];
sx q[3];
rz(-1.8324865) q[3];
sx q[3];
rz(0.28361472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86823157) q[0];
sx q[0];
rz(-1.0366476) q[0];
sx q[0];
rz(-1.114053) q[0];
rz(-0.83203075) q[1];
sx q[1];
rz(-2.6770626) q[1];
sx q[1];
rz(-2.4774036) q[1];
rz(-1.8503415) q[2];
sx q[2];
rz(-2.5646979) q[2];
sx q[2];
rz(-0.23451351) q[2];
rz(2.6408623) q[3];
sx q[3];
rz(-1.9095608) q[3];
sx q[3];
rz(0.6469971) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
