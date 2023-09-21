OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.87941909) q[0];
sx q[0];
rz(4.8882422) q[0];
sx q[0];
rz(12.56765) q[0];
rz(1.4446422) q[1];
sx q[1];
rz(-1.1029707) q[1];
sx q[1];
rz(-0.77494088) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9026731) q[0];
sx q[0];
rz(-0.15107778) q[0];
sx q[0];
rz(0.33688776) q[0];
x q[1];
rz(2.2968729) q[2];
sx q[2];
rz(-1.2436927) q[2];
sx q[2];
rz(-2.6127882) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1577658) q[1];
sx q[1];
rz(-2.9938934) q[1];
sx q[1];
rz(-1.2965135) q[1];
rz(-pi) q[2];
rz(-2.8567021) q[3];
sx q[3];
rz(-1.4244411) q[3];
sx q[3];
rz(2.871606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0455735) q[2];
sx q[2];
rz(-2.0298268) q[2];
sx q[2];
rz(1.9457031) q[2];
rz(1.1536417) q[3];
sx q[3];
rz(-2.3524645) q[3];
sx q[3];
rz(1.7529863) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7213223) q[0];
sx q[0];
rz(-2.7852311) q[0];
sx q[0];
rz(1.2715682) q[0];
rz(2.0416073) q[1];
sx q[1];
rz(-2.0309235) q[1];
sx q[1];
rz(-1.3756479) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7398867) q[0];
sx q[0];
rz(-1.6580083) q[0];
sx q[0];
rz(2.7313822) q[0];
x q[1];
rz(1.7550811) q[2];
sx q[2];
rz(-0.24225907) q[2];
sx q[2];
rz(-0.87528961) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8010506) q[1];
sx q[1];
rz(-1.2060296) q[1];
sx q[1];
rz(-1.3786475) q[1];
rz(-0.35956412) q[3];
sx q[3];
rz(-2.095795) q[3];
sx q[3];
rz(1.1121225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9282844) q[2];
sx q[2];
rz(-1.7958612) q[2];
sx q[2];
rz(2.6518872) q[2];
rz(1.1335763) q[3];
sx q[3];
rz(-2.9462892) q[3];
sx q[3];
rz(1.0666696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60004822) q[0];
sx q[0];
rz(-1.8627889) q[0];
sx q[0];
rz(-2.9009853) q[0];
rz(-2.799017) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(1.906357) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9440207) q[0];
sx q[0];
rz(-0.70957843) q[0];
sx q[0];
rz(2.5463085) q[0];
rz(1.8232934) q[2];
sx q[2];
rz(-1.9324979) q[2];
sx q[2];
rz(-0.30644882) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0864799) q[1];
sx q[1];
rz(-1.8243454) q[1];
sx q[1];
rz(0.26992814) q[1];
x q[2];
rz(-2.0577621) q[3];
sx q[3];
rz(-2.1420797) q[3];
sx q[3];
rz(0.42698241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3774595) q[2];
sx q[2];
rz(-1.5185792) q[2];
sx q[2];
rz(1.4366478) q[2];
rz(1.7403587) q[3];
sx q[3];
rz(-1.876372) q[3];
sx q[3];
rz(2.5333372) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075994611) q[0];
sx q[0];
rz(-1.6125212) q[0];
sx q[0];
rz(2.3377989) q[0];
rz(0.94961387) q[1];
sx q[1];
rz(-1.6751553) q[1];
sx q[1];
rz(0.11985699) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6896497) q[0];
sx q[0];
rz(-1.6065238) q[0];
sx q[0];
rz(2.1769051) q[0];
x q[1];
rz(-1.7159963) q[2];
sx q[2];
rz(-1.928066) q[2];
sx q[2];
rz(-1.6574588) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.90384342) q[1];
sx q[1];
rz(-1.8029873) q[1];
sx q[1];
rz(-0.018149908) q[1];
rz(-pi) q[2];
rz(2.1060137) q[3];
sx q[3];
rz(-1.03973) q[3];
sx q[3];
rz(-2.5553184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63311657) q[2];
sx q[2];
rz(-1.2458331) q[2];
sx q[2];
rz(0.072337739) q[2];
rz(-2.7667601) q[3];
sx q[3];
rz(-2.4852677) q[3];
sx q[3];
rz(1.1981296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6937834) q[0];
sx q[0];
rz(-0.57149514) q[0];
sx q[0];
rz(0.48686349) q[0];
rz(2.411719) q[1];
sx q[1];
rz(-0.90881538) q[1];
sx q[1];
rz(-1.9015076) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70532521) q[0];
sx q[0];
rz(-1.3316532) q[0];
sx q[0];
rz(1.5843841) q[0];
rz(1.8334332) q[2];
sx q[2];
rz(-1.8034168) q[2];
sx q[2];
rz(2.5831985) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.005849) q[1];
sx q[1];
rz(-1.504335) q[1];
sx q[1];
rz(-0.60484109) q[1];
x q[2];
rz(-0.9661071) q[3];
sx q[3];
rz(-1.1629512) q[3];
sx q[3];
rz(0.84433489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1029677) q[2];
sx q[2];
rz(-2.2323148) q[2];
sx q[2];
rz(0.70303482) q[2];
rz(-1.7317584) q[3];
sx q[3];
rz(-1.3491646) q[3];
sx q[3];
rz(2.9838802) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1773961) q[0];
sx q[0];
rz(-1.8494158) q[0];
sx q[0];
rz(-2.1512206) q[0];
rz(3.0888427) q[1];
sx q[1];
rz(-2.2149448) q[1];
sx q[1];
rz(-1.6606768) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71125644) q[0];
sx q[0];
rz(-2.7529786) q[0];
sx q[0];
rz(-1.7646164) q[0];
rz(1.125607) q[2];
sx q[2];
rz(-1.5700969) q[2];
sx q[2];
rz(-0.47847754) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.78696886) q[1];
sx q[1];
rz(-2.0957392) q[1];
sx q[1];
rz(0.24137361) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2495066) q[3];
sx q[3];
rz(-1.6086372) q[3];
sx q[3];
rz(0.91176646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0084373077) q[2];
sx q[2];
rz(-1.6494273) q[2];
sx q[2];
rz(2.5793502) q[2];
rz(-1.0605313) q[3];
sx q[3];
rz(-0.74354592) q[3];
sx q[3];
rz(-0.26091584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9437207) q[0];
sx q[0];
rz(-1.3523538) q[0];
sx q[0];
rz(1.6725756) q[0];
rz(-2.127227) q[1];
sx q[1];
rz(-1.0275774) q[1];
sx q[1];
rz(-1.3247103) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6727407) q[0];
sx q[0];
rz(-1.6981089) q[0];
sx q[0];
rz(-1.5936113) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5092233) q[2];
sx q[2];
rz(-1.0668796) q[2];
sx q[2];
rz(3.0798806) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.054338) q[1];
sx q[1];
rz(-1.9458276) q[1];
sx q[1];
rz(0.43941984) q[1];
rz(-pi) q[2];
rz(-1.0876098) q[3];
sx q[3];
rz(-2.4251068) q[3];
sx q[3];
rz(1.0928175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5614732) q[2];
sx q[2];
rz(-1.3886398) q[2];
sx q[2];
rz(2.6980147) q[2];
rz(0.94868547) q[3];
sx q[3];
rz(-1.1921927) q[3];
sx q[3];
rz(-0.77478066) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.401944) q[0];
sx q[0];
rz(-1.0111324) q[0];
sx q[0];
rz(2.6810714) q[0];
rz(0.1000239) q[1];
sx q[1];
rz(-0.99383751) q[1];
sx q[1];
rz(1.2896279) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3562718) q[0];
sx q[0];
rz(-0.7542146) q[0];
sx q[0];
rz(-0.69099364) q[0];
x q[1];
rz(-0.81823924) q[2];
sx q[2];
rz(-2.1365676) q[2];
sx q[2];
rz(-0.73275369) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3813078) q[1];
sx q[1];
rz(-0.21796255) q[1];
sx q[1];
rz(1.0337619) q[1];
rz(-2.4786948) q[3];
sx q[3];
rz(-1.7015966) q[3];
sx q[3];
rz(0.42636426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.148968) q[2];
sx q[2];
rz(-0.31105369) q[2];
sx q[2];
rz(2.0588493) q[2];
rz(-0.096171245) q[3];
sx q[3];
rz(-1.7008737) q[3];
sx q[3];
rz(1.2307897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37725317) q[0];
sx q[0];
rz(-0.23582533) q[0];
sx q[0];
rz(-0.33426958) q[0];
rz(-1.224068) q[1];
sx q[1];
rz(-1.5806438) q[1];
sx q[1];
rz(-0.28265488) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3960421) q[0];
sx q[0];
rz(-0.96691416) q[0];
sx q[0];
rz(-2.1349483) q[0];
rz(-1.3547782) q[2];
sx q[2];
rz(-1.2470494) q[2];
sx q[2];
rz(0.87460364) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.64687318) q[1];
sx q[1];
rz(-2.6962526) q[1];
sx q[1];
rz(-0.41782197) q[1];
rz(-pi) q[2];
rz(-0.11741365) q[3];
sx q[3];
rz(-2.252929) q[3];
sx q[3];
rz(-2.5168602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9294372) q[2];
sx q[2];
rz(-1.1634469) q[2];
sx q[2];
rz(2.4003417) q[2];
rz(-0.50179982) q[3];
sx q[3];
rz(-1.3391677) q[3];
sx q[3];
rz(0.58399502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0989477) q[0];
sx q[0];
rz(-2.2462923) q[0];
sx q[0];
rz(-2.6328971) q[0];
rz(-0.11518654) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(2.4597816) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1334575) q[0];
sx q[0];
rz(-1.4032161) q[0];
sx q[0];
rz(2.125678) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82294686) q[2];
sx q[2];
rz(-1.032885) q[2];
sx q[2];
rz(-3.092098) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.21916418) q[1];
sx q[1];
rz(-1.0370967) q[1];
sx q[1];
rz(-2.1283172) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82716771) q[3];
sx q[3];
rz(-2.4359772) q[3];
sx q[3];
rz(3.0527715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.29601413) q[2];
sx q[2];
rz(-1.6342376) q[2];
sx q[2];
rz(0.34109035) q[2];
rz(2.0579445) q[3];
sx q[3];
rz(-0.79569474) q[3];
sx q[3];
rz(0.17102374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1682128) q[0];
sx q[0];
rz(-1.4773049) q[0];
sx q[0];
rz(2.1422577) q[0];
rz(0.60733168) q[1];
sx q[1];
rz(-0.9098396) q[1];
sx q[1];
rz(0.20914016) q[1];
rz(-2.4773459) q[2];
sx q[2];
rz(-1.1462117) q[2];
sx q[2];
rz(0.29427634) q[2];
rz(2.7568983) q[3];
sx q[3];
rz(-1.7582498) q[3];
sx q[3];
rz(2.869217) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];