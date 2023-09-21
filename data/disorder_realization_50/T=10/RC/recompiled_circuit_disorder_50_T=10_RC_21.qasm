OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2621736) q[0];
sx q[0];
rz(-1.7466495) q[0];
sx q[0];
rz(3.1403132) q[0];
rz(-1.6969504) q[1];
sx q[1];
rz(4.2445634) q[1];
sx q[1];
rz(7.0581262) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89844184) q[0];
sx q[0];
rz(-1.713322) q[0];
sx q[0];
rz(-1.5205163) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0432267) q[2];
sx q[2];
rz(-2.35765) q[2];
sx q[2];
rz(0.69477615) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.98382681) q[1];
sx q[1];
rz(-0.14769927) q[1];
sx q[1];
rz(1.2965135) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8567021) q[3];
sx q[3];
rz(-1.4244411) q[3];
sx q[3];
rz(-0.26998664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0455735) q[2];
sx q[2];
rz(-2.0298268) q[2];
sx q[2];
rz(1.9457031) q[2];
rz(1.1536417) q[3];
sx q[3];
rz(-0.78912815) q[3];
sx q[3];
rz(1.3886064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4202704) q[0];
sx q[0];
rz(-0.35636154) q[0];
sx q[0];
rz(1.8700245) q[0];
rz(1.0999854) q[1];
sx q[1];
rz(-2.0309235) q[1];
sx q[1];
rz(-1.7659448) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7398867) q[0];
sx q[0];
rz(-1.4835843) q[0];
sx q[0];
rz(0.4102104) q[0];
x q[1];
rz(0.045250821) q[2];
sx q[2];
rz(-1.8088733) q[2];
sx q[2];
rz(-1.0649875) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1609636) q[1];
sx q[1];
rz(-1.3914319) q[1];
sx q[1];
rz(0.37100002) q[1];
rz(1.0248915) q[3];
sx q[3];
rz(-0.62666577) q[3];
sx q[3];
rz(0.46862593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.21330825) q[2];
sx q[2];
rz(-1.7958612) q[2];
sx q[2];
rz(-0.48970547) q[2];
rz(1.1335763) q[3];
sx q[3];
rz(-0.19530345) q[3];
sx q[3];
rz(-1.0666696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.5415444) q[0];
sx q[0];
rz(-1.8627889) q[0];
sx q[0];
rz(-0.24060732) q[0];
rz(-2.799017) q[1];
sx q[1];
rz(-2.1668285) q[1];
sx q[1];
rz(-1.906357) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89861682) q[0];
sx q[0];
rz(-1.1968062) q[0];
sx q[0];
rz(2.5234733) q[0];
rz(-pi) q[1];
rz(-2.7690976) q[2];
sx q[2];
rz(-1.3349581) q[2];
sx q[2];
rz(1.3553938) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8894549) q[1];
sx q[1];
rz(-2.7733907) q[1];
sx q[1];
rz(2.3705269) q[1];
x q[2];
rz(-2.0577621) q[3];
sx q[3];
rz(-0.999513) q[3];
sx q[3];
rz(-0.42698241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.76413313) q[2];
sx q[2];
rz(-1.5185792) q[2];
sx q[2];
rz(1.7049449) q[2];
rz(1.4012339) q[3];
sx q[3];
rz(-1.876372) q[3];
sx q[3];
rz(-2.5333372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.065598) q[0];
sx q[0];
rz(-1.5290715) q[0];
sx q[0];
rz(-2.3377989) q[0];
rz(2.1919788) q[1];
sx q[1];
rz(-1.4664374) q[1];
sx q[1];
rz(0.11985699) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0940995) q[0];
sx q[0];
rz(-0.96512981) q[0];
sx q[0];
rz(3.0981307) q[0];
x q[1];
rz(2.7808431) q[2];
sx q[2];
rz(-1.7067688) q[2];
sx q[2];
rz(-3.0038358) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90384342) q[1];
sx q[1];
rz(-1.8029873) q[1];
sx q[1];
rz(0.018149908) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5991163) q[3];
sx q[3];
rz(-1.1154419) q[3];
sx q[3];
rz(1.276254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5084761) q[2];
sx q[2];
rz(-1.2458331) q[2];
sx q[2];
rz(-0.072337739) q[2];
rz(2.7667601) q[3];
sx q[3];
rz(-0.65632498) q[3];
sx q[3];
rz(1.1981296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-2.4478093) q[0];
sx q[0];
rz(-0.57149514) q[0];
sx q[0];
rz(2.6547292) q[0];
rz(2.411719) q[1];
sx q[1];
rz(-2.2327773) q[1];
sx q[1];
rz(-1.240085) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4935721) q[0];
sx q[0];
rz(-0.23952142) q[0];
sx q[0];
rz(-0.055672107) q[0];
rz(-2.901022) q[2];
sx q[2];
rz(-1.8261989) q[2];
sx q[2];
rz(2.0672928) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6606635) q[1];
sx q[1];
rz(-2.1741121) q[1];
sx q[1];
rz(-1.6515345) q[1];
rz(-pi) q[2];
rz(-0.9661071) q[3];
sx q[3];
rz(-1.9786414) q[3];
sx q[3];
rz(2.2972578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1029677) q[2];
sx q[2];
rz(-2.2323148) q[2];
sx q[2];
rz(2.4385578) q[2];
rz(1.7317584) q[3];
sx q[3];
rz(-1.792428) q[3];
sx q[3];
rz(-0.15771244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96419656) q[0];
sx q[0];
rz(-1.8494158) q[0];
sx q[0];
rz(-2.1512206) q[0];
rz(-3.0888427) q[1];
sx q[1];
rz(-0.92664781) q[1];
sx q[1];
rz(1.4809158) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4617417) q[0];
sx q[0];
rz(-1.6438419) q[0];
sx q[0];
rz(1.1887656) q[0];
rz(1.125607) q[2];
sx q[2];
rz(-1.5700969) q[2];
sx q[2];
rz(-0.47847754) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.33038352) q[1];
sx q[1];
rz(-0.57302176) q[1];
sx q[1];
rz(1.9622383) q[1];
x q[2];
rz(3.092993) q[3];
sx q[3];
rz(-2.2489293) q[3];
sx q[3];
rz(2.5130659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.0084373077) q[2];
sx q[2];
rz(-1.4921654) q[2];
sx q[2];
rz(0.56224242) q[2];
rz(-1.0605313) q[3];
sx q[3];
rz(-2.3980467) q[3];
sx q[3];
rz(0.26091584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9437207) q[0];
sx q[0];
rz(-1.7892388) q[0];
sx q[0];
rz(-1.4690171) q[0];
rz(-2.127227) q[1];
sx q[1];
rz(-1.0275774) q[1];
sx q[1];
rz(1.8168824) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64667386) q[0];
sx q[0];
rz(-3.0122628) q[0];
sx q[0];
rz(-2.9652251) q[0];
rz(-2.1340738) q[2];
sx q[2];
rz(-2.0118606) q[2];
sx q[2];
rz(1.8958467) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.054338) q[1];
sx q[1];
rz(-1.9458276) q[1];
sx q[1];
rz(0.43941984) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0539829) q[3];
sx q[3];
rz(-2.4251068) q[3];
sx q[3];
rz(1.0928175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5614732) q[2];
sx q[2];
rz(-1.7529528) q[2];
sx q[2];
rz(-2.6980147) q[2];
rz(0.94868547) q[3];
sx q[3];
rz(-1.9493999) q[3];
sx q[3];
rz(0.77478066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7396486) q[0];
sx q[0];
rz(-2.1304603) q[0];
sx q[0];
rz(-2.6810714) q[0];
rz(3.0415688) q[1];
sx q[1];
rz(-2.1477551) q[1];
sx q[1];
rz(-1.8519648) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063110654) q[0];
sx q[0];
rz(-2.1266298) q[0];
sx q[0];
rz(2.1102935) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81823924) q[2];
sx q[2];
rz(-2.1365676) q[2];
sx q[2];
rz(-0.73275369) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8045115) q[1];
sx q[1];
rz(-1.6816499) q[1];
sx q[1];
rz(1.7588508) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9309978) q[3];
sx q[3];
rz(-0.67376332) q[3];
sx q[3];
rz(-1.3099561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.9926247) q[2];
sx q[2];
rz(-0.31105369) q[2];
sx q[2];
rz(1.0827433) q[2];
rz(-0.096171245) q[3];
sx q[3];
rz(-1.440719) q[3];
sx q[3];
rz(1.910803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37725317) q[0];
sx q[0];
rz(-0.23582533) q[0];
sx q[0];
rz(0.33426958) q[0];
rz(1.9175247) q[1];
sx q[1];
rz(-1.5806438) q[1];
sx q[1];
rz(-0.28265488) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44336549) q[0];
sx q[0];
rz(-2.3400314) q[0];
sx q[0];
rz(-2.482224) q[0];
x q[1];
rz(0.56844175) q[2];
sx q[2];
rz(-2.754515) q[2];
sx q[2];
rz(-2.8708411) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.104056) q[1];
sx q[1];
rz(-1.1661342) q[1];
sx q[1];
rz(-1.3794823) q[1];
rz(-2.2563124) q[3];
sx q[3];
rz(-1.6618528) q[3];
sx q[3];
rz(-2.1212999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9294372) q[2];
sx q[2];
rz(-1.1634469) q[2];
sx q[2];
rz(0.7412509) q[2];
rz(0.50179982) q[3];
sx q[3];
rz(-1.8024249) q[3];
sx q[3];
rz(-2.5575976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.042645) q[0];
sx q[0];
rz(-0.89530033) q[0];
sx q[0];
rz(-0.50869554) q[0];
rz(0.11518654) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(0.68181109) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6012321) q[0];
sx q[0];
rz(-1.0245748) q[0];
sx q[0];
rz(-0.19646125) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3186458) q[2];
sx q[2];
rz(-2.1087077) q[2];
sx q[2];
rz(3.092098) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9224285) q[1];
sx q[1];
rz(-1.0370967) q[1];
sx q[1];
rz(1.0132755) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6184611) q[3];
sx q[3];
rz(-2.0683859) q[3];
sx q[3];
rz(2.350972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29601413) q[2];
sx q[2];
rz(-1.5073551) q[2];
sx q[2];
rz(-0.34109035) q[2];
rz(1.0836481) q[3];
sx q[3];
rz(-2.3458979) q[3];
sx q[3];
rz(0.17102374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9733799) q[0];
sx q[0];
rz(-1.4773049) q[0];
sx q[0];
rz(2.1422577) q[0];
rz(2.534261) q[1];
sx q[1];
rz(-2.2317531) q[1];
sx q[1];
rz(-2.9324525) q[1];
rz(-2.4773459) q[2];
sx q[2];
rz(-1.1462117) q[2];
sx q[2];
rz(0.29427634) q[2];
rz(-0.38469436) q[3];
sx q[3];
rz(-1.7582498) q[3];
sx q[3];
rz(2.869217) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
