OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10387575) q[0];
sx q[0];
rz(1.2020943) q[0];
sx q[0];
rz(10.572875) q[0];
rz(1.2530874) q[1];
sx q[1];
rz(-2.2009067) q[1];
sx q[1];
rz(1.3936477) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5210261) q[0];
sx q[0];
rz(-1.5687843) q[0];
sx q[0];
rz(-0.30962551) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68140985) q[2];
sx q[2];
rz(-2.133956) q[2];
sx q[2];
rz(-1.5208706) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4204971) q[1];
sx q[1];
rz(-1.8407397) q[1];
sx q[1];
rz(-0.47081486) q[1];
x q[2];
rz(-2.0066891) q[3];
sx q[3];
rz(-2.5385058) q[3];
sx q[3];
rz(3.1325504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.48774886) q[2];
sx q[2];
rz(-1.8493435) q[2];
sx q[2];
rz(0.1208819) q[2];
rz(0.17928784) q[3];
sx q[3];
rz(-2.5458953) q[3];
sx q[3];
rz(-2.9860935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.091846175) q[0];
sx q[0];
rz(-0.76773983) q[0];
sx q[0];
rz(3.0088186) q[0];
rz(-1.6800539) q[1];
sx q[1];
rz(-1.5613873) q[1];
sx q[1];
rz(2.9002088) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7624843) q[0];
sx q[0];
rz(-0.53775162) q[0];
sx q[0];
rz(2.3178029) q[0];
rz(-pi) q[1];
rz(-2.3953305) q[2];
sx q[2];
rz(-2.4733739) q[2];
sx q[2];
rz(-1.2506739) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6893377) q[1];
sx q[1];
rz(-0.66879767) q[1];
sx q[1];
rz(1.7178839) q[1];
rz(-1.7412132) q[3];
sx q[3];
rz(-2.9885871) q[3];
sx q[3];
rz(1.1982329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1277348) q[2];
sx q[2];
rz(-1.3188136) q[2];
sx q[2];
rz(1.1068809) q[2];
rz(-1.3876623) q[3];
sx q[3];
rz(-2.6219086) q[3];
sx q[3];
rz(-1.9096411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.36104193) q[0];
sx q[0];
rz(-1.1013958) q[0];
sx q[0];
rz(-0.74044359) q[0];
rz(-2.6904147) q[1];
sx q[1];
rz(-1.8809044) q[1];
sx q[1];
rz(2.0887451) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8448062) q[0];
sx q[0];
rz(-1.0264945) q[0];
sx q[0];
rz(-1.143572) q[0];
x q[1];
rz(-0.29849507) q[2];
sx q[2];
rz(-0.91909354) q[2];
sx q[2];
rz(-3.0534844) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9991418) q[1];
sx q[1];
rz(-1.72292) q[1];
sx q[1];
rz(-1.4494447) q[1];
rz(-1.5925519) q[3];
sx q[3];
rz(-2.7113911) q[3];
sx q[3];
rz(2.3019626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8911002) q[2];
sx q[2];
rz(-1.5414457) q[2];
sx q[2];
rz(-0.54692522) q[2];
rz(-0.28918239) q[3];
sx q[3];
rz(-0.5368036) q[3];
sx q[3];
rz(0.012399013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96994394) q[0];
sx q[0];
rz(-2.8778853) q[0];
sx q[0];
rz(-1.7893715) q[0];
rz(-0.2098473) q[1];
sx q[1];
rz(-2.4317957) q[1];
sx q[1];
rz(2.9052177) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073233152) q[0];
sx q[0];
rz(-1.7422361) q[0];
sx q[0];
rz(1.3017446) q[0];
x q[1];
rz(-0.9985853) q[2];
sx q[2];
rz(-0.59362715) q[2];
sx q[2];
rz(-1.3015391) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.216815) q[1];
sx q[1];
rz(-1.5773915) q[1];
sx q[1];
rz(1.738468) q[1];
x q[2];
rz(1.4218016) q[3];
sx q[3];
rz(-2.0553203) q[3];
sx q[3];
rz(-1.5414343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9005047) q[2];
sx q[2];
rz(-0.97854486) q[2];
sx q[2];
rz(2.5881055) q[2];
rz(-0.91529804) q[3];
sx q[3];
rz(-1.2585636) q[3];
sx q[3];
rz(0.64490157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(-2.6699162) q[0];
sx q[0];
rz(-1.8138509) q[0];
sx q[0];
rz(0.70415235) q[0];
rz(2.0856805) q[1];
sx q[1];
rz(-0.45509714) q[1];
sx q[1];
rz(-0.59590894) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51673698) q[0];
sx q[0];
rz(-1.7863818) q[0];
sx q[0];
rz(-3.0481824) q[0];
rz(-pi) q[1];
rz(-0.17275177) q[2];
sx q[2];
rz(-0.68095945) q[2];
sx q[2];
rz(1.3576042) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.51922819) q[1];
sx q[1];
rz(-1.16751) q[1];
sx q[1];
rz(1.7720023) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3711575) q[3];
sx q[3];
rz(-1.7426997) q[3];
sx q[3];
rz(2.2558444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.492505) q[2];
sx q[2];
rz(-1.1143782) q[2];
sx q[2];
rz(-2.5904783) q[2];
rz(2.9344432) q[3];
sx q[3];
rz(-1.8271577) q[3];
sx q[3];
rz(1.5392039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9838487) q[0];
sx q[0];
rz(-0.85726964) q[0];
sx q[0];
rz(-0.95170784) q[0];
rz(-2.6668008) q[1];
sx q[1];
rz(-1.1750849) q[1];
sx q[1];
rz(0.29528433) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8808522) q[0];
sx q[0];
rz(-1.7096925) q[0];
sx q[0];
rz(-1.9689346) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6906602) q[2];
sx q[2];
rz(-0.34313289) q[2];
sx q[2];
rz(-0.97230881) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30027896) q[1];
sx q[1];
rz(-2.0804555) q[1];
sx q[1];
rz(-1.0213486) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4548364) q[3];
sx q[3];
rz(-2.6280858) q[3];
sx q[3];
rz(-2.6503369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39020145) q[2];
sx q[2];
rz(-1.7753121) q[2];
sx q[2];
rz(-2.2078216) q[2];
rz(1.4592524) q[3];
sx q[3];
rz(-1.7873584) q[3];
sx q[3];
rz(-1.4565844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(3.0657601) q[0];
sx q[0];
rz(-1.1446784) q[0];
sx q[0];
rz(-0.24205762) q[0];
rz(0.66479713) q[1];
sx q[1];
rz(-1.2882065) q[1];
sx q[1];
rz(2.738293) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86352578) q[0];
sx q[0];
rz(-0.51260548) q[0];
sx q[0];
rz(-2.313732) q[0];
rz(1.2861916) q[2];
sx q[2];
rz(-0.91214123) q[2];
sx q[2];
rz(-1.8919485) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.51822) q[1];
sx q[1];
rz(-1.5726591) q[1];
sx q[1];
rz(-1.3606447) q[1];
x q[2];
rz(-0.025267406) q[3];
sx q[3];
rz(-1.7606887) q[3];
sx q[3];
rz(0.59786284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2287801) q[2];
sx q[2];
rz(-2.0704724) q[2];
sx q[2];
rz(0.43441233) q[2];
rz(-2.1408391) q[3];
sx q[3];
rz(-0.36608168) q[3];
sx q[3];
rz(2.6312857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-3.1338761) q[0];
sx q[0];
rz(-0.0061329734) q[0];
sx q[0];
rz(0.49466053) q[0];
rz(1.6330632) q[1];
sx q[1];
rz(-1.5155019) q[1];
sx q[1];
rz(2.5411434) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52727115) q[0];
sx q[0];
rz(-1.1724171) q[0];
sx q[0];
rz(0.50171731) q[0];
rz(-pi) q[1];
rz(3.1152994) q[2];
sx q[2];
rz(-1.0458938) q[2];
sx q[2];
rz(-2.3412447) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.32158347) q[1];
sx q[1];
rz(-1.4516186) q[1];
sx q[1];
rz(2.9585341) q[1];
rz(1.7793711) q[3];
sx q[3];
rz(-1.8684636) q[3];
sx q[3];
rz(2.0901594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1373458) q[2];
sx q[2];
rz(-0.9937976) q[2];
sx q[2];
rz(-0.11631575) q[2];
rz(-2.7159193) q[3];
sx q[3];
rz(-2.1843572) q[3];
sx q[3];
rz(-11*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15329926) q[0];
sx q[0];
rz(-2.9635552) q[0];
sx q[0];
rz(-1.6631888) q[0];
rz(-0.93961811) q[1];
sx q[1];
rz(-1.3213108) q[1];
sx q[1];
rz(0.41752648) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.674304) q[0];
sx q[0];
rz(-2.1018565) q[0];
sx q[0];
rz(-0.051961016) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1065503) q[2];
sx q[2];
rz(-1.8333734) q[2];
sx q[2];
rz(0.1072466) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.49839834) q[1];
sx q[1];
rz(-2.8865951) q[1];
sx q[1];
rz(-0.71868371) q[1];
rz(-pi) q[2];
rz(2.0017654) q[3];
sx q[3];
rz(-1.7472072) q[3];
sx q[3];
rz(-2.924502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4801165) q[2];
sx q[2];
rz(-2.5051703) q[2];
sx q[2];
rz(-2.647906) q[2];
rz(0.72475973) q[3];
sx q[3];
rz(-1.1586435) q[3];
sx q[3];
rz(0.31989583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9676554) q[0];
sx q[0];
rz(-2.4854361) q[0];
sx q[0];
rz(2.4560112) q[0];
rz(2.8441692) q[1];
sx q[1];
rz(-0.23700266) q[1];
sx q[1];
rz(-1.1313653) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6146651) q[0];
sx q[0];
rz(-2.1242495) q[0];
sx q[0];
rz(1.1620031) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70558833) q[2];
sx q[2];
rz(-1.0552647) q[2];
sx q[2];
rz(2.3956092) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8611421) q[1];
sx q[1];
rz(-0.44499731) q[1];
sx q[1];
rz(0.54235561) q[1];
rz(-pi) q[2];
rz(0.60641842) q[3];
sx q[3];
rz(-1.1675721) q[3];
sx q[3];
rz(1.2244146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3161105) q[2];
sx q[2];
rz(-1.1967412) q[2];
sx q[2];
rz(-2.4342009) q[2];
rz(2.3317544) q[3];
sx q[3];
rz(-0.67088586) q[3];
sx q[3];
rz(2.9530318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0903044) q[0];
sx q[0];
rz(-1.123883) q[0];
sx q[0];
rz(-0.71258769) q[0];
rz(-0.21223016) q[1];
sx q[1];
rz(-1.6925015) q[1];
sx q[1];
rz(-0.5136516) q[1];
rz(-0.60317294) q[2];
sx q[2];
rz(-2.0576253) q[2];
sx q[2];
rz(-1.1228592) q[2];
rz(2.2255696) q[3];
sx q[3];
rz(-2.3229204) q[3];
sx q[3];
rz(-1.7182072) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
