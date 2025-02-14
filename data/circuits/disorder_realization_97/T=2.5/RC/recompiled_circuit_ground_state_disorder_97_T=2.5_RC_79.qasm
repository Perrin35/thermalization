OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.41807732) q[0];
sx q[0];
rz(-0.014208566) q[0];
sx q[0];
rz(3.1007822) q[0];
rz(3.1006676) q[1];
sx q[1];
rz(-2.0857781) q[1];
sx q[1];
rz(-2.5850886) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41123996) q[0];
sx q[0];
rz(-1.7539193) q[0];
sx q[0];
rz(1.9224115) q[0];
rz(0.5719127) q[2];
sx q[2];
rz(-1.3777857) q[2];
sx q[2];
rz(0.62096769) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5750016) q[1];
sx q[1];
rz(-2.0739158) q[1];
sx q[1];
rz(1.9921705) q[1];
x q[2];
rz(-1.1993042) q[3];
sx q[3];
rz(-2.8879734) q[3];
sx q[3];
rz(-2.729148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9875662) q[2];
sx q[2];
rz(-1.3904927) q[2];
sx q[2];
rz(1.7517368) q[2];
rz(-3.0912345) q[3];
sx q[3];
rz(-1.4202838) q[3];
sx q[3];
rz(2.0462947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9339226) q[0];
sx q[0];
rz(-1.727513) q[0];
sx q[0];
rz(-0.06047824) q[0];
rz(-2.1531847) q[1];
sx q[1];
rz(-1.6840839) q[1];
sx q[1];
rz(-2.9284533) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7558862) q[0];
sx q[0];
rz(-1.1637628) q[0];
sx q[0];
rz(-2.4197121) q[0];
rz(-pi) q[1];
rz(1.408281) q[2];
sx q[2];
rz(-1.4569511) q[2];
sx q[2];
rz(1.3448037) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2717183) q[1];
sx q[1];
rz(-1.8228181) q[1];
sx q[1];
rz(1.6099754) q[1];
rz(-1.8993103) q[3];
sx q[3];
rz(-0.51988784) q[3];
sx q[3];
rz(-2.0175135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.1468143) q[2];
sx q[2];
rz(-1.549336) q[2];
sx q[2];
rz(3.1352622) q[2];
rz(-1.7019042) q[3];
sx q[3];
rz(-2.3531745) q[3];
sx q[3];
rz(-2.7009713) q[3];
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
rz(-0.45563662) q[0];
sx q[0];
rz(-0.39304471) q[0];
sx q[0];
rz(-1.8259557) q[0];
rz(1.3872967) q[1];
sx q[1];
rz(-2.5376814) q[1];
sx q[1];
rz(-0.046262892) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10704087) q[0];
sx q[0];
rz(-3.1116293) q[0];
sx q[0];
rz(-3.1065123) q[0];
x q[1];
rz(-2.6753083) q[2];
sx q[2];
rz(-1.2215677) q[2];
sx q[2];
rz(3.0602855) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8475593) q[1];
sx q[1];
rz(-0.98638505) q[1];
sx q[1];
rz(1.1520034) q[1];
rz(0.1667695) q[3];
sx q[3];
rz(-0.96264364) q[3];
sx q[3];
rz(1.2603354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5578654) q[2];
sx q[2];
rz(-1.9670308) q[2];
sx q[2];
rz(-0.91090703) q[2];
rz(1.2298498) q[3];
sx q[3];
rz(-0.76084796) q[3];
sx q[3];
rz(2.7228444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8987027) q[0];
sx q[0];
rz(-1.3985343) q[0];
sx q[0];
rz(-2.718495) q[0];
rz(2.2078919) q[1];
sx q[1];
rz(-1.2010778) q[1];
sx q[1];
rz(2.7074714) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98744623) q[0];
sx q[0];
rz(-2.848296) q[0];
sx q[0];
rz(-3.0213839) q[0];
rz(-pi) q[1];
rz(2.9387596) q[2];
sx q[2];
rz(-2.4713608) q[2];
sx q[2];
rz(0.46593004) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0980743) q[1];
sx q[1];
rz(-1.4427818) q[1];
sx q[1];
rz(-0.56358595) q[1];
rz(1.36606) q[3];
sx q[3];
rz(-1.4582485) q[3];
sx q[3];
rz(1.1074668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.18647974) q[2];
sx q[2];
rz(-0.13804144) q[2];
sx q[2];
rz(-1.0008) q[2];
rz(-2.3872088) q[3];
sx q[3];
rz(-2.3013134) q[3];
sx q[3];
rz(1.3365356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56483018) q[0];
sx q[0];
rz(-2.5046528) q[0];
sx q[0];
rz(-1.7621967) q[0];
rz(2.5881536) q[1];
sx q[1];
rz(-1.484) q[1];
sx q[1];
rz(2.4310506) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64122771) q[0];
sx q[0];
rz(-2.3661748) q[0];
sx q[0];
rz(-1.3432608) q[0];
x q[1];
rz(1.6283054) q[2];
sx q[2];
rz(-1.3293666) q[2];
sx q[2];
rz(1.4032569) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.64046265) q[1];
sx q[1];
rz(-1.340133) q[1];
sx q[1];
rz(-0.035236852) q[1];
rz(-pi) q[2];
rz(0.2749625) q[3];
sx q[3];
rz(-2.5841613) q[3];
sx q[3];
rz(2.9334588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9589299) q[2];
sx q[2];
rz(-1.485433) q[2];
sx q[2];
rz(0.043969285) q[2];
rz(0.26642695) q[3];
sx q[3];
rz(-3.0156342) q[3];
sx q[3];
rz(0.21336475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5037395) q[0];
sx q[0];
rz(-1.6197236) q[0];
sx q[0];
rz(-0.36852401) q[0];
rz(0.35399327) q[1];
sx q[1];
rz(-2.0123672) q[1];
sx q[1];
rz(-1.0634208) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3388004) q[0];
sx q[0];
rz(-1.430585) q[0];
sx q[0];
rz(-3.1288052) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2930128) q[2];
sx q[2];
rz(-1.4627224) q[2];
sx q[2];
rz(1.9982823) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.271558) q[1];
sx q[1];
rz(-1.1599301) q[1];
sx q[1];
rz(-2.690587) q[1];
rz(-2.9408089) q[3];
sx q[3];
rz(-2.3827083) q[3];
sx q[3];
rz(-2.6155265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4949558) q[2];
sx q[2];
rz(-1.8821303) q[2];
sx q[2];
rz(1.7156401) q[2];
rz(0.87092733) q[3];
sx q[3];
rz(-1.0692853) q[3];
sx q[3];
rz(-0.23437414) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17212412) q[0];
sx q[0];
rz(-0.5640465) q[0];
sx q[0];
rz(3.1308351) q[0];
rz(-2.5116008) q[1];
sx q[1];
rz(-1.0456345) q[1];
sx q[1];
rz(-2.2363372) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0765823) q[0];
sx q[0];
rz(-0.53758756) q[0];
sx q[0];
rz(1.697886) q[0];
x q[1];
rz(1.0192746) q[2];
sx q[2];
rz(-2.5576303) q[2];
sx q[2];
rz(1.1332741) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.39054444) q[1];
sx q[1];
rz(-2.6169951) q[1];
sx q[1];
rz(2.8538997) q[1];
rz(-0.91628381) q[3];
sx q[3];
rz(-1.223151) q[3];
sx q[3];
rz(-0.8485444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4524727) q[2];
sx q[2];
rz(-1.8931188) q[2];
sx q[2];
rz(0.11736891) q[2];
rz(-2.4145685) q[3];
sx q[3];
rz(-0.89594642) q[3];
sx q[3];
rz(-2.0015543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3177976) q[0];
sx q[0];
rz(-1.057484) q[0];
sx q[0];
rz(0.29528883) q[0];
rz(-1.7501887) q[1];
sx q[1];
rz(-2.4877805) q[1];
sx q[1];
rz(2.6546435) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21509296) q[0];
sx q[0];
rz(-1.7679259) q[0];
sx q[0];
rz(-1.4678982) q[0];
rz(2.0980663) q[2];
sx q[2];
rz(-2.5621116) q[2];
sx q[2];
rz(-2.5763489) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0978969) q[1];
sx q[1];
rz(-2.0123031) q[1];
sx q[1];
rz(-0.82951464) q[1];
x q[2];
rz(-0.26712931) q[3];
sx q[3];
rz(-1.7448269) q[3];
sx q[3];
rz(1.4166299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9140909) q[2];
sx q[2];
rz(-2.5469683) q[2];
sx q[2];
rz(1.5957069) q[2];
rz(-2.5325736) q[3];
sx q[3];
rz(-2.0273429) q[3];
sx q[3];
rz(2.473089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.2862227) q[0];
sx q[0];
rz(-2.2384221) q[0];
sx q[0];
rz(-1.632593) q[0];
rz(1.1434932) q[1];
sx q[1];
rz(-1.9805464) q[1];
sx q[1];
rz(0.11018363) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78980434) q[0];
sx q[0];
rz(-2.1579701) q[0];
sx q[0];
rz(-1.6950865) q[0];
rz(-pi) q[1];
rz(1.1082591) q[2];
sx q[2];
rz(-1.8018844) q[2];
sx q[2];
rz(-0.05751911) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5852803) q[1];
sx q[1];
rz(-1.9686799) q[1];
sx q[1];
rz(0.75979397) q[1];
x q[2];
rz(1.5694782) q[3];
sx q[3];
rz(-1.073488) q[3];
sx q[3];
rz(0.2494825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1630359) q[2];
sx q[2];
rz(-1.0806934) q[2];
sx q[2];
rz(-0.85868305) q[2];
rz(-1.5052172) q[3];
sx q[3];
rz(-1.3450164) q[3];
sx q[3];
rz(1.6987364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.7929512) q[0];
sx q[0];
rz(-2.0703147) q[0];
sx q[0];
rz(0.68565482) q[0];
rz(0.81949702) q[1];
sx q[1];
rz(-1.6623431) q[1];
sx q[1];
rz(2.562838) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98812719) q[0];
sx q[0];
rz(-2.1442226) q[0];
sx q[0];
rz(0.98116409) q[0];
x q[1];
rz(-2.1376325) q[2];
sx q[2];
rz(-1.3610148) q[2];
sx q[2];
rz(-0.45406859) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.600665) q[1];
sx q[1];
rz(-1.828421) q[1];
sx q[1];
rz(-0.39389807) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41865809) q[3];
sx q[3];
rz(-0.11944977) q[3];
sx q[3];
rz(-0.87629139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3026352) q[2];
sx q[2];
rz(-0.18582782) q[2];
sx q[2];
rz(2.9616984) q[2];
rz(-0.57080615) q[3];
sx q[3];
rz(-2.4195318) q[3];
sx q[3];
rz(-0.65176898) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8410692) q[0];
sx q[0];
rz(-2.4916334) q[0];
sx q[0];
rz(-1.4148225) q[0];
rz(-2.6534446) q[1];
sx q[1];
rz(-2.585325) q[1];
sx q[1];
rz(-1.5811031) q[1];
rz(-2.2019501) q[2];
sx q[2];
rz(-1.2705391) q[2];
sx q[2];
rz(1.4654797) q[2];
rz(1.855424) q[3];
sx q[3];
rz(-1.0502001) q[3];
sx q[3];
rz(-2.735103) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
