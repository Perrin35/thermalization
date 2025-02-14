OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4021969) q[0];
sx q[0];
rz(-0.88391179) q[0];
sx q[0];
rz(-2.28595) q[0];
rz(2.9698676) q[1];
sx q[1];
rz(-3.0260234) q[1];
sx q[1];
rz(0.56245437) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5006258) q[0];
sx q[0];
rz(-1.5175227) q[0];
sx q[0];
rz(-1.9356273) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0777656) q[2];
sx q[2];
rz(-1.4344782) q[2];
sx q[2];
rz(0.16198128) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2515638) q[1];
sx q[1];
rz(-1.9891796) q[1];
sx q[1];
rz(-1.1226001) q[1];
x q[2];
rz(0.39333087) q[3];
sx q[3];
rz(-1.9489904) q[3];
sx q[3];
rz(2.8761169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.90613753) q[2];
sx q[2];
rz(-1.9635341) q[2];
sx q[2];
rz(0.79929024) q[2];
rz(0.47131395) q[3];
sx q[3];
rz(-2.2282232) q[3];
sx q[3];
rz(-0.93588626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1021295) q[0];
sx q[0];
rz(-1.8164182) q[0];
sx q[0];
rz(1.7934196) q[0];
rz(1.3735636) q[1];
sx q[1];
rz(-1.1496239) q[1];
sx q[1];
rz(-1.1522393) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9236167) q[0];
sx q[0];
rz(-2.2647175) q[0];
sx q[0];
rz(0.38461916) q[0];
rz(-2.197108) q[2];
sx q[2];
rz(-0.99530333) q[2];
sx q[2];
rz(1.6863509) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.79346953) q[1];
sx q[1];
rz(-1.4033699) q[1];
sx q[1];
rz(2.4504285) q[1];
rz(-pi) q[2];
rz(-2.1227319) q[3];
sx q[3];
rz(-1.1180919) q[3];
sx q[3];
rz(-0.34358968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3422602) q[2];
sx q[2];
rz(-0.41877425) q[2];
sx q[2];
rz(-2.2104134) q[2];
rz(3.0350507) q[3];
sx q[3];
rz(-2.0276766) q[3];
sx q[3];
rz(0.60025269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0837285) q[0];
sx q[0];
rz(-1.3902384) q[0];
sx q[0];
rz(2.6237543) q[0];
rz(-2.2528516) q[1];
sx q[1];
rz(-0.70777142) q[1];
sx q[1];
rz(-2.6944366) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06181051) q[0];
sx q[0];
rz(-1.634565) q[0];
sx q[0];
rz(0.26173862) q[0];
x q[1];
rz(-2.4141623) q[2];
sx q[2];
rz(-2.0038249) q[2];
sx q[2];
rz(-2.7114781) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4616926) q[1];
sx q[1];
rz(-0.60631207) q[1];
sx q[1];
rz(1.8099422) q[1];
x q[2];
rz(-2.6833862) q[3];
sx q[3];
rz(-0.57469207) q[3];
sx q[3];
rz(1.9783057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37457028) q[2];
sx q[2];
rz(-0.76408237) q[2];
sx q[2];
rz(-0.53263295) q[2];
rz(2.8940708) q[3];
sx q[3];
rz(-0.73740021) q[3];
sx q[3];
rz(-1.0386764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6240876) q[0];
sx q[0];
rz(-3.0563323) q[0];
sx q[0];
rz(0.10661539) q[0];
rz(0.34635776) q[1];
sx q[1];
rz(-0.84500161) q[1];
sx q[1];
rz(1.1603629) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24460565) q[0];
sx q[0];
rz(-1.7443313) q[0];
sx q[0];
rz(0.24994295) q[0];
x q[1];
rz(-2.423942) q[2];
sx q[2];
rz(-2.1776878) q[2];
sx q[2];
rz(-0.90536149) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44768279) q[1];
sx q[1];
rz(-0.47397787) q[1];
sx q[1];
rz(-0.86556566) q[1];
rz(1.0456234) q[3];
sx q[3];
rz(-1.1743288) q[3];
sx q[3];
rz(-0.1016271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0059263) q[2];
sx q[2];
rz(-0.29834193) q[2];
sx q[2];
rz(3.0653817) q[2];
rz(-0.58586079) q[3];
sx q[3];
rz(-1.1548837) q[3];
sx q[3];
rz(-1.3425672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12014408) q[0];
sx q[0];
rz(-1.3055389) q[0];
sx q[0];
rz(0.35010499) q[0];
rz(2.1856951) q[1];
sx q[1];
rz(-1.8616385) q[1];
sx q[1];
rz(1.9116481) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1317539) q[0];
sx q[0];
rz(-1.3309877) q[0];
sx q[0];
rz(1.9480223) q[0];
x q[1];
rz(1.7857871) q[2];
sx q[2];
rz(-0.16915288) q[2];
sx q[2];
rz(0.40119888) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.753641) q[1];
sx q[1];
rz(-1.759583) q[1];
sx q[1];
rz(-1.9061879) q[1];
x q[2];
rz(0.1996207) q[3];
sx q[3];
rz(-0.99310447) q[3];
sx q[3];
rz(2.9132089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2436287) q[2];
sx q[2];
rz(-2.882759) q[2];
sx q[2];
rz(-1.492307) q[2];
rz(2.7367075) q[3];
sx q[3];
rz(-0.76571524) q[3];
sx q[3];
rz(2.305472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.006007) q[0];
sx q[0];
rz(-2.0642991) q[0];
sx q[0];
rz(-2.5591922) q[0];
rz(0.68663418) q[1];
sx q[1];
rz(-1.735894) q[1];
sx q[1];
rz(-1.883421) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5053537) q[0];
sx q[0];
rz(-2.6864144) q[0];
sx q[0];
rz(1.8591465) q[0];
rz(-pi) q[1];
rz(1.011884) q[2];
sx q[2];
rz(-2.5311573) q[2];
sx q[2];
rz(-1.287078) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4542089) q[1];
sx q[1];
rz(-1.5011678) q[1];
sx q[1];
rz(1.5147665) q[1];
rz(-pi) q[2];
rz(1.2518796) q[3];
sx q[3];
rz(-0.61781672) q[3];
sx q[3];
rz(1.2417254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.76576343) q[2];
sx q[2];
rz(-1.148369) q[2];
sx q[2];
rz(1.0180391) q[2];
rz(2.8954519) q[3];
sx q[3];
rz(-1.7459511) q[3];
sx q[3];
rz(1.5181946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(0.70596424) q[0];
sx q[0];
rz(-2.3376597) q[0];
sx q[0];
rz(0.58018082) q[0];
rz(-2.9980581) q[1];
sx q[1];
rz(-0.48964557) q[1];
sx q[1];
rz(-2.9339583) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7589446) q[0];
sx q[0];
rz(-0.84340723) q[0];
sx q[0];
rz(-1.101053) q[0];
rz(-2.1499499) q[2];
sx q[2];
rz(-1.8657547) q[2];
sx q[2];
rz(-2.7911012) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.27397284) q[1];
sx q[1];
rz(-1.311353) q[1];
sx q[1];
rz(-0.8670437) q[1];
x q[2];
rz(-1.4752046) q[3];
sx q[3];
rz(-1.2089653) q[3];
sx q[3];
rz(-1.3528479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2988854) q[2];
sx q[2];
rz(-1.2879813) q[2];
sx q[2];
rz(2.5353954) q[2];
rz(-2.4541564) q[3];
sx q[3];
rz(-1.1339374) q[3];
sx q[3];
rz(0.70639759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6523478) q[0];
sx q[0];
rz(-0.027996538) q[0];
sx q[0];
rz(1.0580753) q[0];
rz(3.0319013) q[1];
sx q[1];
rz(-1.1166162) q[1];
sx q[1];
rz(-1.6995957) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083695166) q[0];
sx q[0];
rz(-2.6287615) q[0];
sx q[0];
rz(-0.67016853) q[0];
rz(-pi) q[1];
rz(-1.2953561) q[2];
sx q[2];
rz(-2.0552962) q[2];
sx q[2];
rz(0.4415919) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.31046384) q[1];
sx q[1];
rz(-1.0894686) q[1];
sx q[1];
rz(2.3075065) q[1];
x q[2];
rz(-1.8423002) q[3];
sx q[3];
rz(-1.3437004) q[3];
sx q[3];
rz(0.61448594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.480392) q[2];
sx q[2];
rz(-2.9471687) q[2];
sx q[2];
rz(-1.9332168) q[2];
rz(-0.66323534) q[3];
sx q[3];
rz(-1.6845208) q[3];
sx q[3];
rz(1.9251582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1709568) q[0];
sx q[0];
rz(-0.96681505) q[0];
sx q[0];
rz(3.048625) q[0];
rz(1.8611106) q[1];
sx q[1];
rz(-2.4130776) q[1];
sx q[1];
rz(0.13519898) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8672392) q[0];
sx q[0];
rz(-1.5177736) q[0];
sx q[0];
rz(-1.5252602) q[0];
x q[1];
rz(-0.5298631) q[2];
sx q[2];
rz(-2.8447897) q[2];
sx q[2];
rz(2.5713845) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1108023) q[1];
sx q[1];
rz(-0.25486481) q[1];
sx q[1];
rz(0.084666208) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0207094) q[3];
sx q[3];
rz(-2.9723047) q[3];
sx q[3];
rz(-1.5811046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.70904237) q[2];
sx q[2];
rz(-1.21864) q[2];
sx q[2];
rz(-0.77152983) q[2];
rz(0.49154526) q[3];
sx q[3];
rz(-1.8621657) q[3];
sx q[3];
rz(0.5947203) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5563357) q[0];
sx q[0];
rz(-2.519643) q[0];
sx q[0];
rz(1.1248032) q[0];
rz(-1.8428165) q[1];
sx q[1];
rz(-0.61538428) q[1];
sx q[1];
rz(-0.76464701) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0394093) q[0];
sx q[0];
rz(-1.656807) q[0];
sx q[0];
rz(-0.10552222) q[0];
x q[1];
rz(2.6420399) q[2];
sx q[2];
rz(-0.95328125) q[2];
sx q[2];
rz(-0.90047405) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5950039) q[1];
sx q[1];
rz(-1.6975132) q[1];
sx q[1];
rz(-1.3664043) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1922791) q[3];
sx q[3];
rz(-2.0928185) q[3];
sx q[3];
rz(1.6917563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7654968) q[2];
sx q[2];
rz(-1.2906047) q[2];
sx q[2];
rz(0.20720227) q[2];
rz(-0.97992212) q[3];
sx q[3];
rz(-1.1332952) q[3];
sx q[3];
rz(2.9492212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1205263) q[0];
sx q[0];
rz(-1.4561894) q[0];
sx q[0];
rz(-0.86984632) q[0];
rz(2.6407241) q[1];
sx q[1];
rz(-0.23575467) q[1];
sx q[1];
rz(-2.2101319) q[1];
rz(1.6899213) q[2];
sx q[2];
rz(-1.7075734) q[2];
sx q[2];
rz(1.3592958) q[2];
rz(1.408314) q[3];
sx q[3];
rz(-1.3357031) q[3];
sx q[3];
rz(1.6817844) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
