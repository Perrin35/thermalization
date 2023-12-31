OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.055846) q[0];
sx q[0];
rz(-3.0598109) q[0];
sx q[0];
rz(2.6401289) q[0];
rz(1.4986829) q[1];
sx q[1];
rz(-2.745435) q[1];
sx q[1];
rz(-0.3224386) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0599521) q[0];
sx q[0];
rz(-1.8917221) q[0];
sx q[0];
rz(-3.0517464) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6104923) q[2];
sx q[2];
rz(-2.1571026) q[2];
sx q[2];
rz(0.58639484) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8391708) q[1];
sx q[1];
rz(-1.9170554) q[1];
sx q[1];
rz(-2.195921) q[1];
rz(-pi) q[2];
x q[2];
rz(2.206771) q[3];
sx q[3];
rz(-1.1879731) q[3];
sx q[3];
rz(-0.32360199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.50513187) q[2];
sx q[2];
rz(-0.59288609) q[2];
sx q[2];
rz(2.5855529) q[2];
rz(2.3089144) q[3];
sx q[3];
rz(-1.4913538) q[3];
sx q[3];
rz(-2.1957943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44822025) q[0];
sx q[0];
rz(-1.4602666) q[0];
sx q[0];
rz(-0.15727501) q[0];
rz(2.8804624) q[1];
sx q[1];
rz(-1.3477247) q[1];
sx q[1];
rz(3.0325586) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.581668) q[0];
sx q[0];
rz(-1.7413057) q[0];
sx q[0];
rz(-1.8018363) q[0];
x q[1];
rz(-2.5017545) q[2];
sx q[2];
rz(-0.32391641) q[2];
sx q[2];
rz(1.4878291) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2687159) q[1];
sx q[1];
rz(-1.1040338) q[1];
sx q[1];
rz(-2.4750535) q[1];
rz(-1.3974959) q[3];
sx q[3];
rz(-1.0565851) q[3];
sx q[3];
rz(1.7931995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9033501) q[2];
sx q[2];
rz(-1.976333) q[2];
sx q[2];
rz(-1.8781352) q[2];
rz(2.8144828) q[3];
sx q[3];
rz(-1.5771022) q[3];
sx q[3];
rz(1.9272778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064421244) q[0];
sx q[0];
rz(-0.049296878) q[0];
sx q[0];
rz(-1.7984614) q[0];
rz(0.24761565) q[1];
sx q[1];
rz(-2.394948) q[1];
sx q[1];
rz(0.48167357) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1789078) q[0];
sx q[0];
rz(-0.6324397) q[0];
sx q[0];
rz(-0.90895598) q[0];
rz(2.0273676) q[2];
sx q[2];
rz(-2.4764428) q[2];
sx q[2];
rz(0.84012023) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3232705) q[1];
sx q[1];
rz(-0.96410492) q[1];
sx q[1];
rz(0.89241772) q[1];
x q[2];
rz(-0.62044575) q[3];
sx q[3];
rz(-2.1944322) q[3];
sx q[3];
rz(0.40943957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8032916) q[2];
sx q[2];
rz(-0.81739601) q[2];
sx q[2];
rz(2.6417007) q[2];
rz(2.5806184) q[3];
sx q[3];
rz(-1.8818972) q[3];
sx q[3];
rz(-1.6803754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19701476) q[0];
sx q[0];
rz(-2.975583) q[0];
sx q[0];
rz(-2.5894077) q[0];
rz(-1.588297) q[1];
sx q[1];
rz(-0.89893666) q[1];
sx q[1];
rz(1.2447371) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8053999) q[0];
sx q[0];
rz(-1.5691783) q[0];
sx q[0];
rz(-2.8021115) q[0];
x q[1];
rz(-1.7719901) q[2];
sx q[2];
rz(-1.6987213) q[2];
sx q[2];
rz(-0.97845562) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1340449) q[1];
sx q[1];
rz(-1.3589077) q[1];
sx q[1];
rz(2.120554) q[1];
x q[2];
rz(0.34282116) q[3];
sx q[3];
rz(-2.5431513) q[3];
sx q[3];
rz(-3.0330021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.84919471) q[2];
sx q[2];
rz(-1.2595824) q[2];
sx q[2];
rz(-1.9909031) q[2];
rz(-1.6644647) q[3];
sx q[3];
rz(-1.5090347) q[3];
sx q[3];
rz(0.48294827) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078159049) q[0];
sx q[0];
rz(-0.76197356) q[0];
sx q[0];
rz(0.081469014) q[0];
rz(-0.062462656) q[1];
sx q[1];
rz(-2.000258) q[1];
sx q[1];
rz(1.6385471) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016184729) q[0];
sx q[0];
rz(-0.59690079) q[0];
sx q[0];
rz(1.7993268) q[0];
x q[1];
rz(-0.91707768) q[2];
sx q[2];
rz(-3.0055025) q[2];
sx q[2];
rz(1.2166785) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7170143) q[1];
sx q[1];
rz(-1.6728405) q[1];
sx q[1];
rz(-1.8841519) q[1];
rz(1.2522069) q[3];
sx q[3];
rz(-1.0831523) q[3];
sx q[3];
rz(-2.3209751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.5082671) q[2];
sx q[2];
rz(-0.94255629) q[2];
sx q[2];
rz(1.903669) q[2];
rz(-1.1226908) q[3];
sx q[3];
rz(-0.676238) q[3];
sx q[3];
rz(-0.52156633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6102819) q[0];
sx q[0];
rz(-2.1370482) q[0];
sx q[0];
rz(2.8748728) q[0];
rz(2.5807014) q[1];
sx q[1];
rz(-1.8436878) q[1];
sx q[1];
rz(2.3430603) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.079433867) q[0];
sx q[0];
rz(-1.845813) q[0];
sx q[0];
rz(-2.9955203) q[0];
rz(-0.49662923) q[2];
sx q[2];
rz(-1.8804669) q[2];
sx q[2];
rz(1.3336099) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8171463) q[1];
sx q[1];
rz(-1.6231316) q[1];
sx q[1];
rz(0.99919341) q[1];
rz(-0.9540326) q[3];
sx q[3];
rz(-2.2700689) q[3];
sx q[3];
rz(1.9642533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9810527) q[2];
sx q[2];
rz(-1.2489698) q[2];
sx q[2];
rz(2.7748761) q[2];
rz(1.8803053) q[3];
sx q[3];
rz(-1.4533318) q[3];
sx q[3];
rz(0.096207531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7325608) q[0];
sx q[0];
rz(-0.92027396) q[0];
sx q[0];
rz(2.5352056) q[0];
rz(-2.9442893) q[1];
sx q[1];
rz(-1.1261255) q[1];
sx q[1];
rz(2.6775449) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7921917) q[0];
sx q[0];
rz(-2.1878625) q[0];
sx q[0];
rz(-2.6130136) q[0];
rz(-pi) q[1];
rz(1.8115933) q[2];
sx q[2];
rz(-1.5511998) q[2];
sx q[2];
rz(2.7407676) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6021767) q[1];
sx q[1];
rz(-0.5836986) q[1];
sx q[1];
rz(-1.0023487) q[1];
rz(1.6528321) q[3];
sx q[3];
rz(-1.335252) q[3];
sx q[3];
rz(-1.9627278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.93418926) q[2];
sx q[2];
rz(-2.1384017) q[2];
sx q[2];
rz(-0.25804538) q[2];
rz(1.9559654) q[3];
sx q[3];
rz(-1.6262755) q[3];
sx q[3];
rz(-0.0035088249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.946452) q[0];
sx q[0];
rz(-1.2807245) q[0];
sx q[0];
rz(2.7602957) q[0];
rz(-0.095245846) q[1];
sx q[1];
rz(-2.1691599) q[1];
sx q[1];
rz(-1.7000748) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1822752) q[0];
sx q[0];
rz(-1.5858985) q[0];
sx q[0];
rz(0.063477593) q[0];
rz(-pi) q[1];
rz(1.6091789) q[2];
sx q[2];
rz(-0.95853171) q[2];
sx q[2];
rz(-1.6910764) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6724832) q[1];
sx q[1];
rz(-1.6961349) q[1];
sx q[1];
rz(-0.13087665) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6137245) q[3];
sx q[3];
rz(-1.2371484) q[3];
sx q[3];
rz(-1.7121401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9986481) q[2];
sx q[2];
rz(-2.7286077) q[2];
sx q[2];
rz(-2.9150035) q[2];
rz(0.44858027) q[3];
sx q[3];
rz(-1.5357163) q[3];
sx q[3];
rz(-0.8297689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7609693) q[0];
sx q[0];
rz(-0.7614823) q[0];
sx q[0];
rz(1.3990078) q[0];
rz(-0.31708583) q[1];
sx q[1];
rz(-1.4750907) q[1];
sx q[1];
rz(-2.1549966) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.164455) q[0];
sx q[0];
rz(-1.7286321) q[0];
sx q[0];
rz(0.49789238) q[0];
rz(-0.95790205) q[2];
sx q[2];
rz(-2.3447678) q[2];
sx q[2];
rz(-0.56607841) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1074859) q[1];
sx q[1];
rz(-1.7120546) q[1];
sx q[1];
rz(-2.5320401) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1986198) q[3];
sx q[3];
rz(-0.44789207) q[3];
sx q[3];
rz(-0.4263634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9265147) q[2];
sx q[2];
rz(-0.72700095) q[2];
sx q[2];
rz(-0.40965664) q[2];
rz(0.26327291) q[3];
sx q[3];
rz(-1.3116838) q[3];
sx q[3];
rz(2.6221361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(0.39919329) q[0];
sx q[0];
rz(-0.078646794) q[0];
sx q[0];
rz(-1.4051399) q[0];
rz(-0.82110226) q[1];
sx q[1];
rz(-0.91870538) q[1];
sx q[1];
rz(1.4155037) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9916519) q[0];
sx q[0];
rz(-1.2399925) q[0];
sx q[0];
rz(-1.3209692) q[0];
rz(2.892898) q[2];
sx q[2];
rz(-1.2173614) q[2];
sx q[2];
rz(0.24662185) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6092097) q[1];
sx q[1];
rz(-1.4192974) q[1];
sx q[1];
rz(-1.7864368) q[1];
rz(-pi) q[2];
rz(2.0831765) q[3];
sx q[3];
rz(-1.4231764) q[3];
sx q[3];
rz(0.83422134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4225509) q[2];
sx q[2];
rz(-2.8406403) q[2];
sx q[2];
rz(-0.12410513) q[2];
rz(0.96578807) q[3];
sx q[3];
rz(-1.4744759) q[3];
sx q[3];
rz(-0.66108274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60349764) q[0];
sx q[0];
rz(-2.8932543) q[0];
sx q[0];
rz(2.2809991) q[0];
rz(2.8339236) q[1];
sx q[1];
rz(-1.2528906) q[1];
sx q[1];
rz(1.2045592) q[1];
rz(2.5406191) q[2];
sx q[2];
rz(-2.8291611) q[2];
sx q[2];
rz(-0.57401382) q[2];
rz(-1.7210759) q[3];
sx q[3];
rz(-2.1223304) q[3];
sx q[3];
rz(0.35029678) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
