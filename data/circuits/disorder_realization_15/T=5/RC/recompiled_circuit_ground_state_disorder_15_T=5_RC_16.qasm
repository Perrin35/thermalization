OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6599967) q[0];
sx q[0];
rz(-0.41756088) q[0];
sx q[0];
rz(-0.041393809) q[0];
rz(0.59250915) q[1];
sx q[1];
rz(-0.23696391) q[1];
sx q[1];
rz(-0.89495975) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1101018) q[0];
sx q[0];
rz(-0.43927017) q[0];
sx q[0];
rz(1.8861559) q[0];
x q[1];
rz(-2.824895) q[2];
sx q[2];
rz(-2.0513655) q[2];
sx q[2];
rz(0.66083913) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0022566) q[1];
sx q[1];
rz(-1.8426241) q[1];
sx q[1];
rz(-0.28693954) q[1];
rz(-1.5271321) q[3];
sx q[3];
rz(-1.6796675) q[3];
sx q[3];
rz(-3.0565232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71901739) q[2];
sx q[2];
rz(-0.97192478) q[2];
sx q[2];
rz(1.0389339) q[2];
rz(0.59764189) q[3];
sx q[3];
rz(-2.938275) q[3];
sx q[3];
rz(1.982127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0036302) q[0];
sx q[0];
rz(-0.85145402) q[0];
sx q[0];
rz(-2.8077069) q[0];
rz(2.4489898) q[1];
sx q[1];
rz(-1.2665766) q[1];
sx q[1];
rz(2.0803221) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6491566) q[0];
sx q[0];
rz(-1.9266495) q[0];
sx q[0];
rz(-2.27764) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40410903) q[2];
sx q[2];
rz(-1.3349018) q[2];
sx q[2];
rz(2.1020232) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1352219) q[1];
sx q[1];
rz(-1.3961482) q[1];
sx q[1];
rz(2.6937301) q[1];
rz(-pi) q[2];
rz(-2.1718086) q[3];
sx q[3];
rz(-2.8967987) q[3];
sx q[3];
rz(-0.56772619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4873203) q[2];
sx q[2];
rz(-1.6795748) q[2];
sx q[2];
rz(0.1839323) q[2];
rz(-3.1390624) q[3];
sx q[3];
rz(-0.80357426) q[3];
sx q[3];
rz(-2.7185503) q[3];
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
rz(2.2976487) q[0];
sx q[0];
rz(-2.9941445) q[0];
sx q[0];
rz(2.3609128) q[0];
rz(2.2536904) q[1];
sx q[1];
rz(-1.0941411) q[1];
sx q[1];
rz(-2.323774) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2116341) q[0];
sx q[0];
rz(-0.51211905) q[0];
sx q[0];
rz(-1.0996487) q[0];
x q[1];
rz(0.79105241) q[2];
sx q[2];
rz(-1.0246236) q[2];
sx q[2];
rz(2.5763022) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1314268) q[1];
sx q[1];
rz(-1.5714116) q[1];
sx q[1];
rz(0.58227957) q[1];
x q[2];
rz(-2.6894916) q[3];
sx q[3];
rz(-1.107405) q[3];
sx q[3];
rz(-0.85833651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4320977) q[2];
sx q[2];
rz(-0.12505394) q[2];
sx q[2];
rz(-1.032426) q[2];
rz(0.89140511) q[3];
sx q[3];
rz(-2.4668507) q[3];
sx q[3];
rz(2.3064822) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29435232) q[0];
sx q[0];
rz(-2.8645741) q[0];
sx q[0];
rz(-2.5936122) q[0];
rz(0.056593865) q[1];
sx q[1];
rz(-1.2013925) q[1];
sx q[1];
rz(-0.024554575) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0191607) q[0];
sx q[0];
rz(-1.1709899) q[0];
sx q[0];
rz(-2.7960577) q[0];
rz(-pi) q[1];
rz(-1.9664832) q[2];
sx q[2];
rz(-1.2597689) q[2];
sx q[2];
rz(2.0360586) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7112507) q[1];
sx q[1];
rz(-0.81626662) q[1];
sx q[1];
rz(-1.0984983) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6484429) q[3];
sx q[3];
rz(-0.34420612) q[3];
sx q[3];
rz(1.5045117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5121439) q[2];
sx q[2];
rz(-0.63261837) q[2];
sx q[2];
rz(2.4150685) q[2];
rz(1.7065382) q[3];
sx q[3];
rz(-1.8972998) q[3];
sx q[3];
rz(-1.4783036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59359819) q[0];
sx q[0];
rz(-1.7578121) q[0];
sx q[0];
rz(-1.7543678) q[0];
rz(-0.45135003) q[1];
sx q[1];
rz(-0.74140048) q[1];
sx q[1];
rz(-0.50419921) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057495385) q[0];
sx q[0];
rz(-1.7404489) q[0];
sx q[0];
rz(-1.6952362) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69602698) q[2];
sx q[2];
rz(-2.2559204) q[2];
sx q[2];
rz(2.7643577) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5686149) q[1];
sx q[1];
rz(-1.0573799) q[1];
sx q[1];
rz(-0.77691369) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80662722) q[3];
sx q[3];
rz(-0.99116814) q[3];
sx q[3];
rz(1.5999853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35895687) q[2];
sx q[2];
rz(-2.2613596) q[2];
sx q[2];
rz(3.1202988) q[2];
rz(3.0535789) q[3];
sx q[3];
rz(-0.26282495) q[3];
sx q[3];
rz(-0.16784856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.39040318) q[0];
sx q[0];
rz(-2.5683371) q[0];
sx q[0];
rz(0.34725749) q[0];
rz(-1.4561397) q[1];
sx q[1];
rz(-0.29734722) q[1];
sx q[1];
rz(0.28486326) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69377855) q[0];
sx q[0];
rz(-1.283857) q[0];
sx q[0];
rz(1.2579495) q[0];
rz(-pi) q[1];
rz(-0.059528298) q[2];
sx q[2];
rz(-0.7452508) q[2];
sx q[2];
rz(-1.1203114) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.59249944) q[1];
sx q[1];
rz(-0.91895267) q[1];
sx q[1];
rz(-2.8054906) q[1];
rz(-2.8480356) q[3];
sx q[3];
rz(-1.2282853) q[3];
sx q[3];
rz(2.5409043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.138729) q[2];
sx q[2];
rz(-1.8314654) q[2];
sx q[2];
rz(-2.1076473) q[2];
rz(-2.406534) q[3];
sx q[3];
rz(-1.4933973) q[3];
sx q[3];
rz(-2.471931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1351778) q[0];
sx q[0];
rz(-2.4786351) q[0];
sx q[0];
rz(-2.2354777) q[0];
rz(-0.16618973) q[1];
sx q[1];
rz(-2.6527185) q[1];
sx q[1];
rz(-0.30950549) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0565694) q[0];
sx q[0];
rz(-2.8583953) q[0];
sx q[0];
rz(-0.92206456) q[0];
rz(0.1997582) q[2];
sx q[2];
rz(-2.3383378) q[2];
sx q[2];
rz(1.8920395) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2848572) q[1];
sx q[1];
rz(-1.3580163) q[1];
sx q[1];
rz(2.2163843) q[1];
x q[2];
rz(1.4847912) q[3];
sx q[3];
rz(-1.5218867) q[3];
sx q[3];
rz(2.870056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6918148) q[2];
sx q[2];
rz(-0.66902995) q[2];
sx q[2];
rz(1.6548033) q[2];
rz(-2.6438223) q[3];
sx q[3];
rz(-2.3064752) q[3];
sx q[3];
rz(2.914079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19787702) q[0];
sx q[0];
rz(-0.72138041) q[0];
sx q[0];
rz(2.873514) q[0];
rz(-0.90191853) q[1];
sx q[1];
rz(-2.0733158) q[1];
sx q[1];
rz(1.8394151) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81277455) q[0];
sx q[0];
rz(-0.98085603) q[0];
sx q[0];
rz(-1.7699715) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5728358) q[2];
sx q[2];
rz(-2.0173948) q[2];
sx q[2];
rz(-2.0799321) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8852325) q[1];
sx q[1];
rz(-2.5361885) q[1];
sx q[1];
rz(-1.3219236) q[1];
x q[2];
rz(-1.2648495) q[3];
sx q[3];
rz(-1.5380203) q[3];
sx q[3];
rz(2.2900801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7997718) q[2];
sx q[2];
rz(-1.9016966) q[2];
sx q[2];
rz(2.5566901) q[2];
rz(1.2029485) q[3];
sx q[3];
rz(-1.3936309) q[3];
sx q[3];
rz(-1.5486307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(1.0840266) q[0];
sx q[0];
rz(-2.6507222) q[0];
sx q[0];
rz(0.67114818) q[0];
rz(0.28823832) q[1];
sx q[1];
rz(-0.75575525) q[1];
sx q[1];
rz(-0.63405687) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6995981) q[0];
sx q[0];
rz(-1.5324679) q[0];
sx q[0];
rz(-3.1257191) q[0];
rz(-pi) q[1];
rz(3.0915758) q[2];
sx q[2];
rz(-1.7361618) q[2];
sx q[2];
rz(0.53420137) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.94407627) q[1];
sx q[1];
rz(-1.895525) q[1];
sx q[1];
rz(-1.5438117) q[1];
rz(-pi) q[2];
rz(-2.4692772) q[3];
sx q[3];
rz(-1.3642715) q[3];
sx q[3];
rz(-0.12609161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5481698) q[2];
sx q[2];
rz(-1.5680743) q[2];
sx q[2];
rz(2.8578952) q[2];
rz(-0.13188322) q[3];
sx q[3];
rz(-0.3564842) q[3];
sx q[3];
rz(2.5186445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.236096) q[0];
sx q[0];
rz(-0.14362366) q[0];
sx q[0];
rz(2.1719601) q[0];
rz(-0.49599221) q[1];
sx q[1];
rz(-2.453936) q[1];
sx q[1];
rz(-2.2775441) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0522033) q[0];
sx q[0];
rz(-1.6411575) q[0];
sx q[0];
rz(-2.2785827) q[0];
x q[1];
rz(2.8882083) q[2];
sx q[2];
rz(-1.9929427) q[2];
sx q[2];
rz(2.6146023) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7590027) q[1];
sx q[1];
rz(-0.99901268) q[1];
sx q[1];
rz(-0.60061213) q[1];
rz(-2.2275023) q[3];
sx q[3];
rz(-1.2063081) q[3];
sx q[3];
rz(-1.2165704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.25829092) q[2];
sx q[2];
rz(-0.51670462) q[2];
sx q[2];
rz(1.0374163) q[2];
rz(-2.9344905) q[3];
sx q[3];
rz(-2.243302) q[3];
sx q[3];
rz(-2.5521539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0265738) q[0];
sx q[0];
rz(-1.4587695) q[0];
sx q[0];
rz(1.2039626) q[0];
rz(-0.013982458) q[1];
sx q[1];
rz(-1.7054066) q[1];
sx q[1];
rz(-1.1048497) q[1];
rz(-1.3399765) q[2];
sx q[2];
rz(-1.6047819) q[2];
sx q[2];
rz(-2.6699382) q[2];
rz(-2.8057475) q[3];
sx q[3];
rz(-2.2686979) q[3];
sx q[3];
rz(3.1115193) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
