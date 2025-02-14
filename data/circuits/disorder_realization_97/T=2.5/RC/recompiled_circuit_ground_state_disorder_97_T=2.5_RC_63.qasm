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
rz(3.1273841) q[0];
sx q[0];
rz(9.4655884) q[0];
rz(-0.040925097) q[1];
sx q[1];
rz(-1.0558145) q[1];
sx q[1];
rz(-0.55650401) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4427821) q[0];
sx q[0];
rz(-0.39467803) q[0];
sx q[0];
rz(1.077432) q[0];
rz(-pi) q[1];
rz(-0.34653133) q[2];
sx q[2];
rz(-0.60014562) q[2];
sx q[2];
rz(-2.481386) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.79136363) q[1];
sx q[1];
rz(-1.204317) q[1];
sx q[1];
rz(0.54270737) q[1];
rz(1.8077778) q[3];
sx q[3];
rz(-1.6620043) q[3];
sx q[3];
rz(2.3438675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9875662) q[2];
sx q[2];
rz(-1.3904927) q[2];
sx q[2];
rz(-1.3898559) q[2];
rz(-3.0912345) q[3];
sx q[3];
rz(-1.7213089) q[3];
sx q[3];
rz(-2.0462947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20767009) q[0];
sx q[0];
rz(-1.727513) q[0];
sx q[0];
rz(0.06047824) q[0];
rz(-2.1531847) q[1];
sx q[1];
rz(-1.6840839) q[1];
sx q[1];
rz(0.21313937) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7558862) q[0];
sx q[0];
rz(-1.1637628) q[0];
sx q[0];
rz(-0.72188053) q[0];
rz(-pi) q[1];
rz(-2.1859841) q[2];
sx q[2];
rz(-0.19813523) q[2];
sx q[2];
rz(2.3096934) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7139529) q[1];
sx q[1];
rz(-2.8866077) q[1];
sx q[1];
rz(-0.1509576) q[1];
x q[2];
rz(1.074269) q[3];
sx q[3];
rz(-1.4098216) q[3];
sx q[3];
rz(2.9825008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9947784) q[2];
sx q[2];
rz(-1.5922567) q[2];
sx q[2];
rz(-0.0063304338) q[2];
rz(-1.7019042) q[3];
sx q[3];
rz(-2.3531745) q[3];
sx q[3];
rz(0.44062135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45563662) q[0];
sx q[0];
rz(-2.7485479) q[0];
sx q[0];
rz(1.8259557) q[0];
rz(1.3872967) q[1];
sx q[1];
rz(-0.60391128) q[1];
sx q[1];
rz(0.046262892) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4286908) q[0];
sx q[0];
rz(-1.5697456) q[0];
sx q[0];
rz(-0.029944972) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4607804) q[2];
sx q[2];
rz(-0.57475427) q[2];
sx q[2];
rz(-2.0865692) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8475593) q[1];
sx q[1];
rz(-2.1552076) q[1];
sx q[1];
rz(-1.1520034) q[1];
x q[2];
rz(2.9748232) q[3];
sx q[3];
rz(-2.178949) q[3];
sx q[3];
rz(1.2603354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5578654) q[2];
sx q[2];
rz(-1.1745619) q[2];
sx q[2];
rz(0.91090703) q[2];
rz(1.2298498) q[3];
sx q[3];
rz(-0.76084796) q[3];
sx q[3];
rz(2.7228444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8987027) q[0];
sx q[0];
rz(-1.3985343) q[0];
sx q[0];
rz(2.718495) q[0];
rz(-0.93370071) q[1];
sx q[1];
rz(-1.9405148) q[1];
sx q[1];
rz(0.43412128) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2796626) q[0];
sx q[0];
rz(-1.2796784) q[0];
sx q[0];
rz(-1.5345957) q[0];
x q[1];
rz(2.4814194) q[2];
sx q[2];
rz(-1.4453363) q[2];
sx q[2];
rz(-1.8769285) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0435183) q[1];
sx q[1];
rz(-1.6988108) q[1];
sx q[1];
rz(2.5780067) q[1];
rz(-pi) q[2];
rz(-2.0781804) q[3];
sx q[3];
rz(-0.23325167) q[3];
sx q[3];
rz(2.1824529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9551129) q[2];
sx q[2];
rz(-3.0035512) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56483018) q[0];
sx q[0];
rz(-0.63693988) q[0];
sx q[0];
rz(1.379396) q[0];
rz(0.55343902) q[1];
sx q[1];
rz(-1.6575927) q[1];
sx q[1];
rz(2.4310506) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64122771) q[0];
sx q[0];
rz(-2.3661748) q[0];
sx q[0];
rz(-1.3432608) q[0];
rz(-2.8997786) q[2];
sx q[2];
rz(-1.514957) q[2];
sx q[2];
rz(-2.9602891) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.64046265) q[1];
sx q[1];
rz(-1.8014596) q[1];
sx q[1];
rz(0.035236852) q[1];
rz(-pi) q[2];
x q[2];
rz(1.738461) q[3];
sx q[3];
rz(-2.1049771) q[3];
sx q[3];
rz(3.0287865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1826627) q[2];
sx q[2];
rz(-1.6561597) q[2];
sx q[2];
rz(0.043969285) q[2];
rz(2.8751657) q[3];
sx q[3];
rz(-0.1259585) q[3];
sx q[3];
rz(-2.9282279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5037395) q[0];
sx q[0];
rz(-1.5218691) q[0];
sx q[0];
rz(0.36852401) q[0];
rz(-2.7875994) q[1];
sx q[1];
rz(-2.0123672) q[1];
sx q[1];
rz(-1.0634208) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3388004) q[0];
sx q[0];
rz(-1.7110076) q[0];
sx q[0];
rz(0.012787435) q[0];
rz(1.2930128) q[2];
sx q[2];
rz(-1.4627224) q[2];
sx q[2];
rz(-1.9982823) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1527271) q[1];
sx q[1];
rz(-0.60040632) q[1];
sx q[1];
rz(0.78515782) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74875777) q[3];
sx q[3];
rz(-1.7084661) q[3];
sx q[3];
rz(-2.2434821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.64663681) q[2];
sx q[2];
rz(-1.2594624) q[2];
sx q[2];
rz(-1.4259526) q[2];
rz(2.2706653) q[3];
sx q[3];
rz(-2.0723074) q[3];
sx q[3];
rz(2.9072185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17212412) q[0];
sx q[0];
rz(-2.5775462) q[0];
sx q[0];
rz(3.1308351) q[0];
rz(2.5116008) q[1];
sx q[1];
rz(-1.0456345) q[1];
sx q[1];
rz(-0.90525544) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0650103) q[0];
sx q[0];
rz(-2.6040051) q[0];
sx q[0];
rz(-1.4437066) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33335061) q[2];
sx q[2];
rz(-2.0596189) q[2];
sx q[2];
rz(-2.6437063) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.39054444) q[1];
sx q[1];
rz(-0.52459756) q[1];
sx q[1];
rz(2.8538997) q[1];
rz(1.0338883) q[3];
sx q[3];
rz(-2.4126588) q[3];
sx q[3];
rz(1.1400853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4524727) q[2];
sx q[2];
rz(-1.8931188) q[2];
sx q[2];
rz(-3.0242237) q[2];
rz(-2.4145685) q[3];
sx q[3];
rz(-2.2456462) q[3];
sx q[3];
rz(2.0015543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3177976) q[0];
sx q[0];
rz(-2.0841086) q[0];
sx q[0];
rz(-0.29528883) q[0];
rz(-1.7501887) q[1];
sx q[1];
rz(-0.65381217) q[1];
sx q[1];
rz(0.48694912) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21509296) q[0];
sx q[0];
rz(-1.7679259) q[0];
sx q[0];
rz(1.6736945) q[0];
rz(-2.0854961) q[2];
sx q[2];
rz(-1.8499415) q[2];
sx q[2];
rz(-2.5893958) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0978969) q[1];
sx q[1];
rz(-1.1292895) q[1];
sx q[1];
rz(-2.312078) q[1];
rz(-pi) q[2];
rz(1.7510909) q[3];
sx q[3];
rz(-1.8337941) q[3];
sx q[3];
rz(-0.20151787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9140909) q[2];
sx q[2];
rz(-2.5469683) q[2];
sx q[2];
rz(1.5957069) q[2];
rz(-2.5325736) q[3];
sx q[3];
rz(-2.0273429) q[3];
sx q[3];
rz(-0.6685037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85536999) q[0];
sx q[0];
rz(-2.2384221) q[0];
sx q[0];
rz(1.5089996) q[0];
rz(-1.1434932) q[1];
sx q[1];
rz(-1.9805464) q[1];
sx q[1];
rz(-0.11018363) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1299898) q[0];
sx q[0];
rz(-2.5429259) q[0];
sx q[0];
rz(-0.18417872) q[0];
rz(-pi) q[1];
rz(2.8844933) q[2];
sx q[2];
rz(-1.1214702) q[2];
sx q[2];
rz(-1.7420235) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.33823925) q[1];
sx q[1];
rz(-0.88283006) q[1];
sx q[1];
rz(-2.0961876) q[1];
rz(-0.0024282495) q[3];
sx q[3];
rz(-2.6442827) q[3];
sx q[3];
rz(-2.8893472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1630359) q[2];
sx q[2];
rz(-1.0806934) q[2];
sx q[2];
rz(-2.2829096) q[2];
rz(-1.5052172) q[3];
sx q[3];
rz(-1.3450164) q[3];
sx q[3];
rz(-1.4428562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34864146) q[0];
sx q[0];
rz(-1.071278) q[0];
sx q[0];
rz(0.68565482) q[0];
rz(2.3220956) q[1];
sx q[1];
rz(-1.6623431) q[1];
sx q[1];
rz(-2.562838) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8775691) q[0];
sx q[0];
rz(-0.79791466) q[0];
sx q[0];
rz(-2.430738) q[0];
rz(-pi) q[1];
rz(-2.8943693) q[2];
sx q[2];
rz(-2.1237388) q[2];
sx q[2];
rz(1.8930815) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.006224) q[1];
sx q[1];
rz(-1.9510165) q[1];
sx q[1];
rz(-1.848741) q[1];
rz(-0.41865809) q[3];
sx q[3];
rz(-3.0221429) q[3];
sx q[3];
rz(2.2653013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8389575) q[2];
sx q[2];
rz(-0.18582782) q[2];
sx q[2];
rz(-0.1798943) q[2];
rz(2.5707865) q[3];
sx q[3];
rz(-2.4195318) q[3];
sx q[3];
rz(-0.65176898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8410692) q[0];
sx q[0];
rz(-0.64995926) q[0];
sx q[0];
rz(1.7267701) q[0];
rz(-0.48814804) q[1];
sx q[1];
rz(-0.5562677) q[1];
sx q[1];
rz(1.5604896) q[1];
rz(2.7753903) q[2];
sx q[2];
rz(-2.1696268) q[2];
sx q[2];
rz(0.10757797) q[2];
rz(-2.6030953) q[3];
sx q[3];
rz(-1.8168728) q[3];
sx q[3];
rz(-1.3088165) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
