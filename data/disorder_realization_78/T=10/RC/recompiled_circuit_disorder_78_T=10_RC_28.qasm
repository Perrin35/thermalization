OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.55943638) q[0];
sx q[0];
rz(3.732582) q[0];
sx q[0];
rz(8.8413722) q[0];
rz(2.9572339) q[1];
sx q[1];
rz(-0.98362041) q[1];
sx q[1];
rz(-0.89259994) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8569782) q[0];
sx q[0];
rz(-2.1224788) q[0];
sx q[0];
rz(-2.7906228) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0775665) q[2];
sx q[2];
rz(-2.1760586) q[2];
sx q[2];
rz(-1.93047) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6751911) q[1];
sx q[1];
rz(-2.2474504) q[1];
sx q[1];
rz(0.40282175) q[1];
rz(-pi) q[2];
rz(-1.5276018) q[3];
sx q[3];
rz(-0.82026635) q[3];
sx q[3];
rz(-2.299472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2261752) q[2];
sx q[2];
rz(-1.6025275) q[2];
sx q[2];
rz(-1.1606476) q[2];
rz(0.21696572) q[3];
sx q[3];
rz(-0.522885) q[3];
sx q[3];
rz(-2.0863566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.083667) q[0];
sx q[0];
rz(-0.96494976) q[0];
sx q[0];
rz(-2.5657186) q[0];
rz(1.8946164) q[1];
sx q[1];
rz(-1.8449102) q[1];
sx q[1];
rz(-1.974568) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8407699) q[0];
sx q[0];
rz(-1.6292028) q[0];
sx q[0];
rz(-1.3129243) q[0];
rz(-pi) q[1];
rz(-2.8950047) q[2];
sx q[2];
rz(-2.0980667) q[2];
sx q[2];
rz(1.2621244) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6845219) q[1];
sx q[1];
rz(-1.6324537) q[1];
sx q[1];
rz(0.80656273) q[1];
rz(-pi) q[2];
rz(2.4715273) q[3];
sx q[3];
rz(-1.1747922) q[3];
sx q[3];
rz(2.0585287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1077659) q[2];
sx q[2];
rz(-1.1788538) q[2];
sx q[2];
rz(-2.9272184) q[2];
rz(0.073444627) q[3];
sx q[3];
rz(-0.44973222) q[3];
sx q[3];
rz(-2.8607821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6630994) q[0];
sx q[0];
rz(-0.61613023) q[0];
sx q[0];
rz(-1.2269155) q[0];
rz(-2.7413209) q[1];
sx q[1];
rz(-1.8881256) q[1];
sx q[1];
rz(-1.0148369) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6944511) q[0];
sx q[0];
rz(-1.6269636) q[0];
sx q[0];
rz(0.64356128) q[0];
rz(0.24200183) q[2];
sx q[2];
rz(-1.1761464) q[2];
sx q[2];
rz(1.5040656) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.44753578) q[1];
sx q[1];
rz(-0.60255614) q[1];
sx q[1];
rz(1.4284929) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.879911) q[3];
sx q[3];
rz(-0.90173429) q[3];
sx q[3];
rz(0.79479587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0959452) q[2];
sx q[2];
rz(-1.0568876) q[2];
sx q[2];
rz(2.2423559) q[2];
rz(0.69747654) q[3];
sx q[3];
rz(-1.8557502) q[3];
sx q[3];
rz(-0.60788679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65524453) q[0];
sx q[0];
rz(-2.0589893) q[0];
sx q[0];
rz(0.43193257) q[0];
rz(-2.5090384) q[1];
sx q[1];
rz(-0.41699854) q[1];
sx q[1];
rz(-2.5057709) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0974554) q[0];
sx q[0];
rz(-0.53253981) q[0];
sx q[0];
rz(-1.9299279) q[0];
rz(-0.16480883) q[2];
sx q[2];
rz(-0.85068446) q[2];
sx q[2];
rz(-1.009843) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1078474) q[1];
sx q[1];
rz(-2.3305232) q[1];
sx q[1];
rz(1.3147522) q[1];
x q[2];
rz(-0.63316791) q[3];
sx q[3];
rz(-2.79106) q[3];
sx q[3];
rz(1.4299973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2137961) q[2];
sx q[2];
rz(-2.9979604) q[2];
sx q[2];
rz(2.4528743) q[2];
rz(-2.8074746) q[3];
sx q[3];
rz(-1.170661) q[3];
sx q[3];
rz(-2.9978602) q[3];
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
rz(-2.9976945) q[0];
sx q[0];
rz(-2.4425638) q[0];
sx q[0];
rz(-1.0744263) q[0];
rz(-2.396446) q[1];
sx q[1];
rz(-1.6586168) q[1];
sx q[1];
rz(-2.863046) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6644088) q[0];
sx q[0];
rz(-1.9946949) q[0];
sx q[0];
rz(0.82188481) q[0];
rz(2.5299046) q[2];
sx q[2];
rz(-1.4809161) q[2];
sx q[2];
rz(2.9181366) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3814195) q[1];
sx q[1];
rz(-0.96254327) q[1];
sx q[1];
rz(0.45254032) q[1];
rz(-pi) q[2];
rz(-0.81551084) q[3];
sx q[3];
rz(-1.9933356) q[3];
sx q[3];
rz(2.5583207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4429861) q[2];
sx q[2];
rz(-1.7437982) q[2];
sx q[2];
rz(0.29423514) q[2];
rz(-0.081929835) q[3];
sx q[3];
rz(-2.6223845) q[3];
sx q[3];
rz(0.036227139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1086403) q[0];
sx q[0];
rz(-0.8529129) q[0];
sx q[0];
rz(3.1325353) q[0];
rz(-2.5065705) q[1];
sx q[1];
rz(-0.68990866) q[1];
sx q[1];
rz(0.10805282) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2436284) q[0];
sx q[0];
rz(-2.0192696) q[0];
sx q[0];
rz(-1.2178221) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4322386) q[2];
sx q[2];
rz(-1.9493305) q[2];
sx q[2];
rz(-1.1059424) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3012645) q[1];
sx q[1];
rz(-2.5003308) q[1];
sx q[1];
rz(0.68323369) q[1];
x q[2];
rz(-1.990854) q[3];
sx q[3];
rz(-2.5541411) q[3];
sx q[3];
rz(2.1849039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1487427) q[2];
sx q[2];
rz(-1.7603346) q[2];
sx q[2];
rz(1.2711058) q[2];
rz(-0.078401119) q[3];
sx q[3];
rz(-1.6631118) q[3];
sx q[3];
rz(2.0097849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.25061297) q[0];
sx q[0];
rz(-1.5654634) q[0];
sx q[0];
rz(-0.65761956) q[0];
rz(1.3972067) q[1];
sx q[1];
rz(-1.9669292) q[1];
sx q[1];
rz(-2.2479642) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69041396) q[0];
sx q[0];
rz(-0.72890857) q[0];
sx q[0];
rz(0.50509392) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0486084) q[2];
sx q[2];
rz(-1.8336979) q[2];
sx q[2];
rz(2.773657) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.5205198) q[1];
sx q[1];
rz(-2.2398661) q[1];
sx q[1];
rz(1.8280162) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4736389) q[3];
sx q[3];
rz(-2.328184) q[3];
sx q[3];
rz(-2.017445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39067337) q[2];
sx q[2];
rz(-0.94736391) q[2];
sx q[2];
rz(-0.52948362) q[2];
rz(2.6654065) q[3];
sx q[3];
rz(-1.3323077) q[3];
sx q[3];
rz(-0.34255323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7628409) q[0];
sx q[0];
rz(-1.6289926) q[0];
sx q[0];
rz(1.09028) q[0];
rz(-0.11225637) q[1];
sx q[1];
rz(-1.1039762) q[1];
sx q[1];
rz(-1.9876678) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40076462) q[0];
sx q[0];
rz(-2.4908998) q[0];
sx q[0];
rz(-0.056218938) q[0];
rz(-pi) q[1];
rz(0.55118982) q[2];
sx q[2];
rz(-1.6973719) q[2];
sx q[2];
rz(-1.7984496) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9824144) q[1];
sx q[1];
rz(-0.7632066) q[1];
sx q[1];
rz(0.089341954) q[1];
rz(-pi) q[2];
rz(-2.7525547) q[3];
sx q[3];
rz(-2.2044047) q[3];
sx q[3];
rz(2.6694359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.551679) q[2];
sx q[2];
rz(-2.7605197) q[2];
sx q[2];
rz(0.38044688) q[2];
rz(-2.0137265) q[3];
sx q[3];
rz(-1.3600072) q[3];
sx q[3];
rz(1.5765223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8001051) q[0];
sx q[0];
rz(-1.8999758) q[0];
sx q[0];
rz(-0.63968101) q[0];
rz(1.2387964) q[1];
sx q[1];
rz(-1.4150554) q[1];
sx q[1];
rz(-1.170084) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6696475) q[0];
sx q[0];
rz(-1.8983316) q[0];
sx q[0];
rz(-1.1963084) q[0];
x q[1];
rz(2.8662888) q[2];
sx q[2];
rz(-2.346056) q[2];
sx q[2];
rz(2.3616633) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.46718405) q[1];
sx q[1];
rz(-1.7519752) q[1];
sx q[1];
rz(0.49023899) q[1];
x q[2];
rz(0.13799237) q[3];
sx q[3];
rz(-2.237461) q[3];
sx q[3];
rz(-0.43825144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8273948) q[2];
sx q[2];
rz(-0.65076995) q[2];
sx q[2];
rz(-0.38273746) q[2];
rz(-0.9283723) q[3];
sx q[3];
rz(-1.9675156) q[3];
sx q[3];
rz(-2.476957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2255573) q[0];
sx q[0];
rz(-1.501843) q[0];
sx q[0];
rz(-2.8826707) q[0];
rz(-2.4312773) q[1];
sx q[1];
rz(-1.0537035) q[1];
sx q[1];
rz(2.6616667) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5537162) q[0];
sx q[0];
rz(-1.9037316) q[0];
sx q[0];
rz(-0.90162189) q[0];
rz(-1.8690228) q[2];
sx q[2];
rz(-0.68694653) q[2];
sx q[2];
rz(0.62703122) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.5727947) q[1];
sx q[1];
rz(-0.63981445) q[1];
sx q[1];
rz(-0.6154284) q[1];
rz(-pi) q[2];
rz(0.89703538) q[3];
sx q[3];
rz(-1.1437136) q[3];
sx q[3];
rz(0.87891146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.84247983) q[2];
sx q[2];
rz(-0.92876902) q[2];
sx q[2];
rz(1.8245565) q[2];
rz(-1.2420098) q[3];
sx q[3];
rz(-0.95359355) q[3];
sx q[3];
rz(0.73808134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9823572) q[0];
sx q[0];
rz(-0.032082162) q[0];
sx q[0];
rz(-1.4627009) q[0];
rz(2.1622529) q[1];
sx q[1];
rz(-2.0420488) q[1];
sx q[1];
rz(2.2534823) q[1];
rz(2.6600044) q[2];
sx q[2];
rz(-1.7272186) q[2];
sx q[2];
rz(-2.061736) q[2];
rz(1.8388207) q[3];
sx q[3];
rz(-1.3025197) q[3];
sx q[3];
rz(-0.11726911) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];