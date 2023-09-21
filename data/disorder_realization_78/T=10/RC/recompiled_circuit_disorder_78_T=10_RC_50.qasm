OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5821563) q[0];
sx q[0];
rz(-0.59098935) q[0];
sx q[0];
rz(0.58340573) q[0];
rz(-0.18435873) q[1];
sx q[1];
rz(4.1252131) q[1];
sx q[1];
rz(10.317378) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6658202) q[0];
sx q[0];
rz(-1.2736397) q[0];
sx q[0];
rz(2.1509403) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66944389) q[2];
sx q[2];
rz(-1.9813683) q[2];
sx q[2];
rz(-2.4759811) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8436369) q[1];
sx q[1];
rz(-1.8814109) q[1];
sx q[1];
rz(-2.288504) q[1];
x q[2];
rz(-3.0953232) q[3];
sx q[3];
rz(-0.75152961) q[3];
sx q[3];
rz(2.3627594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9154174) q[2];
sx q[2];
rz(-1.6025275) q[2];
sx q[2];
rz(-1.9809451) q[2];
rz(-0.21696572) q[3];
sx q[3];
rz(-2.6187077) q[3];
sx q[3];
rz(-2.0863566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0579257) q[0];
sx q[0];
rz(-0.96494976) q[0];
sx q[0];
rz(0.57587409) q[0];
rz(1.2469762) q[1];
sx q[1];
rz(-1.2966825) q[1];
sx q[1];
rz(1.1670246) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28536797) q[0];
sx q[0];
rz(-1.8282187) q[0];
sx q[0];
rz(-0.060398922) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1738271) q[2];
sx q[2];
rz(-2.5644828) q[2];
sx q[2];
rz(-1.7259665) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.45707073) q[1];
sx q[1];
rz(-1.6324537) q[1];
sx q[1];
rz(-2.3350299) q[1];
rz(-2.0608276) q[3];
sx q[3];
rz(-2.1809275) q[3];
sx q[3];
rz(0.78435635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.033826753) q[2];
sx q[2];
rz(-1.1788538) q[2];
sx q[2];
rz(2.9272184) q[2];
rz(0.073444627) q[3];
sx q[3];
rz(-0.44973222) q[3];
sx q[3];
rz(0.28081056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6630994) q[0];
sx q[0];
rz(-2.5254624) q[0];
sx q[0];
rz(-1.2269155) q[0];
rz(-2.7413209) q[1];
sx q[1];
rz(-1.2534671) q[1];
sx q[1];
rz(1.0148369) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0600216) q[0];
sx q[0];
rz(-0.9284174) q[0];
sx q[0];
rz(-1.6409671) q[0];
rz(-2.8995908) q[2];
sx q[2];
rz(-1.9654462) q[2];
sx q[2];
rz(1.637527) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.2753607) q[1];
sx q[1];
rz(-0.97517698) q[1];
sx q[1];
rz(0.097252107) q[1];
rz(-pi) q[2];
rz(2.879911) q[3];
sx q[3];
rz(-2.2398584) q[3];
sx q[3];
rz(-2.3467968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0959452) q[2];
sx q[2];
rz(-1.0568876) q[2];
sx q[2];
rz(-0.8992368) q[2];
rz(-2.4441161) q[3];
sx q[3];
rz(-1.2858425) q[3];
sx q[3];
rz(-2.5337059) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65524453) q[0];
sx q[0];
rz(-1.0826033) q[0];
sx q[0];
rz(2.7096601) q[0];
rz(-0.63255429) q[1];
sx q[1];
rz(-2.7245941) q[1];
sx q[1];
rz(0.63582173) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0441372) q[0];
sx q[0];
rz(-2.6090528) q[0];
sx q[0];
rz(1.2116648) q[0];
rz(-pi) q[1];
rz(-2.9767838) q[2];
sx q[2];
rz(-0.85068446) q[2];
sx q[2];
rz(-2.1317496) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7829261) q[1];
sx q[1];
rz(-1.7554605) q[1];
sx q[1];
rz(-0.77628805) q[1];
rz(-pi) q[2];
rz(-2.5084247) q[3];
sx q[3];
rz(-2.79106) q[3];
sx q[3];
rz(-1.4299973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2137961) q[2];
sx q[2];
rz(-2.9979604) q[2];
sx q[2];
rz(0.68871838) q[2];
rz(-0.33411807) q[3];
sx q[3];
rz(-1.170661) q[3];
sx q[3];
rz(-0.14373246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14389811) q[0];
sx q[0];
rz(-2.4425638) q[0];
sx q[0];
rz(-2.0671663) q[0];
rz(2.396446) q[1];
sx q[1];
rz(-1.4829758) q[1];
sx q[1];
rz(0.27854663) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6313948) q[0];
sx q[0];
rz(-0.83980951) q[0];
sx q[0];
rz(0.98548074) q[0];
rz(1.4611545) q[2];
sx q[2];
rz(-2.1796558) q[2];
sx q[2];
rz(1.4102175) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.46036938) q[1];
sx q[1];
rz(-1.2037828) q[1];
sx q[1];
rz(-2.2296434) q[1];
rz(-pi) q[2];
rz(2.1513125) q[3];
sx q[3];
rz(-2.2970082) q[3];
sx q[3];
rz(-2.5648404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4429861) q[2];
sx q[2];
rz(-1.7437982) q[2];
sx q[2];
rz(0.29423514) q[2];
rz(-0.081929835) q[3];
sx q[3];
rz(-0.51920813) q[3];
sx q[3];
rz(-0.036227139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1086403) q[0];
sx q[0];
rz(-2.2886798) q[0];
sx q[0];
rz(3.1325353) q[0];
rz(-2.5065705) q[1];
sx q[1];
rz(-2.451684) q[1];
sx q[1];
rz(3.0335398) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8979643) q[0];
sx q[0];
rz(-1.122323) q[0];
sx q[0];
rz(-1.2178221) q[0];
rz(-pi) q[1];
rz(-2.8073505) q[2];
sx q[2];
rz(-2.7396482) q[2];
sx q[2];
rz(0.7451171) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6335771) q[1];
sx q[1];
rz(-1.0883691) q[1];
sx q[1];
rz(1.1303933) q[1];
rz(-pi) q[2];
rz(1.990854) q[3];
sx q[3];
rz(-2.5541411) q[3];
sx q[3];
rz(-2.1849039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1487427) q[2];
sx q[2];
rz(-1.7603346) q[2];
sx q[2];
rz(-1.2711058) q[2];
rz(0.078401119) q[3];
sx q[3];
rz(-1.4784808) q[3];
sx q[3];
rz(-1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(0.25061297) q[0];
sx q[0];
rz(-1.5654634) q[0];
sx q[0];
rz(2.4839731) q[0];
rz(1.3972067) q[1];
sx q[1];
rz(-1.9669292) q[1];
sx q[1];
rz(-2.2479642) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4511787) q[0];
sx q[0];
rz(-0.72890857) q[0];
sx q[0];
rz(0.50509392) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29427476) q[2];
sx q[2];
rz(-1.1106967) q[2];
sx q[2];
rz(1.0690881) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0221755) q[1];
sx q[1];
rz(-0.70964538) q[1];
sx q[1];
rz(-2.8303353) q[1];
rz(-pi) q[2];
rz(-2.3818447) q[3];
sx q[3];
rz(-1.500251) q[3];
sx q[3];
rz(2.7618046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.39067337) q[2];
sx q[2];
rz(-2.1942287) q[2];
sx q[2];
rz(-0.52948362) q[2];
rz(-0.47618619) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(0.34255323) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7628409) q[0];
sx q[0];
rz(-1.6289926) q[0];
sx q[0];
rz(-1.09028) q[0];
rz(0.11225637) q[1];
sx q[1];
rz(-2.0376164) q[1];
sx q[1];
rz(1.1539248) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8114416) q[0];
sx q[0];
rz(-2.220287) q[0];
sx q[0];
rz(-1.6135471) q[0];
rz(-pi) q[1];
rz(0.55118982) q[2];
sx q[2];
rz(-1.4442208) q[2];
sx q[2];
rz(1.7984496) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0358419) q[1];
sx q[1];
rz(-0.81139794) q[1];
sx q[1];
rz(1.6559385) q[1];
x q[2];
rz(2.7525547) q[3];
sx q[3];
rz(-2.2044047) q[3];
sx q[3];
rz(-2.6694359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5899137) q[2];
sx q[2];
rz(-2.7605197) q[2];
sx q[2];
rz(-0.38044688) q[2];
rz(-1.1278661) q[3];
sx q[3];
rz(-1.7815855) q[3];
sx q[3];
rz(-1.5650704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8001051) q[0];
sx q[0];
rz(-1.8999758) q[0];
sx q[0];
rz(-0.63968101) q[0];
rz(1.9027963) q[1];
sx q[1];
rz(-1.7265373) q[1];
sx q[1];
rz(-1.170084) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3568048) q[0];
sx q[0];
rz(-2.6492282) q[0];
sx q[0];
rz(2.3193293) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8662888) q[2];
sx q[2];
rz(-0.79553662) q[2];
sx q[2];
rz(0.77992935) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.77800345) q[1];
sx q[1];
rz(-2.6215141) q[1];
sx q[1];
rz(2.7705454) q[1];
x q[2];
rz(1.3977259) q[3];
sx q[3];
rz(-0.67865463) q[3];
sx q[3];
rz(-2.9242587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8273948) q[2];
sx q[2];
rz(-2.4908227) q[2];
sx q[2];
rz(-2.7588552) q[2];
rz(0.9283723) q[3];
sx q[3];
rz(-1.9675156) q[3];
sx q[3];
rz(-0.66463566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2255573) q[0];
sx q[0];
rz(-1.6397497) q[0];
sx q[0];
rz(2.8826707) q[0];
rz(-2.4312773) q[1];
sx q[1];
rz(-1.0537035) q[1];
sx q[1];
rz(-0.47992596) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5537162) q[0];
sx q[0];
rz(-1.237861) q[0];
sx q[0];
rz(-0.90162189) q[0];
rz(-pi) q[1];
rz(2.2357113) q[2];
sx q[2];
rz(-1.3833628) q[2];
sx q[2];
rz(-2.4311709) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2952134) q[1];
sx q[1];
rz(-1.0615674) q[1];
sx q[1];
rz(1.1649811) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5273401) q[3];
sx q[3];
rz(-0.96685997) q[3];
sx q[3];
rz(0.37249836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.84247983) q[2];
sx q[2];
rz(-0.92876902) q[2];
sx q[2];
rz(1.3170362) q[2];
rz(1.2420098) q[3];
sx q[3];
rz(-0.95359355) q[3];
sx q[3];
rz(2.4035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9823572) q[0];
sx q[0];
rz(-3.1095105) q[0];
sx q[0];
rz(1.6788917) q[0];
rz(2.1622529) q[1];
sx q[1];
rz(-2.0420488) q[1];
sx q[1];
rz(2.2534823) q[1];
rz(-1.3946891) q[2];
sx q[2];
rz(-2.0460143) q[2];
sx q[2];
rz(2.569414) q[2];
rz(1.3027719) q[3];
sx q[3];
rz(-1.839073) q[3];
sx q[3];
rz(3.0243235) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];