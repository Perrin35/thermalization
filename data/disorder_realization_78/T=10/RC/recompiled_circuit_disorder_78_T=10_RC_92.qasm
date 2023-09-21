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
rz(-2.5581869) q[0];
rz(2.9572339) q[1];
sx q[1];
rz(-0.98362041) q[1];
sx q[1];
rz(2.2489927) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4757724) q[0];
sx q[0];
rz(-1.8679529) q[0];
sx q[0];
rz(2.1509403) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5298654) q[2];
sx q[2];
rz(-0.76843843) q[2];
sx q[2];
rz(-0.4380463) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.46640151) q[1];
sx q[1];
rz(-0.89414222) q[1];
sx q[1];
rz(-2.7387709) q[1];
rz(-0.75099545) q[3];
sx q[3];
rz(-1.5392116) q[3];
sx q[3];
rz(-0.75814523) q[3];
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
rz(0.21696572) q[3];
sx q[3];
rz(-2.6187077) q[3];
sx q[3];
rz(-1.0552361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.083667) q[0];
sx q[0];
rz(-2.1766429) q[0];
sx q[0];
rz(0.57587409) q[0];
rz(-1.8946164) q[1];
sx q[1];
rz(-1.2966825) q[1];
sx q[1];
rz(-1.974568) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8407699) q[0];
sx q[0];
rz(-1.6292028) q[0];
sx q[0];
rz(1.8286684) q[0];
rz(-pi) q[1];
rz(1.9677656) q[2];
sx q[2];
rz(-2.5644828) q[2];
sx q[2];
rz(-1.7259665) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0547304) q[1];
sx q[1];
rz(-2.3332101) q[1];
sx q[1];
rz(-0.085309172) q[1];
rz(-2.5490709) q[3];
sx q[3];
rz(-2.3791109) q[3];
sx q[3];
rz(-3.1069063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1077659) q[2];
sx q[2];
rz(-1.9627389) q[2];
sx q[2];
rz(0.21437422) q[2];
rz(0.073444627) q[3];
sx q[3];
rz(-0.44973222) q[3];
sx q[3];
rz(-2.8607821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4784933) q[0];
sx q[0];
rz(-2.5254624) q[0];
sx q[0];
rz(-1.9146772) q[0];
rz(0.40027174) q[1];
sx q[1];
rz(-1.2534671) q[1];
sx q[1];
rz(-2.1267557) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0600216) q[0];
sx q[0];
rz(-2.2131753) q[0];
sx q[0];
rz(1.5006256) q[0];
rz(-pi) q[1];
rz(2.0929167) q[2];
sx q[2];
rz(-0.45959696) q[2];
sx q[2];
rz(-0.93333474) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.240757) q[1];
sx q[1];
rz(-1.6512617) q[1];
sx q[1];
rz(2.1686173) q[1];
x q[2];
rz(2.879911) q[3];
sx q[3];
rz(-2.2398584) q[3];
sx q[3];
rz(0.79479587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0456475) q[2];
sx q[2];
rz(-1.0568876) q[2];
sx q[2];
rz(-2.2423559) q[2];
rz(-2.4441161) q[3];
sx q[3];
rz(-1.8557502) q[3];
sx q[3];
rz(-0.60788679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65524453) q[0];
sx q[0];
rz(-1.0826033) q[0];
sx q[0];
rz(-2.7096601) q[0];
rz(0.63255429) q[1];
sx q[1];
rz(-0.41699854) q[1];
sx q[1];
rz(0.63582173) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0974554) q[0];
sx q[0];
rz(-0.53253981) q[0];
sx q[0];
rz(1.2116648) q[0];
rz(-pi) q[1];
rz(0.16480883) q[2];
sx q[2];
rz(-0.85068446) q[2];
sx q[2];
rz(1.009843) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7829261) q[1];
sx q[1];
rz(-1.7554605) q[1];
sx q[1];
rz(2.3653046) q[1];
rz(-pi) q[2];
rz(-0.28663978) q[3];
sx q[3];
rz(-1.3661824) q[3];
sx q[3];
rz(-2.3972547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9277966) q[2];
sx q[2];
rz(-2.9979604) q[2];
sx q[2];
rz(0.68871838) q[2];
rz(2.8074746) q[3];
sx q[3];
rz(-1.170661) q[3];
sx q[3];
rz(-0.14373246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14389811) q[0];
sx q[0];
rz(-2.4425638) q[0];
sx q[0];
rz(1.0744263) q[0];
rz(2.396446) q[1];
sx q[1];
rz(-1.4829758) q[1];
sx q[1];
rz(0.27854663) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8700096) q[0];
sx q[0];
rz(-0.9013114) q[0];
sx q[0];
rz(0.55218009) q[0];
x q[1];
rz(-2.9859221) q[2];
sx q[2];
rz(-0.61741932) q[2];
sx q[2];
rz(1.9215259) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4651471) q[1];
sx q[1];
rz(-2.4009581) q[1];
sx q[1];
rz(2.1315104) q[1];
x q[2];
rz(-0.99028011) q[3];
sx q[3];
rz(-0.84458447) q[3];
sx q[3];
rz(2.5648404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4429861) q[2];
sx q[2];
rz(-1.3977945) q[2];
sx q[2];
rz(0.29423514) q[2];
rz(-3.0596628) q[3];
sx q[3];
rz(-0.51920813) q[3];
sx q[3];
rz(-3.1053655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.0329523) q[0];
sx q[0];
rz(-0.8529129) q[0];
sx q[0];
rz(-0.0090573514) q[0];
rz(0.63502216) q[1];
sx q[1];
rz(-2.451684) q[1];
sx q[1];
rz(-0.10805282) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2436284) q[0];
sx q[0];
rz(-2.0192696) q[0];
sx q[0];
rz(1.2178221) q[0];
rz(-1.7093541) q[2];
sx q[2];
rz(-1.1922622) q[2];
sx q[2];
rz(-1.1059424) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5080155) q[1];
sx q[1];
rz(-1.0883691) q[1];
sx q[1];
rz(2.0111994) q[1];
rz(1.024527) q[3];
sx q[3];
rz(-1.7987935) q[3];
sx q[3];
rz(-2.8834164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.99284995) q[2];
sx q[2];
rz(-1.381258) q[2];
sx q[2];
rz(-1.8704869) q[2];
rz(3.0631915) q[3];
sx q[3];
rz(-1.4784808) q[3];
sx q[3];
rz(-2.0097849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25061297) q[0];
sx q[0];
rz(-1.5761292) q[0];
sx q[0];
rz(-2.4839731) q[0];
rz(1.744386) q[1];
sx q[1];
rz(-1.9669292) q[1];
sx q[1];
rz(2.2479642) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4511787) q[0];
sx q[0];
rz(-0.72890857) q[0];
sx q[0];
rz(-2.6364987) q[0];
rz(-0.29427476) q[2];
sx q[2];
rz(-1.1106967) q[2];
sx q[2];
rz(1.0690881) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5205198) q[1];
sx q[1];
rz(-0.90172651) q[1];
sx q[1];
rz(-1.8280162) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0393533) q[3];
sx q[3];
rz(-2.3792301) q[3];
sx q[3];
rz(-1.2650714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7509193) q[2];
sx q[2];
rz(-0.94736391) q[2];
sx q[2];
rz(0.52948362) q[2];
rz(0.47618619) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(-0.34255323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7628409) q[0];
sx q[0];
rz(-1.6289926) q[0];
sx q[0];
rz(2.0513127) q[0];
rz(-3.0293363) q[1];
sx q[1];
rz(-1.1039762) q[1];
sx q[1];
rz(1.9876678) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9268122) q[0];
sx q[0];
rz(-1.5367537) q[0];
sx q[0];
rz(-2.4916617) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5904028) q[2];
sx q[2];
rz(-1.4442208) q[2];
sx q[2];
rz(1.3431431) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.47626074) q[1];
sx q[1];
rz(-1.5090824) q[1];
sx q[1];
rz(0.76121059) q[1];
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
rz(-0.38107291) q[2];
sx q[2];
rz(2.7611458) q[2];
rz(1.1278661) q[3];
sx q[3];
rz(-1.3600072) q[3];
sx q[3];
rz(-1.5650704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8001051) q[0];
sx q[0];
rz(-1.2416168) q[0];
sx q[0];
rz(-0.63968101) q[0];
rz(-1.2387964) q[1];
sx q[1];
rz(-1.4150554) q[1];
sx q[1];
rz(-1.9715086) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4719452) q[0];
sx q[0];
rz(-1.2432611) q[0];
sx q[0];
rz(-1.1963084) q[0];
rz(-2.8662888) q[2];
sx q[2];
rz(-0.79553662) q[2];
sx q[2];
rz(2.3616633) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3635892) q[1];
sx q[1];
rz(-0.52007857) q[1];
sx q[1];
rz(2.7705454) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2421078) q[3];
sx q[3];
rz(-1.6791108) q[3];
sx q[3];
rz(1.9233821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8273948) q[2];
sx q[2];
rz(-0.65076995) q[2];
sx q[2];
rz(-2.7588552) q[2];
rz(-2.2132204) q[3];
sx q[3];
rz(-1.174077) q[3];
sx q[3];
rz(-2.476957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2255573) q[0];
sx q[0];
rz(-1.501843) q[0];
sx q[0];
rz(-0.25892192) q[0];
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
rz(1.5878764) q[0];
sx q[0];
rz(-1.9037316) q[0];
sx q[0];
rz(0.90162189) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9051022) q[2];
sx q[2];
rz(-0.91954008) q[2];
sx q[2];
rz(1.0054393) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.568798) q[1];
sx q[1];
rz(-0.63981445) q[1];
sx q[1];
rz(-2.5261643) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94060937) q[3];
sx q[3];
rz(-0.77946957) q[3];
sx q[3];
rz(-1.1704695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2991128) q[2];
sx q[2];
rz(-2.2128236) q[2];
sx q[2];
rz(1.8245565) q[2];
rz(-1.2420098) q[3];
sx q[3];
rz(-2.1879991) q[3];
sx q[3];
rz(2.4035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15923545) q[0];
sx q[0];
rz(-0.032082162) q[0];
sx q[0];
rz(-1.4627009) q[0];
rz(-2.1622529) q[1];
sx q[1];
rz(-1.0995438) q[1];
sx q[1];
rz(-0.88811036) q[1];
rz(1.3946891) q[2];
sx q[2];
rz(-1.0955784) q[2];
sx q[2];
rz(-0.57217862) q[2];
rz(-1.8388207) q[3];
sx q[3];
rz(-1.839073) q[3];
sx q[3];
rz(3.0243235) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
