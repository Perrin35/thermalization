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
rz(-0.18435873) q[1];
sx q[1];
rz(-2.1579722) q[1];
sx q[1];
rz(0.89259994) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8569782) q[0];
sx q[0];
rz(-1.0191139) q[0];
sx q[0];
rz(0.3509699) q[0];
x q[1];
rz(2.4721488) q[2];
sx q[2];
rz(-1.9813683) q[2];
sx q[2];
rz(0.66561156) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6751911) q[1];
sx q[1];
rz(-0.89414222) q[1];
sx q[1];
rz(2.7387709) q[1];
rz(-0.75099545) q[3];
sx q[3];
rz(-1.6023811) q[3];
sx q[3];
rz(-2.3834474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9154174) q[2];
sx q[2];
rz(-1.5390652) q[2];
sx q[2];
rz(-1.9809451) q[2];
rz(2.9246269) q[3];
sx q[3];
rz(-0.522885) q[3];
sx q[3];
rz(-1.0552361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0579257) q[0];
sx q[0];
rz(-2.1766429) q[0];
sx q[0];
rz(-2.5657186) q[0];
rz(-1.2469762) q[1];
sx q[1];
rz(-1.8449102) q[1];
sx q[1];
rz(-1.974568) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3008227) q[0];
sx q[0];
rz(-1.6292028) q[0];
sx q[0];
rz(-1.3129243) q[0];
rz(-pi) q[1];
rz(-2.1115233) q[2];
sx q[2];
rz(-1.7833372) q[2];
sx q[2];
rz(2.9589047) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0547304) q[1];
sx q[1];
rz(-0.80838258) q[1];
sx q[1];
rz(3.0562835) q[1];
x q[2];
rz(-2.0608276) q[3];
sx q[3];
rz(-2.1809275) q[3];
sx q[3];
rz(0.78435635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1077659) q[2];
sx q[2];
rz(-1.9627389) q[2];
sx q[2];
rz(-2.9272184) q[2];
rz(-3.068148) q[3];
sx q[3];
rz(-2.6918604) q[3];
sx q[3];
rz(-0.28081056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4784933) q[0];
sx q[0];
rz(-2.5254624) q[0];
sx q[0];
rz(-1.9146772) q[0];
rz(2.7413209) q[1];
sx q[1];
rz(-1.2534671) q[1];
sx q[1];
rz(2.1267557) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0600216) q[0];
sx q[0];
rz(-0.9284174) q[0];
sx q[0];
rz(1.5006256) q[0];
x q[1];
rz(-1.048676) q[2];
sx q[2];
rz(-2.6819957) q[2];
sx q[2];
rz(-2.2082579) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.866232) q[1];
sx q[1];
rz(-0.97517698) q[1];
sx q[1];
rz(3.0443405) q[1];
rz(-0.26168163) q[3];
sx q[3];
rz(-0.90173429) q[3];
sx q[3];
rz(2.3467968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0456475) q[2];
sx q[2];
rz(-1.0568876) q[2];
sx q[2];
rz(-0.8992368) q[2];
rz(2.4441161) q[3];
sx q[3];
rz(-1.2858425) q[3];
sx q[3];
rz(2.5337059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65524453) q[0];
sx q[0];
rz(-2.0589893) q[0];
sx q[0];
rz(-0.43193257) q[0];
rz(2.5090384) q[1];
sx q[1];
rz(-2.7245941) q[1];
sx q[1];
rz(-2.5057709) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0974554) q[0];
sx q[0];
rz(-0.53253981) q[0];
sx q[0];
rz(-1.2116648) q[0];
rz(2.9767838) q[2];
sx q[2];
rz(-0.85068446) q[2];
sx q[2];
rz(-1.009843) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.35866657) q[1];
sx q[1];
rz(-1.7554605) q[1];
sx q[1];
rz(0.77628805) q[1];
rz(-0.63316791) q[3];
sx q[3];
rz(-2.79106) q[3];
sx q[3];
rz(1.4299973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2137961) q[2];
sx q[2];
rz(-0.14363229) q[2];
sx q[2];
rz(-0.68871838) q[2];
rz(-2.8074746) q[3];
sx q[3];
rz(-1.170661) q[3];
sx q[3];
rz(-2.9978602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9976945) q[0];
sx q[0];
rz(-2.4425638) q[0];
sx q[0];
rz(-1.0744263) q[0];
rz(-0.74514666) q[1];
sx q[1];
rz(-1.4829758) q[1];
sx q[1];
rz(-2.863046) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6313948) q[0];
sx q[0];
rz(-2.3017831) q[0];
sx q[0];
rz(-2.1561119) q[0];
rz(1.4611545) q[2];
sx q[2];
rz(-2.1796558) q[2];
sx q[2];
rz(1.4102175) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7601732) q[1];
sx q[1];
rz(-0.96254327) q[1];
sx q[1];
rz(-2.6890523) q[1];
rz(-pi) q[2];
rz(-0.55322247) q[3];
sx q[3];
rz(-2.2461938) q[3];
sx q[3];
rz(1.355987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4429861) q[2];
sx q[2];
rz(-1.3977945) q[2];
sx q[2];
rz(2.8473575) q[2];
rz(3.0596628) q[3];
sx q[3];
rz(-0.51920813) q[3];
sx q[3];
rz(-0.036227139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1086403) q[0];
sx q[0];
rz(-2.2886798) q[0];
sx q[0];
rz(-3.1325353) q[0];
rz(-0.63502216) q[1];
sx q[1];
rz(-0.68990866) q[1];
sx q[1];
rz(-0.10805282) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51443433) q[0];
sx q[0];
rz(-1.8875727) q[0];
sx q[0];
rz(2.6677674) q[0];
rz(-pi) q[1];
rz(0.33424218) q[2];
sx q[2];
rz(-2.7396482) q[2];
sx q[2];
rz(-2.3964756) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5080155) q[1];
sx q[1];
rz(-1.0883691) q[1];
sx q[1];
rz(-2.0111994) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26515682) q[3];
sx q[3];
rz(-1.0401871) q[3];
sx q[3];
rz(1.6924072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.99284995) q[2];
sx q[2];
rz(-1.381258) q[2];
sx q[2];
rz(-1.8704869) q[2];
rz(3.0631915) q[3];
sx q[3];
rz(-1.6631118) q[3];
sx q[3];
rz(-1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8909797) q[0];
sx q[0];
rz(-1.5761292) q[0];
sx q[0];
rz(2.4839731) q[0];
rz(-1.744386) q[1];
sx q[1];
rz(-1.1746635) q[1];
sx q[1];
rz(2.2479642) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69041396) q[0];
sx q[0];
rz(-2.4126841) q[0];
sx q[0];
rz(0.50509392) q[0];
rz(-2.8473179) q[2];
sx q[2];
rz(-2.0308959) q[2];
sx q[2];
rz(1.0690881) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.11941714) q[1];
sx q[1];
rz(-2.4319473) q[1];
sx q[1];
rz(-0.31125734) q[1];
rz(-pi) q[2];
rz(-3.0393533) q[3];
sx q[3];
rz(-2.3792301) q[3];
sx q[3];
rz(1.8765212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7509193) q[2];
sx q[2];
rz(-0.94736391) q[2];
sx q[2];
rz(2.612109) q[2];
rz(0.47618619) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(2.7990394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3787518) q[0];
sx q[0];
rz(-1.6289926) q[0];
sx q[0];
rz(2.0513127) q[0];
rz(3.0293363) q[1];
sx q[1];
rz(-1.1039762) q[1];
sx q[1];
rz(-1.9876678) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.740828) q[0];
sx q[0];
rz(-2.4908998) q[0];
sx q[0];
rz(-3.0853737) q[0];
rz(0.23837337) q[2];
sx q[2];
rz(-2.5775238) q[2];
sx q[2];
rz(-0.43018815) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9824144) q[1];
sx q[1];
rz(-2.3783861) q[1];
sx q[1];
rz(-0.089341954) q[1];
x q[2];
rz(2.2418601) q[3];
sx q[3];
rz(-1.2601488) q[3];
sx q[3];
rz(1.336738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.551679) q[2];
sx q[2];
rz(-2.7605197) q[2];
sx q[2];
rz(-0.38044688) q[2];
rz(1.1278661) q[3];
sx q[3];
rz(-1.3600072) q[3];
sx q[3];
rz(1.5765223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(-0.34148759) q[0];
sx q[0];
rz(-1.8999758) q[0];
sx q[0];
rz(2.5019116) q[0];
rz(1.9027963) q[1];
sx q[1];
rz(-1.7265373) q[1];
sx q[1];
rz(1.9715086) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78478783) q[0];
sx q[0];
rz(-2.6492282) q[0];
sx q[0];
rz(2.3193293) q[0];
x q[1];
rz(-1.3001928) q[2];
sx q[2];
rz(-0.81297183) q[2];
sx q[2];
rz(1.1635309) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.77800345) q[1];
sx q[1];
rz(-0.52007857) q[1];
sx q[1];
rz(0.37104721) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13799237) q[3];
sx q[3];
rz(-2.237461) q[3];
sx q[3];
rz(-0.43825144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8273948) q[2];
sx q[2];
rz(-2.4908227) q[2];
sx q[2];
rz(2.7588552) q[2];
rz(-0.9283723) q[3];
sx q[3];
rz(-1.9675156) q[3];
sx q[3];
rz(-2.476957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9160354) q[0];
sx q[0];
rz(-1.6397497) q[0];
sx q[0];
rz(0.25892192) q[0];
rz(-2.4312773) q[1];
sx q[1];
rz(-1.0537035) q[1];
sx q[1];
rz(-0.47992596) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23586789) q[0];
sx q[0];
rz(-0.94434443) q[0];
sx q[0];
rz(-2.7263374) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2357113) q[2];
sx q[2];
rz(-1.3833628) q[2];
sx q[2];
rz(-2.4311709) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5727947) q[1];
sx q[1];
rz(-0.63981445) q[1];
sx q[1];
rz(-0.6154284) q[1];
rz(-pi) q[2];
rz(-0.89703538) q[3];
sx q[3];
rz(-1.1437136) q[3];
sx q[3];
rz(2.2626812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.84247983) q[2];
sx q[2];
rz(-2.2128236) q[2];
sx q[2];
rz(1.3170362) q[2];
rz(1.2420098) q[3];
sx q[3];
rz(-0.95359355) q[3];
sx q[3];
rz(-0.73808134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
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
rz(-1.7469035) q[2];
sx q[2];
rz(-1.0955784) q[2];
sx q[2];
rz(-0.57217862) q[2];
rz(-0.27771523) q[3];
sx q[3];
rz(-1.3125827) q[3];
sx q[3];
rz(-1.7607341) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];