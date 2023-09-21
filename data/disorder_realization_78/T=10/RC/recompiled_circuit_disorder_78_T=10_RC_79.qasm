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
rz(-2.5506033) q[0];
sx q[0];
rz(-0.58340573) q[0];
rz(2.9572339) q[1];
sx q[1];
rz(-0.98362041) q[1];
sx q[1];
rz(-0.89259994) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6658202) q[0];
sx q[0];
rz(-1.2736397) q[0];
sx q[0];
rz(-0.99065234) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61172723) q[2];
sx q[2];
rz(-0.76843843) q[2];
sx q[2];
rz(-0.4380463) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6751911) q[1];
sx q[1];
rz(-2.2474504) q[1];
sx q[1];
rz(2.7387709) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3905972) q[3];
sx q[3];
rz(-1.5392116) q[3];
sx q[3];
rz(-0.75814523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2261752) q[2];
sx q[2];
rz(-1.5390652) q[2];
sx q[2];
rz(-1.9809451) q[2];
rz(-2.9246269) q[3];
sx q[3];
rz(-2.6187077) q[3];
sx q[3];
rz(2.0863566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3008227) q[0];
sx q[0];
rz(-1.6292028) q[0];
sx q[0];
rz(1.8286684) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1115233) q[2];
sx q[2];
rz(-1.3582555) q[2];
sx q[2];
rz(2.9589047) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.963672) q[1];
sx q[1];
rz(-2.3753787) q[1];
sx q[1];
rz(-1.4818165) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67006536) q[3];
sx q[3];
rz(-1.1747922) q[3];
sx q[3];
rz(2.0585287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4784933) q[0];
sx q[0];
rz(-0.61613023) q[0];
sx q[0];
rz(1.9146772) q[0];
rz(-2.7413209) q[1];
sx q[1];
rz(-1.2534671) q[1];
sx q[1];
rz(-2.1267557) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0815711) q[0];
sx q[0];
rz(-2.2131753) q[0];
sx q[0];
rz(1.6409671) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8995908) q[2];
sx q[2];
rz(-1.1761464) q[2];
sx q[2];
rz(-1.5040656) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9008357) q[1];
sx q[1];
rz(-1.6512617) q[1];
sx q[1];
rz(0.97297538) q[1];
rz(-pi) q[2];
x q[2];
rz(0.26168163) q[3];
sx q[3];
rz(-0.90173429) q[3];
sx q[3];
rz(-2.3467968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0456475) q[2];
sx q[2];
rz(-1.0568876) q[2];
sx q[2];
rz(0.8992368) q[2];
rz(0.69747654) q[3];
sx q[3];
rz(-1.2858425) q[3];
sx q[3];
rz(-2.5337059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4863481) q[0];
sx q[0];
rz(-1.0826033) q[0];
sx q[0];
rz(-2.7096601) q[0];
rz(0.63255429) q[1];
sx q[1];
rz(-2.7245941) q[1];
sx q[1];
rz(-0.63582173) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0441372) q[0];
sx q[0];
rz(-0.53253981) q[0];
sx q[0];
rz(1.9299279) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7556778) q[2];
sx q[2];
rz(-2.4061678) q[2];
sx q[2];
rz(2.3787969) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7446049) q[1];
sx q[1];
rz(-0.79345353) q[1];
sx q[1];
rz(2.8810487) q[1];
rz(0.28663978) q[3];
sx q[3];
rz(-1.3661824) q[3];
sx q[3];
rz(-0.74433792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2137961) q[2];
sx q[2];
rz(-2.9979604) q[2];
sx q[2];
rz(0.68871838) q[2];
rz(-2.8074746) q[3];
sx q[3];
rz(-1.170661) q[3];
sx q[3];
rz(0.14373246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14389811) q[0];
sx q[0];
rz(-0.69902885) q[0];
sx q[0];
rz(1.0744263) q[0];
rz(-2.396446) q[1];
sx q[1];
rz(-1.6586168) q[1];
sx q[1];
rz(-2.863046) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8700096) q[0];
sx q[0];
rz(-2.2402813) q[0];
sx q[0];
rz(2.5894126) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61168806) q[2];
sx q[2];
rz(-1.4809161) q[2];
sx q[2];
rz(-2.9181366) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6812233) q[1];
sx q[1];
rz(-1.2037828) q[1];
sx q[1];
rz(-2.2296434) q[1];
x q[2];
rz(-0.55322247) q[3];
sx q[3];
rz(-0.8953989) q[3];
sx q[3];
rz(1.7856057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6986065) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1086403) q[0];
sx q[0];
rz(-2.2886798) q[0];
sx q[0];
rz(-3.1325353) q[0];
rz(0.63502216) q[1];
sx q[1];
rz(-2.451684) q[1];
sx q[1];
rz(3.0335398) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51443433) q[0];
sx q[0];
rz(-1.8875727) q[0];
sx q[0];
rz(-0.47382521) q[0];
rz(-pi) q[1];
rz(0.38185264) q[2];
sx q[2];
rz(-1.4421041) q[2];
sx q[2];
rz(0.51634386) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9891226) q[1];
sx q[1];
rz(-1.9580541) q[1];
sx q[1];
rz(0.5247922) q[1];
rz(-pi) q[2];
rz(1.990854) q[3];
sx q[3];
rz(-2.5541411) q[3];
sx q[3];
rz(-2.1849039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.99284995) q[2];
sx q[2];
rz(-1.7603346) q[2];
sx q[2];
rz(1.2711058) q[2];
rz(-3.0631915) q[3];
sx q[3];
rz(-1.6631118) q[3];
sx q[3];
rz(1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25061297) q[0];
sx q[0];
rz(-1.5654634) q[0];
sx q[0];
rz(0.65761956) q[0];
rz(-1.3972067) q[1];
sx q[1];
rz(-1.9669292) q[1];
sx q[1];
rz(-0.89362842) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.87003) q[0];
sx q[0];
rz(-1.2426408) q[0];
sx q[0];
rz(0.66332711) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29427476) q[2];
sx q[2];
rz(-1.1106967) q[2];
sx q[2];
rz(2.0725046) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0221755) q[1];
sx q[1];
rz(-0.70964538) q[1];
sx q[1];
rz(0.31125734) q[1];
rz(-3.0393533) q[3];
sx q[3];
rz(-2.3792301) q[3];
sx q[3];
rz(1.8765212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7509193) q[2];
sx q[2];
rz(-2.1942287) q[2];
sx q[2];
rz(2.612109) q[2];
rz(0.47618619) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(-0.34255323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.3787518) q[0];
sx q[0];
rz(-1.5126001) q[0];
sx q[0];
rz(1.09028) q[0];
rz(3.0293363) q[1];
sx q[1];
rz(-2.0376164) q[1];
sx q[1];
rz(-1.1539248) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2147804) q[0];
sx q[0];
rz(-1.6048389) q[0];
sx q[0];
rz(0.649931) q[0];
rz(-2.5904028) q[2];
sx q[2];
rz(-1.4442208) q[2];
sx q[2];
rz(-1.3431431) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1591783) q[1];
sx q[1];
rz(-0.7632066) q[1];
sx q[1];
rz(0.089341954) q[1];
rz(2.7525547) q[3];
sx q[3];
rz(-2.2044047) q[3];
sx q[3];
rz(-2.6694359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5899137) q[2];
sx q[2];
rz(-2.7605197) q[2];
sx q[2];
rz(0.38044688) q[2];
rz(-1.1278661) q[3];
sx q[3];
rz(-1.3600072) q[3];
sx q[3];
rz(-1.5765223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.34148759) q[0];
sx q[0];
rz(-1.8999758) q[0];
sx q[0];
rz(2.5019116) q[0];
rz(-1.9027963) q[1];
sx q[1];
rz(-1.4150554) q[1];
sx q[1];
rz(-1.170084) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4719452) q[0];
sx q[0];
rz(-1.8983316) q[0];
sx q[0];
rz(-1.1963084) q[0];
rz(-0.77634546) q[2];
sx q[2];
rz(-1.3753969) q[2];
sx q[2];
rz(2.5459144) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.46718405) q[1];
sx q[1];
rz(-1.3896175) q[1];
sx q[1];
rz(-0.49023899) q[1];
x q[2];
rz(-0.89948489) q[3];
sx q[3];
rz(-1.6791108) q[3];
sx q[3];
rz(1.9233821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3141979) q[2];
sx q[2];
rz(-2.4908227) q[2];
sx q[2];
rz(-2.7588552) q[2];
rz(0.9283723) q[3];
sx q[3];
rz(-1.174077) q[3];
sx q[3];
rz(0.66463566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9160354) q[0];
sx q[0];
rz(-1.6397497) q[0];
sx q[0];
rz(-2.8826707) q[0];
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
rz(-2.9057248) q[0];
sx q[0];
rz(-0.94434443) q[0];
sx q[0];
rz(0.41525526) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23649044) q[2];
sx q[2];
rz(-0.91954008) q[2];
sx q[2];
rz(-2.1361534) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5727947) q[1];
sx q[1];
rz(-0.63981445) q[1];
sx q[1];
rz(0.6154284) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2445573) q[3];
sx q[3];
rz(-1.1437136) q[3];
sx q[3];
rz(0.87891146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.84247983) q[2];
sx q[2];
rz(-0.92876902) q[2];
sx q[2];
rz(1.3170362) q[2];
rz(-1.2420098) q[3];
sx q[3];
rz(-2.1879991) q[3];
sx q[3];
rz(2.4035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9823572) q[0];
sx q[0];
rz(-3.1095105) q[0];
sx q[0];
rz(1.6788917) q[0];
rz(-2.1622529) q[1];
sx q[1];
rz(-1.0995438) q[1];
sx q[1];
rz(-0.88811036) q[1];
rz(-2.813415) q[2];
sx q[2];
rz(-0.50445088) q[2];
sx q[2];
rz(2.9403461) q[2];
rz(0.76673037) q[3];
sx q[3];
rz(-2.7646716) q[3];
sx q[3];
rz(-0.9203831) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];