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
rz(2.2489927) q[1];
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
x q[1];
rz(-2.0775665) q[2];
sx q[2];
rz(-2.1760586) q[2];
sx q[2];
rz(1.2111226) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6751911) q[1];
sx q[1];
rz(-2.2474504) q[1];
sx q[1];
rz(-2.7387709) q[1];
rz(-pi) q[2];
rz(-3.0953232) q[3];
sx q[3];
rz(-0.75152961) q[3];
sx q[3];
rz(2.3627594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9154174) q[2];
sx q[2];
rz(-1.5390652) q[2];
sx q[2];
rz(1.9809451) q[2];
rz(2.9246269) q[3];
sx q[3];
rz(-2.6187077) q[3];
sx q[3];
rz(1.0552361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.083667) q[0];
sx q[0];
rz(-0.96494976) q[0];
sx q[0];
rz(2.5657186) q[0];
rz(-1.2469762) q[1];
sx q[1];
rz(-1.2966825) q[1];
sx q[1];
rz(1.974568) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.089433) q[0];
sx q[0];
rz(-2.8773327) q[0];
sx q[0];
rz(-1.3454076) q[0];
rz(-pi) q[1];
rz(-1.1738271) q[2];
sx q[2];
rz(-0.57710986) q[2];
sx q[2];
rz(1.7259665) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.963672) q[1];
sx q[1];
rz(-2.3753787) q[1];
sx q[1];
rz(-1.4818165) q[1];
rz(-pi) q[2];
rz(0.59252177) q[3];
sx q[3];
rz(-2.3791109) q[3];
sx q[3];
rz(0.034686397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1077659) q[2];
sx q[2];
rz(-1.9627389) q[2];
sx q[2];
rz(-0.21437422) q[2];
rz(0.073444627) q[3];
sx q[3];
rz(-0.44973222) q[3];
sx q[3];
rz(0.28081056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6630994) q[0];
sx q[0];
rz(-0.61613023) q[0];
sx q[0];
rz(-1.2269155) q[0];
rz(2.7413209) q[1];
sx q[1];
rz(-1.8881256) q[1];
sx q[1];
rz(1.0148369) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6944511) q[0];
sx q[0];
rz(-1.5146291) q[0];
sx q[0];
rz(-2.4980314) q[0];
x q[1];
rz(0.24200183) q[2];
sx q[2];
rz(-1.1761464) q[2];
sx q[2];
rz(-1.637527) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.240757) q[1];
sx q[1];
rz(-1.6512617) q[1];
sx q[1];
rz(-0.97297538) q[1];
rz(1.8869927) q[3];
sx q[3];
rz(-2.4305775) q[3];
sx q[3];
rz(-2.7544114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0959452) q[2];
sx q[2];
rz(-1.0568876) q[2];
sx q[2];
rz(-2.2423559) q[2];
rz(-2.4441161) q[3];
sx q[3];
rz(-1.2858425) q[3];
sx q[3];
rz(0.60788679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.4863481) q[0];
sx q[0];
rz(-2.0589893) q[0];
sx q[0];
rz(0.43193257) q[0];
rz(2.5090384) q[1];
sx q[1];
rz(-2.7245941) q[1];
sx q[1];
rz(-2.5057709) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0441372) q[0];
sx q[0];
rz(-2.6090528) q[0];
sx q[0];
rz(-1.2116648) q[0];
rz(2.9767838) q[2];
sx q[2];
rz(-0.85068446) q[2];
sx q[2];
rz(-1.009843) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1078474) q[1];
sx q[1];
rz(-0.81106942) q[1];
sx q[1];
rz(-1.3147522) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5084247) q[3];
sx q[3];
rz(-0.35053262) q[3];
sx q[3];
rz(1.7115953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9277966) q[2];
sx q[2];
rz(-2.9979604) q[2];
sx q[2];
rz(-2.4528743) q[2];
rz(2.8074746) q[3];
sx q[3];
rz(-1.170661) q[3];
sx q[3];
rz(-0.14373246) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9976945) q[0];
sx q[0];
rz(-0.69902885) q[0];
sx q[0];
rz(-2.0671663) q[0];
rz(0.74514666) q[1];
sx q[1];
rz(-1.6586168) q[1];
sx q[1];
rz(-2.863046) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2715831) q[0];
sx q[0];
rz(-2.2402813) q[0];
sx q[0];
rz(-0.55218009) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15567051) q[2];
sx q[2];
rz(-0.61741932) q[2];
sx q[2];
rz(1.9215259) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7601732) q[1];
sx q[1];
rz(-2.1790494) q[1];
sx q[1];
rz(2.6890523) q[1];
rz(2.3260818) q[3];
sx q[3];
rz(-1.9933356) q[3];
sx q[3];
rz(-0.58327196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6986065) q[2];
sx q[2];
rz(-1.3977945) q[2];
sx q[2];
rz(0.29423514) q[2];
rz(3.0596628) q[3];
sx q[3];
rz(-0.51920813) q[3];
sx q[3];
rz(3.1053655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1086403) q[0];
sx q[0];
rz(-0.8529129) q[0];
sx q[0];
rz(-3.1325353) q[0];
rz(0.63502216) q[1];
sx q[1];
rz(-2.451684) q[1];
sx q[1];
rz(-0.10805282) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2436284) q[0];
sx q[0];
rz(-1.122323) q[0];
sx q[0];
rz(-1.9237706) q[0];
rz(-pi) q[1];
rz(-1.4322386) q[2];
sx q[2];
rz(-1.9493305) q[2];
sx q[2];
rz(2.0356503) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5080155) q[1];
sx q[1];
rz(-2.0532236) q[1];
sx q[1];
rz(-2.0111994) q[1];
rz(1.1507387) q[3];
sx q[3];
rz(-0.58745158) q[3];
sx q[3];
rz(-2.1849039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.99284995) q[2];
sx q[2];
rz(-1.381258) q[2];
sx q[2];
rz(1.8704869) q[2];
rz(0.078401119) q[3];
sx q[3];
rz(-1.4784808) q[3];
sx q[3];
rz(-1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25061297) q[0];
sx q[0];
rz(-1.5654634) q[0];
sx q[0];
rz(-2.4839731) q[0];
rz(1.3972067) q[1];
sx q[1];
rz(-1.1746635) q[1];
sx q[1];
rz(2.2479642) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2715627) q[0];
sx q[0];
rz(-1.2426408) q[0];
sx q[0];
rz(0.66332711) q[0];
rz(-2.8473179) q[2];
sx q[2];
rz(-1.1106967) q[2];
sx q[2];
rz(-1.0690881) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9295846) q[1];
sx q[1];
rz(-1.3699023) q[1];
sx q[1];
rz(2.4561873) q[1];
x q[2];
rz(2.3818447) q[3];
sx q[3];
rz(-1.6413416) q[3];
sx q[3];
rz(-0.37978803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7509193) q[2];
sx q[2];
rz(-0.94736391) q[2];
sx q[2];
rz(-2.612109) q[2];
rz(-0.47618619) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(0.34255323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(-2.0513127) q[0];
rz(-0.11225637) q[1];
sx q[1];
rz(-2.0376164) q[1];
sx q[1];
rz(1.9876678) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40076462) q[0];
sx q[0];
rz(-0.65069288) q[0];
sx q[0];
rz(3.0853737) q[0];
rz(-pi) q[1];
rz(-1.7190785) q[2];
sx q[2];
rz(-1.0245171) q[2];
sx q[2];
rz(0.15020457) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9824144) q[1];
sx q[1];
rz(-0.7632066) q[1];
sx q[1];
rz(-3.0522507) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89973255) q[3];
sx q[3];
rz(-1.2601488) q[3];
sx q[3];
rz(-1.8048546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5899137) q[2];
sx q[2];
rz(-0.38107291) q[2];
sx q[2];
rz(2.7611458) q[2];
rz(-2.0137265) q[3];
sx q[3];
rz(-1.7815855) q[3];
sx q[3];
rz(1.5650704) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8001051) q[0];
sx q[0];
rz(-1.8999758) q[0];
sx q[0];
rz(2.5019116) q[0];
rz(1.9027963) q[1];
sx q[1];
rz(-1.7265373) q[1];
sx q[1];
rz(-1.170084) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4719452) q[0];
sx q[0];
rz(-1.8983316) q[0];
sx q[0];
rz(1.9452842) q[0];
x q[1];
rz(1.8413999) q[2];
sx q[2];
rz(-0.81297183) q[2];
sx q[2];
rz(1.1635309) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.46718405) q[1];
sx q[1];
rz(-1.3896175) q[1];
sx q[1];
rz(0.49023899) q[1];
rz(-0.89948489) q[3];
sx q[3];
rz(-1.6791108) q[3];
sx q[3];
rz(1.9233821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8273948) q[2];
sx q[2];
rz(-2.4908227) q[2];
sx q[2];
rz(-0.38273746) q[2];
rz(0.9283723) q[3];
sx q[3];
rz(-1.174077) q[3];
sx q[3];
rz(0.66463566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2255573) q[0];
sx q[0];
rz(-1.501843) q[0];
sx q[0];
rz(-0.25892192) q[0];
rz(2.4312773) q[1];
sx q[1];
rz(-1.0537035) q[1];
sx q[1];
rz(0.47992596) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23586789) q[0];
sx q[0];
rz(-0.94434443) q[0];
sx q[0];
rz(-2.7263374) q[0];
rz(1.2725699) q[2];
sx q[2];
rz(-0.68694653) q[2];
sx q[2];
rz(0.62703122) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.5727947) q[1];
sx q[1];
rz(-2.5017782) q[1];
sx q[1];
rz(2.5261643) q[1];
rz(-pi) q[2];
rz(0.89703538) q[3];
sx q[3];
rz(-1.1437136) q[3];
sx q[3];
rz(-2.2626812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.84247983) q[2];
sx q[2];
rz(-0.92876902) q[2];
sx q[2];
rz(-1.3170362) q[2];
rz(-1.2420098) q[3];
sx q[3];
rz(-2.1879991) q[3];
sx q[3];
rz(-0.73808134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15923545) q[0];
sx q[0];
rz(-0.032082162) q[0];
sx q[0];
rz(-1.4627009) q[0];
rz(0.97933979) q[1];
sx q[1];
rz(-1.0995438) q[1];
sx q[1];
rz(-0.88811036) q[1];
rz(0.48158823) q[2];
sx q[2];
rz(-1.4143741) q[2];
sx q[2];
rz(1.0798567) q[2];
rz(2.8638774) q[3];
sx q[3];
rz(-1.3125827) q[3];
sx q[3];
rz(-1.7607341) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];