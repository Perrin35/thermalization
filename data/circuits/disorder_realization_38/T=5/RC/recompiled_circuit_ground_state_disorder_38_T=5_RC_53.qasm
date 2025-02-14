OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.225086) q[0];
sx q[0];
rz(-0.14296159) q[0];
sx q[0];
rz(11.077865) q[0];
rz(1.1324963) q[1];
sx q[1];
rz(-2.7096665) q[1];
sx q[1];
rz(-0.63634029) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4637488) q[0];
sx q[0];
rz(-1.1882458) q[0];
sx q[0];
rz(2.2027459) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77677988) q[2];
sx q[2];
rz(-2.619639) q[2];
sx q[2];
rz(1.7372434) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8752567) q[1];
sx q[1];
rz(-2.7260029) q[1];
sx q[1];
rz(-1.8824091) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79145517) q[3];
sx q[3];
rz(-0.88578445) q[3];
sx q[3];
rz(-0.58429694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1119614) q[2];
sx q[2];
rz(-2.0165069) q[2];
sx q[2];
rz(-2.4911528) q[2];
rz(1.8247617) q[3];
sx q[3];
rz(-1.3464758) q[3];
sx q[3];
rz(-3.0397547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0874852) q[0];
sx q[0];
rz(-2.5547145) q[0];
sx q[0];
rz(2.6869539) q[0];
rz(-3.1184323) q[1];
sx q[1];
rz(-1.6975479) q[1];
sx q[1];
rz(0.32618162) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21837337) q[0];
sx q[0];
rz(-1.0207286) q[0];
sx q[0];
rz(-2.0327507) q[0];
x q[1];
rz(-2.1670879) q[2];
sx q[2];
rz(-2.6428004) q[2];
sx q[2];
rz(1.6287664) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5311969) q[1];
sx q[1];
rz(-2.463974) q[1];
sx q[1];
rz(-0.68790959) q[1];
rz(-pi) q[2];
x q[2];
rz(1.420191) q[3];
sx q[3];
rz(-2.4753627) q[3];
sx q[3];
rz(-3.1325983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0252016) q[2];
sx q[2];
rz(-1.9390743) q[2];
sx q[2];
rz(0.95428673) q[2];
rz(1.8918234) q[3];
sx q[3];
rz(-0.3229177) q[3];
sx q[3];
rz(0.0013110411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8702451) q[0];
sx q[0];
rz(-1.988669) q[0];
sx q[0];
rz(1.1767607) q[0];
rz(2.7765043) q[1];
sx q[1];
rz(-1.8389713) q[1];
sx q[1];
rz(-1.6129859) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9823526) q[0];
sx q[0];
rz(-1.5499347) q[0];
sx q[0];
rz(0.013289159) q[0];
x q[1];
rz(3.0625383) q[2];
sx q[2];
rz(-2.0131265) q[2];
sx q[2];
rz(-2.4238264) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6512201) q[1];
sx q[1];
rz(-1.4191322) q[1];
sx q[1];
rz(3.1107799) q[1];
x q[2];
rz(-3.0729483) q[3];
sx q[3];
rz(-1.513283) q[3];
sx q[3];
rz(2.7781093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.24885808) q[2];
sx q[2];
rz(-2.5881519) q[2];
sx q[2];
rz(1.9888606) q[2];
rz(1.8966127) q[3];
sx q[3];
rz(-2.1608519) q[3];
sx q[3];
rz(2.8583756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1588441) q[0];
sx q[0];
rz(-pi/12) q[0];
sx q[0];
rz(-0.68956476) q[0];
rz(-0.025029643) q[1];
sx q[1];
rz(-1.6233147) q[1];
sx q[1];
rz(-1.4208581) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0293504) q[0];
sx q[0];
rz(-0.94932244) q[0];
sx q[0];
rz(1.6981324) q[0];
rz(0.17949149) q[2];
sx q[2];
rz(-1.9948261) q[2];
sx q[2];
rz(-1.9551203) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5031227) q[1];
sx q[1];
rz(-1.5431738) q[1];
sx q[1];
rz(-0.50325985) q[1];
x q[2];
rz(-2.8042364) q[3];
sx q[3];
rz(-1.8280262) q[3];
sx q[3];
rz(1.4816062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20874061) q[2];
sx q[2];
rz(-2.9661861) q[2];
sx q[2];
rz(1.1669195) q[2];
rz(2.6973727) q[3];
sx q[3];
rz(-1.9053713) q[3];
sx q[3];
rz(-0.3096295) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8532448) q[0];
sx q[0];
rz(-1.7720368) q[0];
sx q[0];
rz(-0.41611588) q[0];
rz(-0.5101282) q[1];
sx q[1];
rz(-0.77113873) q[1];
sx q[1];
rz(-0.77521926) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3907523) q[0];
sx q[0];
rz(-1.3715944) q[0];
sx q[0];
rz(-2.5222523) q[0];
rz(-pi) q[1];
rz(1.3890319) q[2];
sx q[2];
rz(-1.1715638) q[2];
sx q[2];
rz(2.4295074) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.21501274) q[1];
sx q[1];
rz(-1.4583419) q[1];
sx q[1];
rz(1.6968898) q[1];
rz(-2.5465762) q[3];
sx q[3];
rz(-1.2538246) q[3];
sx q[3];
rz(1.9645312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9267209) q[2];
sx q[2];
rz(-2.1472609) q[2];
sx q[2];
rz(1.035824) q[2];
rz(0.18946798) q[3];
sx q[3];
rz(-0.47179705) q[3];
sx q[3];
rz(-2.6389879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8871317) q[0];
sx q[0];
rz(-3.0939565) q[0];
sx q[0];
rz(2.7857842) q[0];
rz(1.4954781) q[1];
sx q[1];
rz(-0.9551841) q[1];
sx q[1];
rz(-1.1411508) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4116652) q[0];
sx q[0];
rz(-0.91331867) q[0];
sx q[0];
rz(2.8748926) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95924033) q[2];
sx q[2];
rz(-1.578176) q[2];
sx q[2];
rz(0.5317229) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0867978) q[1];
sx q[1];
rz(-1.7112477) q[1];
sx q[1];
rz(-1.4649006) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5445956) q[3];
sx q[3];
rz(-2.7075534) q[3];
sx q[3];
rz(-1.7569913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5100539) q[2];
sx q[2];
rz(-2.4271991) q[2];
sx q[2];
rz(1.4368524) q[2];
rz(2.0884183) q[3];
sx q[3];
rz(-1.6097693) q[3];
sx q[3];
rz(0.50301445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.961504) q[0];
sx q[0];
rz(-2.8453974) q[0];
sx q[0];
rz(3.0190813) q[0];
rz(-2.4854614) q[1];
sx q[1];
rz(-1.2686814) q[1];
sx q[1];
rz(1.8001385) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91514523) q[0];
sx q[0];
rz(-1.6065803) q[0];
sx q[0];
rz(0.045374845) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35690633) q[2];
sx q[2];
rz(-0.30115899) q[2];
sx q[2];
rz(-0.44277175) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0813876) q[1];
sx q[1];
rz(-1.837018) q[1];
sx q[1];
rz(-0.93764825) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44388598) q[3];
sx q[3];
rz(-1.6735014) q[3];
sx q[3];
rz(0.041180276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5682257) q[2];
sx q[2];
rz(-1.919701) q[2];
sx q[2];
rz(-1.6778256) q[2];
rz(-0.73417869) q[3];
sx q[3];
rz(-2.8575183) q[3];
sx q[3];
rz(-1.3045788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3463335) q[0];
sx q[0];
rz(-1.7743552) q[0];
sx q[0];
rz(2.6829868) q[0];
rz(-2.1268225) q[1];
sx q[1];
rz(-1.9472803) q[1];
sx q[1];
rz(-1.0626622) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5073481) q[0];
sx q[0];
rz(-1.5131823) q[0];
sx q[0];
rz(0.16648023) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1193211) q[2];
sx q[2];
rz(-2.3000225) q[2];
sx q[2];
rz(-2.4265576) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1717689) q[1];
sx q[1];
rz(-2.9967628) q[1];
sx q[1];
rz(-2.0098196) q[1];
rz(-pi) q[2];
rz(0.32286947) q[3];
sx q[3];
rz(-1.0738651) q[3];
sx q[3];
rz(0.54847713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8899272) q[2];
sx q[2];
rz(-1.7532316) q[2];
sx q[2];
rz(-1.3512705) q[2];
rz(-1.2190367) q[3];
sx q[3];
rz(-0.42870298) q[3];
sx q[3];
rz(-0.8333227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625921) q[0];
sx q[0];
rz(-1.7575678) q[0];
sx q[0];
rz(-0.90743995) q[0];
rz(-0.22706789) q[1];
sx q[1];
rz(-2.0957969) q[1];
sx q[1];
rz(-0.89920941) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.718841) q[0];
sx q[0];
rz(-0.97373527) q[0];
sx q[0];
rz(0.62946749) q[0];
rz(-1.4883409) q[2];
sx q[2];
rz(-1.4926148) q[2];
sx q[2];
rz(2.9417335) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.21630741) q[1];
sx q[1];
rz(-1.3853964) q[1];
sx q[1];
rz(2.0865296) q[1];
x q[2];
rz(0.050961625) q[3];
sx q[3];
rz(-2.3182541) q[3];
sx q[3];
rz(-2.2281102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1850618) q[2];
sx q[2];
rz(-1.9756292) q[2];
sx q[2];
rz(-0.60066191) q[2];
rz(2.0447842) q[3];
sx q[3];
rz(-0.59022248) q[3];
sx q[3];
rz(-2.2685331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3908865) q[0];
sx q[0];
rz(-1.0989256) q[0];
sx q[0];
rz(0.98439687) q[0];
rz(-0.4920494) q[1];
sx q[1];
rz(-1.4130519) q[1];
sx q[1];
rz(1.1955998) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0988783) q[0];
sx q[0];
rz(-0.98473583) q[0];
sx q[0];
rz(1.1436966) q[0];
rz(-pi) q[1];
rz(2.1629647) q[2];
sx q[2];
rz(-0.76972967) q[2];
sx q[2];
rz(-0.20537381) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3477563) q[1];
sx q[1];
rz(-1.6518005) q[1];
sx q[1];
rz(-1.9543314) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77287425) q[3];
sx q[3];
rz(-1.2307271) q[3];
sx q[3];
rz(-0.085207663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.055858) q[2];
sx q[2];
rz(-1.6207638) q[2];
sx q[2];
rz(-2.1585507) q[2];
rz(2.8899657) q[3];
sx q[3];
rz(-2.7139137) q[3];
sx q[3];
rz(-1.0035286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6428103) q[0];
sx q[0];
rz(-0.97220535) q[0];
sx q[0];
rz(-2.4844949) q[0];
rz(-1.2596754) q[1];
sx q[1];
rz(-1.7301662) q[1];
sx q[1];
rz(0.87283254) q[1];
rz(2.2839584) q[2];
sx q[2];
rz(-1.7339755) q[2];
sx q[2];
rz(-0.82790464) q[2];
rz(1.891741) q[3];
sx q[3];
rz(-1.8886011) q[3];
sx q[3];
rz(-0.26816751) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
