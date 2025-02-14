OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4966709) q[0];
sx q[0];
rz(-2.6725197) q[0];
sx q[0];
rz(0.99035779) q[0];
rz(-1.1733836) q[1];
sx q[1];
rz(2.2819509) q[1];
sx q[1];
rz(11.717164) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25739057) q[0];
sx q[0];
rz(-1.5369731) q[0];
sx q[0];
rz(-3.1387563) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8013632) q[2];
sx q[2];
rz(-2.2417667) q[2];
sx q[2];
rz(-0.036341993) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9473331) q[1];
sx q[1];
rz(-0.68945706) q[1];
sx q[1];
rz(1.1442776) q[1];
x q[2];
rz(1.0484656) q[3];
sx q[3];
rz(-1.4712787) q[3];
sx q[3];
rz(-2.8985648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2363756) q[2];
sx q[2];
rz(-1.507501) q[2];
sx q[2];
rz(-0.53984731) q[2];
rz(-1.5652462) q[3];
sx q[3];
rz(-2.7326475) q[3];
sx q[3];
rz(1.4096155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0419615) q[0];
sx q[0];
rz(-0.67622447) q[0];
sx q[0];
rz(1.3270295) q[0];
rz(-2.3565893) q[1];
sx q[1];
rz(-1.4229341) q[1];
sx q[1];
rz(2.6410417) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0324361) q[0];
sx q[0];
rz(-1.5457499) q[0];
sx q[0];
rz(-2.5902201) q[0];
rz(-pi) q[1];
rz(-2.0346257) q[2];
sx q[2];
rz(-2.7358449) q[2];
sx q[2];
rz(1.7057989) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.079502843) q[1];
sx q[1];
rz(-0.70598733) q[1];
sx q[1];
rz(-1.0393618) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0359382) q[3];
sx q[3];
rz(-0.81722471) q[3];
sx q[3];
rz(1.5508088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0565722) q[2];
sx q[2];
rz(-2.4041921) q[2];
sx q[2];
rz(-0.48193398) q[2];
rz(0.36007544) q[3];
sx q[3];
rz(-1.1432546) q[3];
sx q[3];
rz(-0.20865194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
rz(2.8423186) q[0];
sx q[0];
rz(-2.0330918) q[0];
sx q[0];
rz(0.47766787) q[0];
rz(2.281669) q[1];
sx q[1];
rz(-2.0932902) q[1];
sx q[1];
rz(-0.28712505) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.053517) q[0];
sx q[0];
rz(-1.6968124) q[0];
sx q[0];
rz(-3.0749223) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7874009) q[2];
sx q[2];
rz(-1.2239211) q[2];
sx q[2];
rz(-1.5225449) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9707253) q[1];
sx q[1];
rz(-1.5620462) q[1];
sx q[1];
rz(-0.88717242) q[1];
rz(-pi) q[2];
rz(0.085137376) q[3];
sx q[3];
rz(-1.9508024) q[3];
sx q[3];
rz(0.054084965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3890248) q[2];
sx q[2];
rz(-1.3501046) q[2];
sx q[2];
rz(-0.019850578) q[2];
rz(0.94414532) q[3];
sx q[3];
rz(-2.4343334) q[3];
sx q[3];
rz(-0.39785644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1860564) q[0];
sx q[0];
rz(-2.5113386) q[0];
sx q[0];
rz(-1.3007042) q[0];
rz(-0.4862673) q[1];
sx q[1];
rz(-1.9673037) q[1];
sx q[1];
rz(-3.1062612) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0606421) q[0];
sx q[0];
rz(-1.544702) q[0];
sx q[0];
rz(1.9992725) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8043121) q[2];
sx q[2];
rz(-1.1832675) q[2];
sx q[2];
rz(-0.0079517297) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7481193) q[1];
sx q[1];
rz(-1.5474209) q[1];
sx q[1];
rz(1.5764025) q[1];
rz(2.1256623) q[3];
sx q[3];
rz(-1.0175704) q[3];
sx q[3];
rz(-0.31496668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.47762927) q[2];
sx q[2];
rz(-2.5783381) q[2];
sx q[2];
rz(-0.25838724) q[2];
rz(-2.431331) q[3];
sx q[3];
rz(-0.60451549) q[3];
sx q[3];
rz(-1.5593504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6579987) q[0];
sx q[0];
rz(-0.76376629) q[0];
sx q[0];
rz(2.4972231) q[0];
rz(0.98980347) q[1];
sx q[1];
rz(-2.3376696) q[1];
sx q[1];
rz(-1.1092626) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2671701) q[0];
sx q[0];
rz(-2.5760898) q[0];
sx q[0];
rz(1.9419975) q[0];
x q[1];
rz(0.039626683) q[2];
sx q[2];
rz(-1.0575231) q[2];
sx q[2];
rz(2.4921592) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.39711827) q[1];
sx q[1];
rz(-0.79192415) q[1];
sx q[1];
rz(0.20361118) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8568138) q[3];
sx q[3];
rz(-2.3030465) q[3];
sx q[3];
rz(0.58457182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.41205078) q[2];
sx q[2];
rz(-1.5645626) q[2];
sx q[2];
rz(2.5315419) q[2];
rz(-0.080502056) q[3];
sx q[3];
rz(-2.9952315) q[3];
sx q[3];
rz(-3.1350873) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7215111) q[0];
sx q[0];
rz(-0.93669909) q[0];
sx q[0];
rz(-1.6625846) q[0];
rz(-1.9526019) q[1];
sx q[1];
rz(-0.34591302) q[1];
sx q[1];
rz(-1.6993274) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95992559) q[0];
sx q[0];
rz(-0.64347351) q[0];
sx q[0];
rz(1.305278) q[0];
rz(-0.10043721) q[2];
sx q[2];
rz(-0.58375508) q[2];
sx q[2];
rz(-0.69086087) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.14076041) q[1];
sx q[1];
rz(-0.78703431) q[1];
sx q[1];
rz(-1.0972904) q[1];
rz(-1.9154432) q[3];
sx q[3];
rz(-2.0807869) q[3];
sx q[3];
rz(0.6450212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6815765) q[2];
sx q[2];
rz(-2.4676393) q[2];
sx q[2];
rz(-0.67503929) q[2];
rz(-2.8071844) q[3];
sx q[3];
rz(-1.4760009) q[3];
sx q[3];
rz(2.7736751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48208958) q[0];
sx q[0];
rz(-0.98564321) q[0];
sx q[0];
rz(3.0159045) q[0];
rz(-2.9478759) q[1];
sx q[1];
rz(-0.88231641) q[1];
sx q[1];
rz(1.151459) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73456681) q[0];
sx q[0];
rz(-2.1687963) q[0];
sx q[0];
rz(2.1728188) q[0];
x q[1];
rz(1.2000243) q[2];
sx q[2];
rz(-1.3091012) q[2];
sx q[2];
rz(-0.38693025) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1035422) q[1];
sx q[1];
rz(-1.0443496) q[1];
sx q[1];
rz(-1.0564694) q[1];
x q[2];
rz(-1.098387) q[3];
sx q[3];
rz(-2.189872) q[3];
sx q[3];
rz(-0.38451871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8573528) q[2];
sx q[2];
rz(-1.3606631) q[2];
sx q[2];
rz(2.8947158) q[2];
rz(0.85865584) q[3];
sx q[3];
rz(-2.1571428) q[3];
sx q[3];
rz(-2.8624559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7849279) q[0];
sx q[0];
rz(-0.44342884) q[0];
sx q[0];
rz(-0.32522935) q[0];
rz(-0.52109703) q[1];
sx q[1];
rz(-2.5083713) q[1];
sx q[1];
rz(0.31347832) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8756008) q[0];
sx q[0];
rz(-0.82339215) q[0];
sx q[0];
rz(2.4317047) q[0];
rz(-pi) q[1];
rz(1.8584537) q[2];
sx q[2];
rz(-2.9069942) q[2];
sx q[2];
rz(-1.2755659) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9814138) q[1];
sx q[1];
rz(-1.2048843) q[1];
sx q[1];
rz(-1.2104079) q[1];
rz(-pi) q[2];
rz(0.20200396) q[3];
sx q[3];
rz(-1.0583377) q[3];
sx q[3];
rz(-0.92512586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7079805) q[2];
sx q[2];
rz(-1.2182451) q[2];
sx q[2];
rz(-1.8617967) q[2];
rz(-0.72569877) q[3];
sx q[3];
rz(-1.1596707) q[3];
sx q[3];
rz(2.3316135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63876605) q[0];
sx q[0];
rz(-1.1932729) q[0];
sx q[0];
rz(2.9916812) q[0];
rz(1.2431078) q[1];
sx q[1];
rz(-1.5538235) q[1];
sx q[1];
rz(-1.4814203) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96326665) q[0];
sx q[0];
rz(-1.5218922) q[0];
sx q[0];
rz(0.045128926) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3301351) q[2];
sx q[2];
rz(-1.6208555) q[2];
sx q[2];
rz(-0.71739774) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8259623) q[1];
sx q[1];
rz(-0.44188979) q[1];
sx q[1];
rz(0.52200861) q[1];
rz(-1.0970988) q[3];
sx q[3];
rz(-0.27333958) q[3];
sx q[3];
rz(-0.36318446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.82568613) q[2];
sx q[2];
rz(-1.76182) q[2];
sx q[2];
rz(0.97563499) q[2];
rz(-1.1286831) q[3];
sx q[3];
rz(-0.55207878) q[3];
sx q[3];
rz(-0.25477195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.16633701) q[0];
sx q[0];
rz(-1.0003426) q[0];
sx q[0];
rz(2.9794203) q[0];
rz(1.3316679) q[1];
sx q[1];
rz(-0.63260308) q[1];
sx q[1];
rz(0.046028927) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9105658) q[0];
sx q[0];
rz(-1.5767158) q[0];
sx q[0];
rz(-1.5529412) q[0];
rz(-pi) q[1];
rz(-0.68387939) q[2];
sx q[2];
rz(-2.5683142) q[2];
sx q[2];
rz(-1.8901907) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.11124736) q[1];
sx q[1];
rz(-1.971389) q[1];
sx q[1];
rz(-2.9811893) q[1];
rz(0.45451136) q[3];
sx q[3];
rz(-0.3743518) q[3];
sx q[3];
rz(0.020937048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.33596805) q[2];
sx q[2];
rz(-0.38463548) q[2];
sx q[2];
rz(1.9920721) q[2];
rz(2.6286821) q[3];
sx q[3];
rz(-1.3866813) q[3];
sx q[3];
rz(0.30470595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49160663) q[0];
sx q[0];
rz(-0.91616022) q[0];
sx q[0];
rz(1.1471163) q[0];
rz(-0.06123771) q[1];
sx q[1];
rz(-1.1920659) q[1];
sx q[1];
rz(1.7658284) q[1];
rz(2.3394924) q[2];
sx q[2];
rz(-1.987793) q[2];
sx q[2];
rz(-1.1588617) q[2];
rz(0.41027222) q[3];
sx q[3];
rz(-0.82220746) q[3];
sx q[3];
rz(-2.0668277) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
