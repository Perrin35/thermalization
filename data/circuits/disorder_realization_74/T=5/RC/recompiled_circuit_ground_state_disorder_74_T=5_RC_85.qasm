OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.57269078) q[0];
sx q[0];
rz(-2.1532018) q[0];
sx q[0];
rz(0.44901499) q[0];
rz(-0.16945893) q[1];
sx q[1];
rz(-3.0100477) q[1];
sx q[1];
rz(-2.0102672) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5857475) q[0];
sx q[0];
rz(-0.97363421) q[0];
sx q[0];
rz(-2.7859429) q[0];
x q[1];
rz(0.63458459) q[2];
sx q[2];
rz(-1.1401083) q[2];
sx q[2];
rz(2.787295) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0180964) q[1];
sx q[1];
rz(-0.61401788) q[1];
sx q[1];
rz(1.2662925) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9692124) q[3];
sx q[3];
rz(-2.1076492) q[3];
sx q[3];
rz(-2.7011724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.955287) q[2];
sx q[2];
rz(-2.7489642) q[2];
sx q[2];
rz(3.0379831) q[2];
rz(-0.3668395) q[3];
sx q[3];
rz(-1.6678526) q[3];
sx q[3];
rz(1.8240671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1845301) q[0];
sx q[0];
rz(-0.4758895) q[0];
sx q[0];
rz(2.6153508) q[0];
rz(1.8980252) q[1];
sx q[1];
rz(-1.4105816) q[1];
sx q[1];
rz(-0.19613656) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2945902) q[0];
sx q[0];
rz(-0.46200141) q[0];
sx q[0];
rz(-1.7209956) q[0];
x q[1];
rz(-2.9483825) q[2];
sx q[2];
rz(-0.71882283) q[2];
sx q[2];
rz(-0.012398243) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.96207075) q[1];
sx q[1];
rz(-1.7129494) q[1];
sx q[1];
rz(2.0217871) q[1];
x q[2];
rz(1.5603746) q[3];
sx q[3];
rz(-1.5589514) q[3];
sx q[3];
rz(-2.7238977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5474995) q[2];
sx q[2];
rz(-0.27442351) q[2];
sx q[2];
rz(0.37626949) q[2];
rz(-0.88092342) q[3];
sx q[3];
rz(-1.3353525) q[3];
sx q[3];
rz(-2.3295565) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6905717) q[0];
sx q[0];
rz(-0.32830992) q[0];
sx q[0];
rz(2.1300533) q[0];
rz(-0.66863376) q[1];
sx q[1];
rz(-1.9790383) q[1];
sx q[1];
rz(2.5943601) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6100734) q[0];
sx q[0];
rz(-1.7721869) q[0];
sx q[0];
rz(1.5511484) q[0];
x q[1];
rz(1.2195361) q[2];
sx q[2];
rz(-0.81163663) q[2];
sx q[2];
rz(-1.3724302) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.95067945) q[1];
sx q[1];
rz(-0.70531323) q[1];
sx q[1];
rz(-1.9478134) q[1];
rz(-2.3272334) q[3];
sx q[3];
rz(-1.7118771) q[3];
sx q[3];
rz(0.94179487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4554567) q[2];
sx q[2];
rz(-0.28527173) q[2];
sx q[2];
rz(-1.5020465) q[2];
rz(-0.57602588) q[3];
sx q[3];
rz(-1.4833769) q[3];
sx q[3];
rz(2.7801133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5107002) q[0];
sx q[0];
rz(-0.80131131) q[0];
sx q[0];
rz(-0.33962387) q[0];
rz(-0.54061186) q[1];
sx q[1];
rz(-0.70053354) q[1];
sx q[1];
rz(0.9300173) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.697111) q[0];
sx q[0];
rz(-0.15938317) q[0];
sx q[0];
rz(-1.1109933) q[0];
x q[1];
rz(2.5581237) q[2];
sx q[2];
rz(-1.2440878) q[2];
sx q[2];
rz(2.938156) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8992638) q[1];
sx q[1];
rz(-1.5643969) q[1];
sx q[1];
rz(1.0914299) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2131683) q[3];
sx q[3];
rz(-1.9898212) q[3];
sx q[3];
rz(0.85538543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0129619) q[2];
sx q[2];
rz(-1.7065115) q[2];
sx q[2];
rz(-1.5738515) q[2];
rz(-0.74357998) q[3];
sx q[3];
rz(-0.79135528) q[3];
sx q[3];
rz(-2.1775406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6745233) q[0];
sx q[0];
rz(-2.2860797) q[0];
sx q[0];
rz(-0.056644406) q[0];
rz(1.665834) q[1];
sx q[1];
rz(-2.00878) q[1];
sx q[1];
rz(-1.7830361) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9525653) q[0];
sx q[0];
rz(-1.0882821) q[0];
sx q[0];
rz(-2.5607462) q[0];
rz(-pi) q[1];
rz(1.3681202) q[2];
sx q[2];
rz(-0.92331013) q[2];
sx q[2];
rz(1.7062239) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.1093692) q[1];
sx q[1];
rz(-1.0894766) q[1];
sx q[1];
rz(0.66415031) q[1];
x q[2];
rz(-0.18601619) q[3];
sx q[3];
rz(-0.56944344) q[3];
sx q[3];
rz(-1.525477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8654827) q[2];
sx q[2];
rz(-1.411974) q[2];
sx q[2];
rz(1.126368) q[2];
rz(-1.8051091) q[3];
sx q[3];
rz(-2.0650605) q[3];
sx q[3];
rz(-2.0520463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.69018501) q[0];
sx q[0];
rz(-2.1162338) q[0];
sx q[0];
rz(-3.0080646) q[0];
rz(0.97767699) q[1];
sx q[1];
rz(-0.97594273) q[1];
sx q[1];
rz(-0.28688637) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27659076) q[0];
sx q[0];
rz(-1.6185068) q[0];
sx q[0];
rz(-2.2722831) q[0];
x q[1];
rz(-3.098549) q[2];
sx q[2];
rz(-1.2428987) q[2];
sx q[2];
rz(0.69059935) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.15039728) q[1];
sx q[1];
rz(-2.5270473) q[1];
sx q[1];
rz(-1.9646364) q[1];
rz(-pi) q[2];
x q[2];
rz(1.828015) q[3];
sx q[3];
rz(-2.1495499) q[3];
sx q[3];
rz(2.677315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7202683) q[2];
sx q[2];
rz(-1.6608394) q[2];
sx q[2];
rz(-1.6555017) q[2];
rz(-2.7739575) q[3];
sx q[3];
rz(-1.0272107) q[3];
sx q[3];
rz(2.0462842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4246282) q[0];
sx q[0];
rz(-2.3724738) q[0];
sx q[0];
rz(-0.47384438) q[0];
rz(-2.4773856) q[1];
sx q[1];
rz(-1.217548) q[1];
sx q[1];
rz(1.7582105) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5067399) q[0];
sx q[0];
rz(-2.228745) q[0];
sx q[0];
rz(0.94292504) q[0];
rz(-pi) q[1];
rz(-1.9888617) q[2];
sx q[2];
rz(-1.5961313) q[2];
sx q[2];
rz(2.096614) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.7935087) q[1];
sx q[1];
rz(-1.7941107) q[1];
sx q[1];
rz(-2.9053754) q[1];
rz(2.9630205) q[3];
sx q[3];
rz(-1.7925486) q[3];
sx q[3];
rz(-0.54820433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.15709269) q[2];
sx q[2];
rz(-1.7326771) q[2];
sx q[2];
rz(0.3248997) q[2];
rz(-2.7609008) q[3];
sx q[3];
rz(-0.84563962) q[3];
sx q[3];
rz(-2.0438173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070504524) q[0];
sx q[0];
rz(-0.66362137) q[0];
sx q[0];
rz(2.2027503) q[0];
rz(2.4414869) q[1];
sx q[1];
rz(-0.63799262) q[1];
sx q[1];
rz(-0.28900388) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5225713) q[0];
sx q[0];
rz(-2.8998525) q[0];
sx q[0];
rz(-1.2416583) q[0];
rz(-pi) q[1];
rz(-2.2174066) q[2];
sx q[2];
rz(-1.1614387) q[2];
sx q[2];
rz(1.5034624) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6470312) q[1];
sx q[1];
rz(-0.53003487) q[1];
sx q[1];
rz(-1.7575592) q[1];
x q[2];
rz(-2.7606332) q[3];
sx q[3];
rz(-1.6376312) q[3];
sx q[3];
rz(1.9240954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2603904) q[2];
sx q[2];
rz(-1.5217047) q[2];
sx q[2];
rz(0.41268665) q[2];
rz(-2.2212501) q[3];
sx q[3];
rz(-0.50655443) q[3];
sx q[3];
rz(1.9119561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5145787) q[0];
sx q[0];
rz(-0.98015061) q[0];
sx q[0];
rz(-0.29943109) q[0];
rz(-1.1224271) q[1];
sx q[1];
rz(-1.2010682) q[1];
sx q[1];
rz(1.0312414) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64410463) q[0];
sx q[0];
rz(-1.9369164) q[0];
sx q[0];
rz(3.1300504) q[0];
rz(2.756713) q[2];
sx q[2];
rz(-2.6483734) q[2];
sx q[2];
rz(0.2156336) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.448062) q[1];
sx q[1];
rz(-1.281257) q[1];
sx q[1];
rz(0.84622835) q[1];
rz(-pi) q[2];
rz(2.5599285) q[3];
sx q[3];
rz(-1.0137034) q[3];
sx q[3];
rz(2.795199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1624182) q[2];
sx q[2];
rz(-1.3434429) q[2];
sx q[2];
rz(-2.8988885) q[2];
rz(1.3001214) q[3];
sx q[3];
rz(-2.0827115) q[3];
sx q[3];
rz(-2.9259031) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8278787) q[0];
sx q[0];
rz(-0.21610459) q[0];
sx q[0];
rz(1.4208273) q[0];
rz(0.31006649) q[1];
sx q[1];
rz(-2.2646751) q[1];
sx q[1];
rz(-1.5975331) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2313854) q[0];
sx q[0];
rz(-0.69150309) q[0];
sx q[0];
rz(0.67759902) q[0];
rz(1.9126911) q[2];
sx q[2];
rz(-1.4871305) q[2];
sx q[2];
rz(0.25632206) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76031715) q[1];
sx q[1];
rz(-2.7997428) q[1];
sx q[1];
rz(-2.443497) q[1];
x q[2];
rz(-1.710911) q[3];
sx q[3];
rz(-1.7286674) q[3];
sx q[3];
rz(-0.50413361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1342423) q[2];
sx q[2];
rz(-1.4126567) q[2];
sx q[2];
rz(-1.3664112) q[2];
rz(2.2640696) q[3];
sx q[3];
rz(-1.4656504) q[3];
sx q[3];
rz(2.2519978) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9919745) q[0];
sx q[0];
rz(-1.5852954) q[0];
sx q[0];
rz(-1.4854767) q[0];
rz(0.86434518) q[1];
sx q[1];
rz(-1.0135916) q[1];
sx q[1];
rz(2.7085173) q[1];
rz(-2.6411459) q[2];
sx q[2];
rz(-1.7519578) q[2];
sx q[2];
rz(2.0448207) q[2];
rz(1.0917615) q[3];
sx q[3];
rz(-1.5216266) q[3];
sx q[3];
rz(-2.1306642) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
