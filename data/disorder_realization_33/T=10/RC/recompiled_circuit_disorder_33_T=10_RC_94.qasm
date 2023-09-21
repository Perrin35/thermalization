OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6361976) q[0];
sx q[0];
rz(-0.27591053) q[0];
sx q[0];
rz(1.3077868) q[0];
rz(-2.0055327) q[1];
sx q[1];
rz(4.0772822) q[1];
sx q[1];
rz(4.7128591) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9309064) q[0];
sx q[0];
rz(-1.9084198) q[0];
sx q[0];
rz(-2.7772285) q[0];
rz(2.3762796) q[2];
sx q[2];
rz(-2.1047449) q[2];
sx q[2];
rz(3.090976) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5174487) q[1];
sx q[1];
rz(-1.9074719) q[1];
sx q[1];
rz(1.8271853) q[1];
rz(-pi) q[2];
rz(1.9340431) q[3];
sx q[3];
rz(-1.7540635) q[3];
sx q[3];
rz(2.7024384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2661665) q[2];
sx q[2];
rz(-0.29310075) q[2];
sx q[2];
rz(-2.0092633) q[2];
rz(1.6752361) q[3];
sx q[3];
rz(-1.8050067) q[3];
sx q[3];
rz(2.1291389) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9448626) q[0];
sx q[0];
rz(-2.9319627) q[0];
sx q[0];
rz(0.18584132) q[0];
rz(-0.56022412) q[1];
sx q[1];
rz(-1.8461684) q[1];
sx q[1];
rz(2.9247608) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8502055) q[0];
sx q[0];
rz(-2.4225525) q[0];
sx q[0];
rz(-1.1262116) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0160604) q[2];
sx q[2];
rz(-1.1604571) q[2];
sx q[2];
rz(-2.1330657) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.195897) q[1];
sx q[1];
rz(-0.90598124) q[1];
sx q[1];
rz(2.871454) q[1];
rz(-0.64036815) q[3];
sx q[3];
rz(-2.159517) q[3];
sx q[3];
rz(1.1191739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8314787) q[2];
sx q[2];
rz(-0.82565132) q[2];
sx q[2];
rz(1.2878093) q[2];
rz(0.76256049) q[3];
sx q[3];
rz(-1.1688787) q[3];
sx q[3];
rz(0.30502239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.6771616) q[0];
sx q[0];
rz(-2.7966249) q[0];
sx q[0];
rz(-2.537354) q[0];
rz(-1.3263946) q[1];
sx q[1];
rz(-1.3605958) q[1];
sx q[1];
rz(-0.93260971) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55721012) q[0];
sx q[0];
rz(-0.32214468) q[0];
sx q[0];
rz(1.4003217) q[0];
rz(-pi) q[1];
rz(-1.3435752) q[2];
sx q[2];
rz(-1.3933239) q[2];
sx q[2];
rz(0.15572671) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.697726) q[1];
sx q[1];
rz(-1.4523456) q[1];
sx q[1];
rz(-2.5812134) q[1];
rz(-1.014939) q[3];
sx q[3];
rz(-1.9564637) q[3];
sx q[3];
rz(-2.8876497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9937667) q[2];
sx q[2];
rz(-1.0819165) q[2];
sx q[2];
rz(2.0489342) q[2];
rz(-2.5993733) q[3];
sx q[3];
rz(-1.0850302) q[3];
sx q[3];
rz(-2.1742163) q[3];
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
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7820691) q[0];
sx q[0];
rz(-0.096465915) q[0];
sx q[0];
rz(2.6413667) q[0];
rz(-2.3362828) q[1];
sx q[1];
rz(-1.9814682) q[1];
sx q[1];
rz(-1.6436228) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.362975) q[0];
sx q[0];
rz(-0.59016363) q[0];
sx q[0];
rz(-2.1182563) q[0];
rz(-pi) q[1];
rz(-0.056604071) q[2];
sx q[2];
rz(-1.3805693) q[2];
sx q[2];
rz(0.20400001) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.56470358) q[1];
sx q[1];
rz(-2.7895045) q[1];
sx q[1];
rz(2.4114154) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1180531) q[3];
sx q[3];
rz(-1.9948043) q[3];
sx q[3];
rz(0.31183576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3952289) q[2];
sx q[2];
rz(-2.5791898) q[2];
sx q[2];
rz(2.4397819) q[2];
rz(-0.83135215) q[3];
sx q[3];
rz(-2.1777007) q[3];
sx q[3];
rz(-2.519616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2410626) q[0];
sx q[0];
rz(-0.59589544) q[0];
sx q[0];
rz(-0.81533122) q[0];
rz(1.5218081) q[1];
sx q[1];
rz(-0.83414572) q[1];
sx q[1];
rz(1.048208) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4002776) q[0];
sx q[0];
rz(-1.8198697) q[0];
sx q[0];
rz(1.2988017) q[0];
rz(-pi) q[1];
rz(-1.5830718) q[2];
sx q[2];
rz(-0.91931146) q[2];
sx q[2];
rz(0.94142454) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9427467) q[1];
sx q[1];
rz(-1.2578576) q[1];
sx q[1];
rz(-1.8206157) q[1];
rz(2.1861595) q[3];
sx q[3];
rz(-1.9121998) q[3];
sx q[3];
rz(1.1036901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52577019) q[2];
sx q[2];
rz(-2.5746391) q[2];
sx q[2];
rz(1.099951) q[2];
rz(-0.82529092) q[3];
sx q[3];
rz(-2.0992978) q[3];
sx q[3];
rz(0.88551372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0681756) q[0];
sx q[0];
rz(-2.5475579) q[0];
sx q[0];
rz(-2.2391879) q[0];
rz(1.0166608) q[1];
sx q[1];
rz(-2.0817751) q[1];
sx q[1];
rz(-3.0117603) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31045612) q[0];
sx q[0];
rz(-1.2015011) q[0];
sx q[0];
rz(2.5184758) q[0];
rz(-pi) q[1];
rz(1.0340704) q[2];
sx q[2];
rz(-0.64646361) q[2];
sx q[2];
rz(1.549364) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2643913) q[1];
sx q[1];
rz(-2.1346722) q[1];
sx q[1];
rz(-1.1370204) q[1];
rz(-pi) q[2];
rz(-0.54883212) q[3];
sx q[3];
rz(-1.4026814) q[3];
sx q[3];
rz(2.7041534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8292024) q[2];
sx q[2];
rz(-0.94909334) q[2];
sx q[2];
rz(0.20425805) q[2];
rz(-1.9355109) q[3];
sx q[3];
rz(-1.6198502) q[3];
sx q[3];
rz(0.23541418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4181353) q[0];
sx q[0];
rz(-1.8122939) q[0];
sx q[0];
rz(-1.6947421) q[0];
rz(1.2591259) q[1];
sx q[1];
rz(-2.1513758) q[1];
sx q[1];
rz(-0.68626219) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9335564) q[0];
sx q[0];
rz(-2.4718923) q[0];
sx q[0];
rz(0.74525381) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83001901) q[2];
sx q[2];
rz(-2.0247211) q[2];
sx q[2];
rz(-1.9759076) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8441019) q[1];
sx q[1];
rz(-0.19980783) q[1];
sx q[1];
rz(-1.954345) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8920184) q[3];
sx q[3];
rz(-1.1676844) q[3];
sx q[3];
rz(-1.5461127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4454322) q[2];
sx q[2];
rz(-1.7636718) q[2];
sx q[2];
rz(0.0017722842) q[2];
rz(2.5799675) q[3];
sx q[3];
rz(-2.2300945) q[3];
sx q[3];
rz(-1.5047489) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5381662) q[0];
sx q[0];
rz(-0.68646938) q[0];
sx q[0];
rz(-1.6954533) q[0];
rz(-0.7810477) q[1];
sx q[1];
rz(-1.8361517) q[1];
sx q[1];
rz(-1.6400281) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7088889) q[0];
sx q[0];
rz(-2.6484657) q[0];
sx q[0];
rz(1.6119484) q[0];
rz(-1.5090452) q[2];
sx q[2];
rz(-1.5973063) q[2];
sx q[2];
rz(2.3960631) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0280684) q[1];
sx q[1];
rz(-1.6599732) q[1];
sx q[1];
rz(-2.813617) q[1];
rz(-0.73236671) q[3];
sx q[3];
rz(-2.3673956) q[3];
sx q[3];
rz(-2.458651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7897196) q[2];
sx q[2];
rz(-1.3616273) q[2];
sx q[2];
rz(-1.8224576) q[2];
rz(-1.9296648) q[3];
sx q[3];
rz(-1.8550248) q[3];
sx q[3];
rz(-0.31931988) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8050352) q[0];
sx q[0];
rz(-0.55258495) q[0];
sx q[0];
rz(1.2040899) q[0];
rz(-2.7583292) q[1];
sx q[1];
rz(-2.6158694) q[1];
sx q[1];
rz(-0.35167545) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29031819) q[0];
sx q[0];
rz(-1.9418678) q[0];
sx q[0];
rz(-1.0822269) q[0];
x q[1];
rz(0.69182379) q[2];
sx q[2];
rz(-1.8080538) q[2];
sx q[2];
rz(-1.1292063) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0723567) q[1];
sx q[1];
rz(-1.623739) q[1];
sx q[1];
rz(-2.3318021) q[1];
rz(-3.0319801) q[3];
sx q[3];
rz(-1.5188367) q[3];
sx q[3];
rz(2.6250641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7982771) q[2];
sx q[2];
rz(-2.0337992) q[2];
sx q[2];
rz(1.8593672) q[2];
rz(1.6451689) q[3];
sx q[3];
rz(-1.6069501) q[3];
sx q[3];
rz(-1.055868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6431817) q[0];
sx q[0];
rz(-1.8739941) q[0];
sx q[0];
rz(2.9472651) q[0];
rz(2.1037897) q[1];
sx q[1];
rz(-2.5732645) q[1];
sx q[1];
rz(-1.0338354) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68034222) q[0];
sx q[0];
rz(-0.67078062) q[0];
sx q[0];
rz(-1.781342) q[0];
rz(-pi) q[1];
rz(-2.2718272) q[2];
sx q[2];
rz(-1.9734205) q[2];
sx q[2];
rz(-1.0588156) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.49627134) q[1];
sx q[1];
rz(-1.4693345) q[1];
sx q[1];
rz(-0.99781499) q[1];
x q[2];
rz(0.25063534) q[3];
sx q[3];
rz(-0.72892979) q[3];
sx q[3];
rz(-0.38754101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0795435) q[2];
sx q[2];
rz(-0.94576183) q[2];
sx q[2];
rz(-2.5058084) q[2];
rz(-0.27030269) q[3];
sx q[3];
rz(-0.79939866) q[3];
sx q[3];
rz(1.5283782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4476267) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(1.7059965) q[1];
sx q[1];
rz(-1.5626848) q[1];
sx q[1];
rz(-2.3609153) q[1];
rz(0.96799093) q[2];
sx q[2];
rz(-1.5445166) q[2];
sx q[2];
rz(-1.9894285) q[2];
rz(-1.900832) q[3];
sx q[3];
rz(-1.6374554) q[3];
sx q[3];
rz(-1.0367254) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];