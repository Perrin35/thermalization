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
rz(-4.0487657) q[0];
sx q[0];
rz(2.2450759) q[0];
sx q[0];
rz(9.4656691) q[0];
rz(1.8812802) q[1];
sx q[1];
rz(-0.44668302) q[1];
sx q[1];
rz(3.0466411) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68317682) q[0];
sx q[0];
rz(-1.9703034) q[0];
sx q[0];
rz(1.3241121) q[0];
rz(-1.3407732) q[2];
sx q[2];
rz(-2.1162976) q[2];
sx q[2];
rz(0.70822424) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5168001) q[1];
sx q[1];
rz(-1.9408556) q[1];
sx q[1];
rz(0.79536749) q[1];
x q[2];
rz(0.818472) q[3];
sx q[3];
rz(-2.5305037) q[3];
sx q[3];
rz(-2.1095534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.2824668) q[2];
sx q[2];
rz(-0.60255113) q[2];
sx q[2];
rz(1.1653384) q[2];
rz(2.2230542) q[3];
sx q[3];
rz(-3.0375752) q[3];
sx q[3];
rz(1.2134086) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4024046) q[0];
sx q[0];
rz(-0.63888752) q[0];
sx q[0];
rz(2.8366587) q[0];
rz(2.1091499) q[1];
sx q[1];
rz(-2.86125) q[1];
sx q[1];
rz(-0.84411821) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4798571) q[0];
sx q[0];
rz(-0.2847372) q[0];
sx q[0];
rz(0.99951331) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0424625) q[2];
sx q[2];
rz(-0.24046637) q[2];
sx q[2];
rz(3.0853396) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6355953) q[1];
sx q[1];
rz(-1.330101) q[1];
sx q[1];
rz(-2.7240915) q[1];
rz(-pi) q[2];
rz(1.5553383) q[3];
sx q[3];
rz(-0.94383701) q[3];
sx q[3];
rz(2.5687237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.9742763) q[2];
sx q[2];
rz(-2.4582489) q[2];
sx q[2];
rz(2.5565476) q[2];
rz(-2.282418) q[3];
sx q[3];
rz(-1.1480568) q[3];
sx q[3];
rz(-1.1332716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0228731) q[0];
sx q[0];
rz(-2.3206503) q[0];
sx q[0];
rz(-3.0430479) q[0];
rz(0.56602829) q[1];
sx q[1];
rz(-0.84605828) q[1];
sx q[1];
rz(-0.50618323) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0289259) q[0];
sx q[0];
rz(-0.62917626) q[0];
sx q[0];
rz(2.5553246) q[0];
rz(-pi) q[1];
rz(-2.3213941) q[2];
sx q[2];
rz(-1.9026105) q[2];
sx q[2];
rz(0.92727509) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4892046) q[1];
sx q[1];
rz(-1.4474157) q[1];
sx q[1];
rz(1.3493933) q[1];
rz(-0.99870317) q[3];
sx q[3];
rz(-0.65352189) q[3];
sx q[3];
rz(1.3746978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.11518662) q[2];
sx q[2];
rz(-1.2135442) q[2];
sx q[2];
rz(-3.0565267) q[2];
rz(2.3885942) q[3];
sx q[3];
rz(-2.4716061) q[3];
sx q[3];
rz(1.5299214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.531115) q[0];
sx q[0];
rz(-1.6401289) q[0];
sx q[0];
rz(-2.6016972) q[0];
rz(-0.029622948) q[1];
sx q[1];
rz(-2.3952775) q[1];
sx q[1];
rz(1.7866887) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8919075) q[0];
sx q[0];
rz(-2.2985795) q[0];
sx q[0];
rz(-0.2951727) q[0];
rz(-pi) q[1];
rz(-2.7672655) q[2];
sx q[2];
rz(-2.3951963) q[2];
sx q[2];
rz(1.0837006) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.190092) q[1];
sx q[1];
rz(-1.8279549) q[1];
sx q[1];
rz(2.5567104) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6605234) q[3];
sx q[3];
rz(-2.3589241) q[3];
sx q[3];
rz(1.6663446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4386091) q[2];
sx q[2];
rz(-0.97038022) q[2];
sx q[2];
rz(-2.2461829) q[2];
rz(-1.6926951) q[3];
sx q[3];
rz(-1.1619032) q[3];
sx q[3];
rz(2.7321775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92622906) q[0];
sx q[0];
rz(-2.138593) q[0];
sx q[0];
rz(-3.1410425) q[0];
rz(-2.6161361) q[1];
sx q[1];
rz(-0.84392396) q[1];
sx q[1];
rz(-2.1702683) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7346974) q[0];
sx q[0];
rz(-0.79888035) q[0];
sx q[0];
rz(-0.44476923) q[0];
rz(2.7883456) q[2];
sx q[2];
rz(-1.9611036) q[2];
sx q[2];
rz(2.3164904) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3974053) q[1];
sx q[1];
rz(-2.1402855) q[1];
sx q[1];
rz(-2.2389212) q[1];
rz(0.21495238) q[3];
sx q[3];
rz(-1.8773517) q[3];
sx q[3];
rz(-2.7988899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9690669) q[2];
sx q[2];
rz(-1.0774287) q[2];
sx q[2];
rz(1.7871008) q[2];
rz(-0.6012249) q[3];
sx q[3];
rz(-1.0433334) q[3];
sx q[3];
rz(-1.5662947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.4300267) q[0];
sx q[0];
rz(-2.8373748) q[0];
sx q[0];
rz(-0.29997224) q[0];
rz(-0.9777588) q[1];
sx q[1];
rz(-0.58039665) q[1];
sx q[1];
rz(2.207644) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2364607) q[0];
sx q[0];
rz(-1.2398402) q[0];
sx q[0];
rz(-1.4230498) q[0];
x q[1];
rz(0.63542346) q[2];
sx q[2];
rz(-0.72180702) q[2];
sx q[2];
rz(1.3363163) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.63864795) q[1];
sx q[1];
rz(-0.80728045) q[1];
sx q[1];
rz(2.817383) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4277064) q[3];
sx q[3];
rz(-0.51091226) q[3];
sx q[3];
rz(0.78277222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0705491) q[2];
sx q[2];
rz(-1.0087548) q[2];
sx q[2];
rz(1.1318995) q[2];
rz(-1.2567358) q[3];
sx q[3];
rz(-2.3102424) q[3];
sx q[3];
rz(1.4164213) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8202332) q[0];
sx q[0];
rz(-1.8978523) q[0];
sx q[0];
rz(-2.5309122) q[0];
rz(-0.75675476) q[1];
sx q[1];
rz(-2.7532531) q[1];
sx q[1];
rz(1.6441708) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1377378) q[0];
sx q[0];
rz(-0.12659368) q[0];
sx q[0];
rz(-2.9403482) q[0];
x q[1];
rz(0.87074222) q[2];
sx q[2];
rz(-0.53360924) q[2];
sx q[2];
rz(-2.6117532) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.25631902) q[1];
sx q[1];
rz(-0.65769629) q[1];
sx q[1];
rz(0.17430507) q[1];
rz(-0.11388643) q[3];
sx q[3];
rz(-0.72491693) q[3];
sx q[3];
rz(0.33958401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.18621592) q[2];
sx q[2];
rz(-0.66408855) q[2];
sx q[2];
rz(-0.60626283) q[2];
rz(2.9210505) q[3];
sx q[3];
rz(-1.8044148) q[3];
sx q[3];
rz(2.7748599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99899387) q[0];
sx q[0];
rz(-1.4936916) q[0];
sx q[0];
rz(2.2338474) q[0];
rz(-1.6084464) q[1];
sx q[1];
rz(-2.1870859) q[1];
sx q[1];
rz(-3.122701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7291537) q[0];
sx q[0];
rz(-1.6866418) q[0];
sx q[0];
rz(-2.6687268) q[0];
rz(-pi) q[1];
rz(2.9540145) q[2];
sx q[2];
rz(-0.30454985) q[2];
sx q[2];
rz(1.9001324) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6396913) q[1];
sx q[1];
rz(-2.0328201) q[1];
sx q[1];
rz(1.6056745) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3109679) q[3];
sx q[3];
rz(-1.1234094) q[3];
sx q[3];
rz(0.4448286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3733526) q[2];
sx q[2];
rz(-1.5529996) q[2];
sx q[2];
rz(-0.24660435) q[2];
rz(1.8433579) q[3];
sx q[3];
rz(-1.4539098) q[3];
sx q[3];
rz(-1.8728144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3835555) q[0];
sx q[0];
rz(-1.0724496) q[0];
sx q[0];
rz(-2.2776336) q[0];
rz(1.9268688) q[1];
sx q[1];
rz(-2.0236969) q[1];
sx q[1];
rz(2.5606959) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37947919) q[0];
sx q[0];
rz(-2.223928) q[0];
sx q[0];
rz(0.94193052) q[0];
x q[1];
rz(3.1031392) q[2];
sx q[2];
rz(-0.83476257) q[2];
sx q[2];
rz(2.9308386) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.89199663) q[1];
sx q[1];
rz(-2.432715) q[1];
sx q[1];
rz(-0.53541553) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8472927) q[3];
sx q[3];
rz(-2.3307335) q[3];
sx q[3];
rz(3.1408059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2555799) q[2];
sx q[2];
rz(-2.6801127) q[2];
sx q[2];
rz(-0.18730051) q[2];
rz(-2.1646132) q[3];
sx q[3];
rz(-1.6410442) q[3];
sx q[3];
rz(1.2788844) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5115857) q[0];
sx q[0];
rz(-0.66660175) q[0];
sx q[0];
rz(-1.8774207) q[0];
rz(0.97201792) q[1];
sx q[1];
rz(-2.3897901) q[1];
sx q[1];
rz(-1.972563) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3384333) q[0];
sx q[0];
rz(-1.2375323) q[0];
sx q[0];
rz(-2.9985371) q[0];
rz(-pi) q[1];
rz(0.75836597) q[2];
sx q[2];
rz(-2.6347199) q[2];
sx q[2];
rz(0.015794347) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.15195512) q[1];
sx q[1];
rz(-2.440976) q[1];
sx q[1];
rz(1.408427) q[1];
rz(-2.8253978) q[3];
sx q[3];
rz(-1.6796111) q[3];
sx q[3];
rz(-1.8121882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.4740037) q[2];
sx q[2];
rz(-1.7550125) q[2];
sx q[2];
rz(-2.7645195) q[2];
rz(0.22370473) q[3];
sx q[3];
rz(-0.5880028) q[3];
sx q[3];
rz(-1.2288176) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2558462) q[0];
sx q[0];
rz(-1.3930014) q[0];
sx q[0];
rz(-0.2572671) q[0];
rz(-0.59355758) q[1];
sx q[1];
rz(-0.7919574) q[1];
sx q[1];
rz(1.6603464) q[1];
rz(-2.42057) q[2];
sx q[2];
rz(-2.1446054) q[2];
sx q[2];
rz(-2.9531324) q[2];
rz(-1.1720584) q[3];
sx q[3];
rz(-1.8137365) q[3];
sx q[3];
rz(0.89567281) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
