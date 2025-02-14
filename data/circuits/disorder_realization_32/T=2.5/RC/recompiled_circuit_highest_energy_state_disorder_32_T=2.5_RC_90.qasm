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
rz(0.12953144) q[0];
sx q[0];
rz(3.01053) q[0];
sx q[0];
rz(13.170903) q[0];
rz(-0.52855748) q[1];
sx q[1];
rz(2.8837535) q[1];
sx q[1];
rz(16.937994) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2211232) q[0];
sx q[0];
rz(-1.9029494) q[0];
sx q[0];
rz(1.8355379) q[0];
x q[1];
rz(2.5276353) q[2];
sx q[2];
rz(-1.3836622) q[2];
sx q[2];
rz(0.89636114) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.9046272) q[1];
sx q[1];
rz(-0.89348999) q[1];
sx q[1];
rz(-0.64579247) q[1];
rz(-pi) q[2];
rz(2.0666201) q[3];
sx q[3];
rz(-1.6832441) q[3];
sx q[3];
rz(-1.5898286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0872385) q[2];
sx q[2];
rz(-1.9638991) q[2];
sx q[2];
rz(-0.13620201) q[2];
rz(-0.44998351) q[3];
sx q[3];
rz(-2.4583702) q[3];
sx q[3];
rz(2.3581678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0403274) q[0];
sx q[0];
rz(-2.3430921) q[0];
sx q[0];
rz(-2.878317) q[0];
rz(1.0385665) q[1];
sx q[1];
rz(-2.7949605) q[1];
sx q[1];
rz(2.3668049) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8717125) q[0];
sx q[0];
rz(-2.1795666) q[0];
sx q[0];
rz(1.3193498) q[0];
rz(2.3900142) q[2];
sx q[2];
rz(-0.84779352) q[2];
sx q[2];
rz(0.64848775) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0475743) q[1];
sx q[1];
rz(-1.9500004) q[1];
sx q[1];
rz(-0.039888558) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4382884) q[3];
sx q[3];
rz(-1.9580132) q[3];
sx q[3];
rz(-2.4037698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7210377) q[2];
sx q[2];
rz(-1.5169531) q[2];
sx q[2];
rz(-0.59845412) q[2];
rz(-1.362644) q[3];
sx q[3];
rz(-0.88038954) q[3];
sx q[3];
rz(-0.27073282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8552928) q[0];
sx q[0];
rz(-2.6787651) q[0];
sx q[0];
rz(-1.6337974) q[0];
rz(0.15448054) q[1];
sx q[1];
rz(-1.721761) q[1];
sx q[1];
rz(2.3522164) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95473991) q[0];
sx q[0];
rz(-0.53760872) q[0];
sx q[0];
rz(-2.9718999) q[0];
x q[1];
rz(-2.6811483) q[2];
sx q[2];
rz(-2.8297462) q[2];
sx q[2];
rz(-2.0058035) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4142609) q[1];
sx q[1];
rz(-1.8751029) q[1];
sx q[1];
rz(2.2060966) q[1];
rz(-pi) q[2];
rz(-0.49279883) q[3];
sx q[3];
rz(-0.86233739) q[3];
sx q[3];
rz(0.81564834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.35884759) q[2];
sx q[2];
rz(-2.2300356) q[2];
sx q[2];
rz(1.6645128) q[2];
rz(-1.9278256) q[3];
sx q[3];
rz(-1.8002847) q[3];
sx q[3];
rz(-1.7267797) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3690935) q[0];
sx q[0];
rz(-1.0558244) q[0];
sx q[0];
rz(2.5308727) q[0];
rz(0.67131132) q[1];
sx q[1];
rz(-1.6339615) q[1];
sx q[1];
rz(-1.7866561) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5089534) q[0];
sx q[0];
rz(-2.0193978) q[0];
sx q[0];
rz(-2.6030356) q[0];
rz(1.1995537) q[2];
sx q[2];
rz(-1.4535731) q[2];
sx q[2];
rz(-0.10088149) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0860885) q[1];
sx q[1];
rz(-1.633606) q[1];
sx q[1];
rz(1.8163866) q[1];
rz(-pi) q[2];
rz(0.98396222) q[3];
sx q[3];
rz(-1.4915183) q[3];
sx q[3];
rz(-2.6266499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9352202) q[2];
sx q[2];
rz(-2.361203) q[2];
sx q[2];
rz(-2.8727403) q[2];
rz(-2.3937461) q[3];
sx q[3];
rz(-2.3220389) q[3];
sx q[3];
rz(-1.4486754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.185323) q[0];
sx q[0];
rz(-1.7522426) q[0];
sx q[0];
rz(-2.4638033) q[0];
rz(0.14450821) q[1];
sx q[1];
rz(-2.3542207) q[1];
sx q[1];
rz(1.4871303) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6078454) q[0];
sx q[0];
rz(-0.13838875) q[0];
sx q[0];
rz(-1.3800623) q[0];
rz(-pi) q[1];
rz(3.0838548) q[2];
sx q[2];
rz(-1.1694093) q[2];
sx q[2];
rz(-0.58338469) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.9095535) q[1];
sx q[1];
rz(-2.0318446) q[1];
sx q[1];
rz(-0.98566815) q[1];
rz(-pi) q[2];
rz(-1.4492338) q[3];
sx q[3];
rz(-2.7633177) q[3];
sx q[3];
rz(0.88751436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3378478) q[2];
sx q[2];
rz(-2.8779112) q[2];
sx q[2];
rz(-2.7860876) q[2];
rz(2.096094) q[3];
sx q[3];
rz(-1.239536) q[3];
sx q[3];
rz(-0.56226468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1657555) q[0];
sx q[0];
rz(-0.58335692) q[0];
sx q[0];
rz(3.052886) q[0];
rz(1.6070131) q[1];
sx q[1];
rz(-1.3101703) q[1];
sx q[1];
rz(-3.065899) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3189118) q[0];
sx q[0];
rz(-1.7633798) q[0];
sx q[0];
rz(-0.67355021) q[0];
rz(-pi) q[1];
rz(2.3024302) q[2];
sx q[2];
rz(-1.5328836) q[2];
sx q[2];
rz(-2.1597852) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5548426) q[1];
sx q[1];
rz(-1.4289218) q[1];
sx q[1];
rz(2.1290522) q[1];
rz(-pi) q[2];
x q[2];
rz(2.663732) q[3];
sx q[3];
rz(-0.63242542) q[3];
sx q[3];
rz(-2.6101108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.82464108) q[2];
sx q[2];
rz(-1.3836766) q[2];
sx q[2];
rz(-1.0392044) q[2];
rz(-2.9292987) q[3];
sx q[3];
rz(-2.0391235) q[3];
sx q[3];
rz(-1.645741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067326389) q[0];
sx q[0];
rz(-0.75263158) q[0];
sx q[0];
rz(-2.0972032) q[0];
rz(0.77230612) q[1];
sx q[1];
rz(-0.65985313) q[1];
sx q[1];
rz(0.73371249) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.808107) q[0];
sx q[0];
rz(-0.86577387) q[0];
sx q[0];
rz(0.66410983) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15007054) q[2];
sx q[2];
rz(-1.9378621) q[2];
sx q[2];
rz(2.31524) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1024688) q[1];
sx q[1];
rz(-0.3141292) q[1];
sx q[1];
rz(0.38508319) q[1];
x q[2];
rz(2.3055607) q[3];
sx q[3];
rz(-1.4905018) q[3];
sx q[3];
rz(-1.3463595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.61651984) q[2];
sx q[2];
rz(-2.6313582) q[2];
sx q[2];
rz(0.60866848) q[2];
rz(1.0673374) q[3];
sx q[3];
rz(-2.267024) q[3];
sx q[3];
rz(2.5909891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6390425) q[0];
sx q[0];
rz(-0.35824963) q[0];
sx q[0];
rz(1.6453561) q[0];
rz(1.7773588) q[1];
sx q[1];
rz(-0.61449209) q[1];
sx q[1];
rz(-0.38633698) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57418121) q[0];
sx q[0];
rz(-1.8339388) q[0];
sx q[0];
rz(-2.4469923) q[0];
x q[1];
rz(0.069938439) q[2];
sx q[2];
rz(-1.3856263) q[2];
sx q[2];
rz(2.1500146) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.92703351) q[1];
sx q[1];
rz(-0.3541358) q[1];
sx q[1];
rz(2.781809) q[1];
rz(-pi) q[2];
x q[2];
rz(1.531485) q[3];
sx q[3];
rz(-2.1609801) q[3];
sx q[3];
rz(-3.085768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6284457) q[2];
sx q[2];
rz(-1.3725504) q[2];
sx q[2];
rz(-0.15303843) q[2];
rz(-0.89093527) q[3];
sx q[3];
rz(-2.1864083) q[3];
sx q[3];
rz(-0.74131596) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9855758) q[0];
sx q[0];
rz(-2.9122536) q[0];
sx q[0];
rz(-2.7942221) q[0];
rz(-3.103745) q[1];
sx q[1];
rz(-1.9719351) q[1];
sx q[1];
rz(-2.6191424) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7670531) q[0];
sx q[0];
rz(-1.660083) q[0];
sx q[0];
rz(-1.1566628) q[0];
rz(3.1407479) q[2];
sx q[2];
rz(-0.36549308) q[2];
sx q[2];
rz(-1.5216684) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47679893) q[1];
sx q[1];
rz(-1.069624) q[1];
sx q[1];
rz(2.0528) q[1];
x q[2];
rz(-0.96790989) q[3];
sx q[3];
rz(-2.4914352) q[3];
sx q[3];
rz(-3.0126743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.6309506) q[2];
sx q[2];
rz(-2.6783671) q[2];
sx q[2];
rz(1.9412712) q[2];
rz(0.55780324) q[3];
sx q[3];
rz(-1.5188981) q[3];
sx q[3];
rz(0.38448486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8484304) q[0];
sx q[0];
rz(-0.98386216) q[0];
sx q[0];
rz(-1.4137319) q[0];
rz(-0.14104715) q[1];
sx q[1];
rz(-1.7660564) q[1];
sx q[1];
rz(-2.486855) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.065718) q[0];
sx q[0];
rz(-1.362934) q[0];
sx q[0];
rz(1.8130568) q[0];
rz(0.59487307) q[2];
sx q[2];
rz(-2.0030177) q[2];
sx q[2];
rz(3.1405666) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6605986) q[1];
sx q[1];
rz(-2.8428322) q[1];
sx q[1];
rz(-2.1887652) q[1];
rz(-pi) q[2];
x q[2];
rz(1.039417) q[3];
sx q[3];
rz(-2.1221042) q[3];
sx q[3];
rz(-1.8816063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.35141382) q[2];
sx q[2];
rz(-1.979579) q[2];
sx q[2];
rz(-0.72511017) q[2];
rz(-2.2509947) q[3];
sx q[3];
rz(-2.2913439) q[3];
sx q[3];
rz(-2.5543673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9643322) q[0];
sx q[0];
rz(-1.5252508) q[0];
sx q[0];
rz(-1.9510212) q[0];
rz(1.1217077) q[1];
sx q[1];
rz(-1.5737166) q[1];
sx q[1];
rz(2.166688) q[1];
rz(2.4779392) q[2];
sx q[2];
rz(-1.2812231) q[2];
sx q[2];
rz(-0.44322586) q[2];
rz(2.4287281) q[3];
sx q[3];
rz(-1.6920857) q[3];
sx q[3];
rz(-1.6209775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
