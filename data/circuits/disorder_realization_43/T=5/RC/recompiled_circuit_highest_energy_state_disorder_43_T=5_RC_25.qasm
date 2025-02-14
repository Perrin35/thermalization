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
rz(-0.22696683) q[0];
sx q[0];
rz(-2.5300955) q[0];
sx q[0];
rz(-0.55846941) q[0];
rz(0.59977579) q[1];
sx q[1];
rz(1.37473) q[1];
sx q[1];
rz(9.3461499) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0315379) q[0];
sx q[0];
rz(-1.0194091) q[0];
sx q[0];
rz(0.64161513) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4457277) q[2];
sx q[2];
rz(-2.2457221) q[2];
sx q[2];
rz(1.9927911) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.78863555) q[1];
sx q[1];
rz(-2.87559) q[1];
sx q[1];
rz(2.8255844) q[1];
x q[2];
rz(-2.7897415) q[3];
sx q[3];
rz(-1.8738865) q[3];
sx q[3];
rz(1.6325655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.89995304) q[2];
sx q[2];
rz(-0.048154801) q[2];
sx q[2];
rz(2.107035) q[2];
rz(0.11040802) q[3];
sx q[3];
rz(-1.4805099) q[3];
sx q[3];
rz(-2.835527) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0050874) q[0];
sx q[0];
rz(-1.445048) q[0];
sx q[0];
rz(-1.9553631) q[0];
rz(-1.2677445) q[1];
sx q[1];
rz(-2.6823951) q[1];
sx q[1];
rz(-1.7523821) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2572309) q[0];
sx q[0];
rz(-0.62998929) q[0];
sx q[0];
rz(-1.0347741) q[0];
rz(-0.21135881) q[2];
sx q[2];
rz(-2.6508001) q[2];
sx q[2];
rz(0.15891128) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.6450883) q[1];
sx q[1];
rz(-1.8594655) q[1];
sx q[1];
rz(1.225807) q[1];
rz(-pi) q[2];
rz(1.4382382) q[3];
sx q[3];
rz(-2.5131559) q[3];
sx q[3];
rz(2.4579163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.17239751) q[2];
sx q[2];
rz(-1.1149422) q[2];
sx q[2];
rz(1.7051075) q[2];
rz(-1.8819594) q[3];
sx q[3];
rz(-2.6862222) q[3];
sx q[3];
rz(3.1148541) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36419511) q[0];
sx q[0];
rz(-1.7561678) q[0];
sx q[0];
rz(-1.9245603) q[0];
rz(-2.9265535) q[1];
sx q[1];
rz(-1.8937078) q[1];
sx q[1];
rz(1.4395813) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3453044) q[0];
sx q[0];
rz(-2.469696) q[0];
sx q[0];
rz(-2.94126) q[0];
rz(-pi) q[1];
rz(-1.7363988) q[2];
sx q[2];
rz(-2.6244727) q[2];
sx q[2];
rz(1.0696326) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.176217) q[1];
sx q[1];
rz(-1.4439519) q[1];
sx q[1];
rz(-1.3398257) q[1];
rz(2.2398364) q[3];
sx q[3];
rz(-0.34325156) q[3];
sx q[3];
rz(2.1829829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3415459) q[2];
sx q[2];
rz(-1.6365106) q[2];
sx q[2];
rz(2.3254507) q[2];
rz(-0.0084361313) q[3];
sx q[3];
rz(-0.71788994) q[3];
sx q[3];
rz(-1.5754023) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7208045) q[0];
sx q[0];
rz(-1.9705462) q[0];
sx q[0];
rz(0.26914832) q[0];
rz(1.3828269) q[1];
sx q[1];
rz(-2.5803284) q[1];
sx q[1];
rz(-0.47163481) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29971805) q[0];
sx q[0];
rz(-2.4025859) q[0];
sx q[0];
rz(1.7790545) q[0];
x q[1];
rz(-0.40377517) q[2];
sx q[2];
rz(-1.6036183) q[2];
sx q[2];
rz(0.36640247) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9593352) q[1];
sx q[1];
rz(-1.253009) q[1];
sx q[1];
rz(2.5885568) q[1];
x q[2];
rz(-2.7876623) q[3];
sx q[3];
rz(-2.2056863) q[3];
sx q[3];
rz(-0.77087444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.65583324) q[2];
sx q[2];
rz(-1.4071608) q[2];
sx q[2];
rz(1.2539697) q[2];
rz(-2.8433825) q[3];
sx q[3];
rz(-1.8658274) q[3];
sx q[3];
rz(-1.5378753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1042079) q[0];
sx q[0];
rz(-2.2315114) q[0];
sx q[0];
rz(-0.75697672) q[0];
rz(-1.0589927) q[1];
sx q[1];
rz(-1.4451566) q[1];
sx q[1];
rz(2.7991378) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9152074) q[0];
sx q[0];
rz(-1.7614375) q[0];
sx q[0];
rz(1.1676428) q[0];
rz(-pi) q[1];
rz(1.9497237) q[2];
sx q[2];
rz(-2.4034326) q[2];
sx q[2];
rz(-1.9110695) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7360644) q[1];
sx q[1];
rz(-1.2031415) q[1];
sx q[1];
rz(0.50761948) q[1];
rz(-pi) q[2];
rz(-2.5352372) q[3];
sx q[3];
rz(-2.6487051) q[3];
sx q[3];
rz(2.3541401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4464438) q[2];
sx q[2];
rz(-1.6764287) q[2];
sx q[2];
rz(-0.62015074) q[2];
rz(-2.5946963) q[3];
sx q[3];
rz(-2.9345025) q[3];
sx q[3];
rz(1.0774405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3243489) q[0];
sx q[0];
rz(-0.20340782) q[0];
sx q[0];
rz(1.9203583) q[0];
rz(2.0381894) q[1];
sx q[1];
rz(-1.9788479) q[1];
sx q[1];
rz(2.3308636) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9098783) q[0];
sx q[0];
rz(-2.8389769) q[0];
sx q[0];
rz(-0.55993863) q[0];
rz(-0.76580183) q[2];
sx q[2];
rz(-0.81379902) q[2];
sx q[2];
rz(-1.9794996) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74720464) q[1];
sx q[1];
rz(-1.6252246) q[1];
sx q[1];
rz(-1.7096976) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8437988) q[3];
sx q[3];
rz(-1.9341365) q[3];
sx q[3];
rz(-2.9693574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3409884) q[2];
sx q[2];
rz(-1.1977414) q[2];
sx q[2];
rz(1.9761337) q[2];
rz(0.10945877) q[3];
sx q[3];
rz(-1.0215003) q[3];
sx q[3];
rz(0.060733184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8282181) q[0];
sx q[0];
rz(-0.010638588) q[0];
sx q[0];
rz(-0.64205545) q[0];
rz(0.24318801) q[1];
sx q[1];
rz(-2.7311192) q[1];
sx q[1];
rz(0.64116716) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.084448) q[0];
sx q[0];
rz(-2.085745) q[0];
sx q[0];
rz(-0.65339974) q[0];
x q[1];
rz(-0.75058241) q[2];
sx q[2];
rz(-1.8403862) q[2];
sx q[2];
rz(-1.3113111) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7542759) q[1];
sx q[1];
rz(-2.2914304) q[1];
sx q[1];
rz(-0.44447036) q[1];
rz(1.9690597) q[3];
sx q[3];
rz(-1.7855111) q[3];
sx q[3];
rz(-2.1977148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0082561) q[2];
sx q[2];
rz(-1.9828321) q[2];
sx q[2];
rz(2.1978417) q[2];
rz(0.64822316) q[3];
sx q[3];
rz(-0.74130487) q[3];
sx q[3];
rz(-0.6664204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7159395) q[0];
sx q[0];
rz(-0.81099993) q[0];
sx q[0];
rz(-2.7124523) q[0];
rz(2.7375713) q[1];
sx q[1];
rz(-2.7262913) q[1];
sx q[1];
rz(-2.2228352) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63916535) q[0];
sx q[0];
rz(-2.5326061) q[0];
sx q[0];
rz(-0.43242411) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8317675) q[2];
sx q[2];
rz(-1.6918285) q[2];
sx q[2];
rz(2.3132035) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1240108) q[1];
sx q[1];
rz(-2.7582844) q[1];
sx q[1];
rz(0.95335754) q[1];
rz(-pi) q[2];
rz(-1.8848583) q[3];
sx q[3];
rz(-1.6410284) q[3];
sx q[3];
rz(-2.8633871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.018844133) q[2];
sx q[2];
rz(-0.65905535) q[2];
sx q[2];
rz(-1.3502236) q[2];
rz(-0.33251479) q[3];
sx q[3];
rz(-2.2612488) q[3];
sx q[3];
rz(1.3688709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3933082) q[0];
sx q[0];
rz(-0.64993334) q[0];
sx q[0];
rz(0.45676029) q[0];
rz(1.7204174) q[1];
sx q[1];
rz(-1.2799809) q[1];
sx q[1];
rz(-0.054904003) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47726705) q[0];
sx q[0];
rz(-2.3033963) q[0];
sx q[0];
rz(3.102836) q[0];
rz(-pi) q[1];
rz(2.3308067) q[2];
sx q[2];
rz(-0.56466757) q[2];
sx q[2];
rz(1.1155903) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0753638) q[1];
sx q[1];
rz(-1.7875331) q[1];
sx q[1];
rz(1.073887) q[1];
x q[2];
rz(-0.53590448) q[3];
sx q[3];
rz(-1.8632938) q[3];
sx q[3];
rz(-1.3939987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2533337) q[2];
sx q[2];
rz(-2.4312225) q[2];
sx q[2];
rz(-0.19676512) q[2];
rz(-2.0678068) q[3];
sx q[3];
rz(-2.00627) q[3];
sx q[3];
rz(0.73608583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1573023) q[0];
sx q[0];
rz(-0.8534011) q[0];
sx q[0];
rz(2.2421457) q[0];
rz(-2.8304214) q[1];
sx q[1];
rz(-1.2093465) q[1];
sx q[1];
rz(-1.7412294) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084564805) q[0];
sx q[0];
rz(-1.0175144) q[0];
sx q[0];
rz(0.24648376) q[0];
rz(-pi) q[1];
rz(-1.1323523) q[2];
sx q[2];
rz(-0.52784) q[2];
sx q[2];
rz(2.6486625) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.77620095) q[1];
sx q[1];
rz(-1.1836056) q[1];
sx q[1];
rz(-2.4939818) q[1];
x q[2];
rz(-2.4813249) q[3];
sx q[3];
rz(-2.3077624) q[3];
sx q[3];
rz(1.6136606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.060999) q[2];
sx q[2];
rz(-0.62005764) q[2];
sx q[2];
rz(1.9798123) q[2];
rz(-1.0787841) q[3];
sx q[3];
rz(-0.99083841) q[3];
sx q[3];
rz(0.88206464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5494736) q[0];
sx q[0];
rz(-1.5220806) q[0];
sx q[0];
rz(1.4694389) q[0];
rz(3.0081765) q[1];
sx q[1];
rz(-1.3795556) q[1];
sx q[1];
rz(2.1605927) q[1];
rz(0.0010112671) q[2];
sx q[2];
rz(-1.5173923) q[2];
sx q[2];
rz(-1.5569729) q[2];
rz(1.6422938) q[3];
sx q[3];
rz(-2.0713948) q[3];
sx q[3];
rz(-1.757302) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
