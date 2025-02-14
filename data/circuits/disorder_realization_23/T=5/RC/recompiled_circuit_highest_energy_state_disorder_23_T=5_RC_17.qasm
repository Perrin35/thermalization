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
rz(-2.5858606) q[0];
sx q[0];
rz(-1.2795804) q[0];
sx q[0];
rz(-2.8137141) q[0];
rz(-2.9887587) q[1];
sx q[1];
rz(-2.6522377) q[1];
sx q[1];
rz(2.1305003) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7231434) q[0];
sx q[0];
rz(-1.7307502) q[0];
sx q[0];
rz(2.5470995) q[0];
x q[1];
rz(1.3921803) q[2];
sx q[2];
rz(-1.9733323) q[2];
sx q[2];
rz(-2.2312763) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9823237) q[1];
sx q[1];
rz(-1.9344653) q[1];
sx q[1];
rz(2.3669142) q[1];
x q[2];
rz(-1.1527083) q[3];
sx q[3];
rz(-1.3527414) q[3];
sx q[3];
rz(-2.5421028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.27935394) q[2];
sx q[2];
rz(-0.78247672) q[2];
sx q[2];
rz(-1.5834825) q[2];
rz(0.33485788) q[3];
sx q[3];
rz(-1.0840651) q[3];
sx q[3];
rz(0.28086942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25664499) q[0];
sx q[0];
rz(-2.5769951) q[0];
sx q[0];
rz(0.33194342) q[0];
rz(-2.7815107) q[1];
sx q[1];
rz(-1.304345) q[1];
sx q[1];
rz(-2.8538381) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028988722) q[0];
sx q[0];
rz(-1.3224012) q[0];
sx q[0];
rz(-1.7775787) q[0];
x q[1];
rz(-1.8910725) q[2];
sx q[2];
rz(-0.8647635) q[2];
sx q[2];
rz(1.4780188) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0313686) q[1];
sx q[1];
rz(-1.2249631) q[1];
sx q[1];
rz(2.7909173) q[1];
rz(-pi) q[2];
rz(-2.5310181) q[3];
sx q[3];
rz(-2.861851) q[3];
sx q[3];
rz(0.66204643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.45431367) q[2];
sx q[2];
rz(-1.5328898) q[2];
sx q[2];
rz(1.3909371) q[2];
rz(0.85919291) q[3];
sx q[3];
rz(-1.8682559) q[3];
sx q[3];
rz(0.19101645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.745382) q[0];
sx q[0];
rz(-1.5353545) q[0];
sx q[0];
rz(-0.23042738) q[0];
rz(2.6241809) q[1];
sx q[1];
rz(-2.0473862) q[1];
sx q[1];
rz(0.28883019) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17265564) q[0];
sx q[0];
rz(-0.41303493) q[0];
sx q[0];
rz(-1.719219) q[0];
x q[1];
rz(-2.9952461) q[2];
sx q[2];
rz(-0.47292559) q[2];
sx q[2];
rz(0.84640102) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3687262) q[1];
sx q[1];
rz(-1.9802367) q[1];
sx q[1];
rz(-3.0837497) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9656319) q[3];
sx q[3];
rz(-1.7909414) q[3];
sx q[3];
rz(-2.1484321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1032224) q[2];
sx q[2];
rz(-0.73870814) q[2];
sx q[2];
rz(-0.040741097) q[2];
rz(0.1013969) q[3];
sx q[3];
rz(-1.8056168) q[3];
sx q[3];
rz(2.0749157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3876225) q[0];
sx q[0];
rz(-0.0059703537) q[0];
sx q[0];
rz(0.23400865) q[0];
rz(-0.19800828) q[1];
sx q[1];
rz(-1.0375236) q[1];
sx q[1];
rz(-2.1580946) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4490383) q[0];
sx q[0];
rz(-0.48547599) q[0];
sx q[0];
rz(-1.6621157) q[0];
x q[1];
rz(1.6621236) q[2];
sx q[2];
rz(-1.0589561) q[2];
sx q[2];
rz(2.3607766) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9069104) q[1];
sx q[1];
rz(-2.249472) q[1];
sx q[1];
rz(-2.3386416) q[1];
rz(-pi) q[2];
rz(1.4928905) q[3];
sx q[3];
rz(-0.39038218) q[3];
sx q[3];
rz(-1.4161033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46821758) q[2];
sx q[2];
rz(-0.28926352) q[2];
sx q[2];
rz(2.3667228) q[2];
rz(-1.8153927) q[3];
sx q[3];
rz(-1.9522791) q[3];
sx q[3];
rz(0.51138043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4020017) q[0];
sx q[0];
rz(-0.87257659) q[0];
sx q[0];
rz(0.032489754) q[0];
rz(1.2292817) q[1];
sx q[1];
rz(-0.928343) q[1];
sx q[1];
rz(-0.083018735) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0262526) q[0];
sx q[0];
rz(-1.3466382) q[0];
sx q[0];
rz(2.6782533) q[0];
rz(-pi) q[1];
rz(-1.4013702) q[2];
sx q[2];
rz(-1.3063141) q[2];
sx q[2];
rz(-1.5096017) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2526892) q[1];
sx q[1];
rz(-2.0344043) q[1];
sx q[1];
rz(0.75485922) q[1];
rz(-1.549785) q[3];
sx q[3];
rz(-2.5046231) q[3];
sx q[3];
rz(-1.4786428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3518389) q[2];
sx q[2];
rz(-0.80078501) q[2];
sx q[2];
rz(-0.16072533) q[2];
rz(1.1927346) q[3];
sx q[3];
rz(-0.21218097) q[3];
sx q[3];
rz(0.76782697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47650325) q[0];
sx q[0];
rz(-2.022321) q[0];
sx q[0];
rz(0.29092586) q[0];
rz(-1.7581958) q[1];
sx q[1];
rz(-1.8427126) q[1];
sx q[1];
rz(-2.6417522) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4668149) q[0];
sx q[0];
rz(-1.6143911) q[0];
sx q[0];
rz(-1.528426) q[0];
rz(-pi) q[1];
rz(-0.55576365) q[2];
sx q[2];
rz(-2.2190208) q[2];
sx q[2];
rz(-2.6157635) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9466797) q[1];
sx q[1];
rz(-0.99338594) q[1];
sx q[1];
rz(1.0275847) q[1];
x q[2];
rz(2.1736702) q[3];
sx q[3];
rz(-2.335066) q[3];
sx q[3];
rz(1.6953118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2073888) q[2];
sx q[2];
rz(-2.5298205) q[2];
sx q[2];
rz(2.440051) q[2];
rz(0.97405854) q[3];
sx q[3];
rz(-0.8166703) q[3];
sx q[3];
rz(-0.90132236) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3153673) q[0];
sx q[0];
rz(-0.81804818) q[0];
sx q[0];
rz(2.4819964) q[0];
rz(-1.9024128) q[1];
sx q[1];
rz(-1.0737123) q[1];
sx q[1];
rz(1.0741796) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65339966) q[0];
sx q[0];
rz(-1.5852514) q[0];
sx q[0];
rz(-3.1281592) q[0];
rz(-pi) q[1];
rz(0.72788179) q[2];
sx q[2];
rz(-0.54833503) q[2];
sx q[2];
rz(-0.20073433) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2905137) q[1];
sx q[1];
rz(-2.004262) q[1];
sx q[1];
rz(-0.23833935) q[1];
x q[2];
rz(2.6857008) q[3];
sx q[3];
rz(-0.8178725) q[3];
sx q[3];
rz(-2.4867833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2913975) q[2];
sx q[2];
rz(-0.75866282) q[2];
sx q[2];
rz(-0.58471739) q[2];
rz(1.0662474) q[3];
sx q[3];
rz(-1.2277579) q[3];
sx q[3];
rz(-0.40759459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.19121118) q[0];
sx q[0];
rz(-1.1703015) q[0];
sx q[0];
rz(0.48072746) q[0];
rz(-1.2113384) q[1];
sx q[1];
rz(-1.5296661) q[1];
sx q[1];
rz(2.5239677) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.077972875) q[0];
sx q[0];
rz(-1.5750066) q[0];
sx q[0];
rz(1.2456889) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60725339) q[2];
sx q[2];
rz(-1.1285845) q[2];
sx q[2];
rz(-2.2970125) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.27965266) q[1];
sx q[1];
rz(-1.8649001) q[1];
sx q[1];
rz(2.9169464) q[1];
x q[2];
rz(-1.6061574) q[3];
sx q[3];
rz(-1.7628944) q[3];
sx q[3];
rz(-2.0336322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9628613) q[2];
sx q[2];
rz(-1.7469254) q[2];
sx q[2];
rz(0.89087957) q[2];
rz(-2.1122872) q[3];
sx q[3];
rz(-1.3903214) q[3];
sx q[3];
rz(-1.5503957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35128281) q[0];
sx q[0];
rz(-2.2064378) q[0];
sx q[0];
rz(-1.1736322) q[0];
rz(-1.1890746) q[1];
sx q[1];
rz(-0.74780858) q[1];
sx q[1];
rz(-0.48114166) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3404589) q[0];
sx q[0];
rz(-1.0164355) q[0];
sx q[0];
rz(-1.3075243) q[0];
rz(-pi) q[1];
rz(0.5598716) q[2];
sx q[2];
rz(-1.6563043) q[2];
sx q[2];
rz(-0.97031236) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.12901017) q[1];
sx q[1];
rz(-0.6193739) q[1];
sx q[1];
rz(-1.8257837) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2523426) q[3];
sx q[3];
rz(-2.4947531) q[3];
sx q[3];
rz(-1.3726039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9159307) q[2];
sx q[2];
rz(-1.91232) q[2];
sx q[2];
rz(1.4813598) q[2];
rz(-3.0952752) q[3];
sx q[3];
rz(-1.9527083) q[3];
sx q[3];
rz(1.9143298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7939821) q[0];
sx q[0];
rz(-1.0337669) q[0];
sx q[0];
rz(0.069742918) q[0];
rz(0.28930411) q[1];
sx q[1];
rz(-1.2846839) q[1];
sx q[1];
rz(1.0522254) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4493794) q[0];
sx q[0];
rz(-1.2468296) q[0];
sx q[0];
rz(-0.37925668) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9751189) q[2];
sx q[2];
rz(-2.2262339) q[2];
sx q[2];
rz(2.2827998) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7888513) q[1];
sx q[1];
rz(-0.53376275) q[1];
sx q[1];
rz(-0.8064005) q[1];
x q[2];
rz(1.1484365) q[3];
sx q[3];
rz(-1.1760654) q[3];
sx q[3];
rz(-1.0843866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5436486) q[2];
sx q[2];
rz(-1.1755377) q[2];
sx q[2];
rz(0.0607461) q[2];
rz(2.8525823) q[3];
sx q[3];
rz(-1.3649536) q[3];
sx q[3];
rz(-1.2750767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.65171) q[0];
sx q[0];
rz(-1.4986421) q[0];
sx q[0];
rz(-1.0810252) q[0];
rz(1.0501077) q[1];
sx q[1];
rz(-1.0625912) q[1];
sx q[1];
rz(-1.2407632) q[1];
rz(-1.1012668) q[2];
sx q[2];
rz(-1.7055837) q[2];
sx q[2];
rz(-2.3619426) q[2];
rz(2.1331187) q[3];
sx q[3];
rz(-0.49839603) q[3];
sx q[3];
rz(-2.8965542) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
