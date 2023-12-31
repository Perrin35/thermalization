OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.0032089631) q[0];
sx q[0];
rz(-0.15455833) q[0];
sx q[0];
rz(0.69252339) q[0];
rz(1.9321631) q[1];
sx q[1];
rz(-1.2485319) q[1];
sx q[1];
rz(1.7564397) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9557578) q[0];
sx q[0];
rz(-1.8896566) q[0];
sx q[0];
rz(2.3715109) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55180438) q[2];
sx q[2];
rz(-1.3349512) q[2];
sx q[2];
rz(-2.9896196) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3424211) q[1];
sx q[1];
rz(-0.28563269) q[1];
sx q[1];
rz(-0.61113961) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6119421) q[3];
sx q[3];
rz(-0.7512593) q[3];
sx q[3];
rz(0.9179759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2549071) q[2];
sx q[2];
rz(-0.79780769) q[2];
sx q[2];
rz(-2.936426) q[2];
rz(0.77130476) q[3];
sx q[3];
rz(-2.3588534) q[3];
sx q[3];
rz(1.1024968) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7339864) q[0];
sx q[0];
rz(-0.74626958) q[0];
sx q[0];
rz(-0.45390391) q[0];
rz(-2.1167963) q[1];
sx q[1];
rz(-2.7339934) q[1];
sx q[1];
rz(1.227238) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5204029) q[0];
sx q[0];
rz(-1.495201) q[0];
sx q[0];
rz(1.441342) q[0];
rz(-1.111755) q[2];
sx q[2];
rz(-1.0978205) q[2];
sx q[2];
rz(-0.74726653) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2482359) q[1];
sx q[1];
rz(-1.8289939) q[1];
sx q[1];
rz(2.8398819) q[1];
x q[2];
rz(2.482588) q[3];
sx q[3];
rz(-0.25203029) q[3];
sx q[3];
rz(0.27740955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1002905) q[2];
sx q[2];
rz(-1.1854478) q[2];
sx q[2];
rz(-0.56742898) q[2];
rz(-2.7764017) q[3];
sx q[3];
rz(-1.4130211) q[3];
sx q[3];
rz(-0.96810961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48297468) q[0];
sx q[0];
rz(-2.5768319) q[0];
sx q[0];
rz(-0.89865249) q[0];
rz(-0.99575106) q[1];
sx q[1];
rz(-1.5834705) q[1];
sx q[1];
rz(0.333289) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5211398) q[0];
sx q[0];
rz(-0.95690173) q[0];
sx q[0];
rz(-0.99434538) q[0];
rz(-pi) q[1];
rz(-1.0684418) q[2];
sx q[2];
rz(-0.2012673) q[2];
sx q[2];
rz(-1.7114491) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7283199) q[1];
sx q[1];
rz(-1.7227168) q[1];
sx q[1];
rz(2.1876213) q[1];
x q[2];
rz(1.0201449) q[3];
sx q[3];
rz(-1.2216976) q[3];
sx q[3];
rz(0.89494866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4553392) q[2];
sx q[2];
rz(-1.3505961) q[2];
sx q[2];
rz(-1.9906445) q[2];
rz(2.3006556) q[3];
sx q[3];
rz(-1.1154113) q[3];
sx q[3];
rz(1.1545198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6999917) q[0];
sx q[0];
rz(-1.561152) q[0];
sx q[0];
rz(0.68471318) q[0];
rz(2.1060064) q[1];
sx q[1];
rz(-0.5077478) q[1];
sx q[1];
rz(-1.205014) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.372615) q[0];
sx q[0];
rz(-2.2616771) q[0];
sx q[0];
rz(-3.0983503) q[0];
x q[1];
rz(-0.20385216) q[2];
sx q[2];
rz(-0.9365558) q[2];
sx q[2];
rz(-2.7472592) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1531096) q[1];
sx q[1];
rz(-1.8433237) q[1];
sx q[1];
rz(-0.86221282) q[1];
x q[2];
rz(-1.6550001) q[3];
sx q[3];
rz(-0.70662543) q[3];
sx q[3];
rz(3.1228309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.24923199) q[2];
sx q[2];
rz(-1.7148596) q[2];
sx q[2];
rz(-2.7704346) q[2];
rz(-1.4012198) q[3];
sx q[3];
rz(-0.6597844) q[3];
sx q[3];
rz(-2.0223117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.086833) q[0];
sx q[0];
rz(-2.355447) q[0];
sx q[0];
rz(0.13312419) q[0];
rz(2.1482824) q[1];
sx q[1];
rz(-1.3860093) q[1];
sx q[1];
rz(-2.5865119) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7147987) q[0];
sx q[0];
rz(-1.5580651) q[0];
sx q[0];
rz(1.2554332) q[0];
rz(-pi) q[1];
x q[1];
rz(1.416989) q[2];
sx q[2];
rz(-0.4193192) q[2];
sx q[2];
rz(0.16659444) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.75025573) q[1];
sx q[1];
rz(-1.5899982) q[1];
sx q[1];
rz(-1.125543) q[1];
rz(-0.57226945) q[3];
sx q[3];
rz(-0.096147691) q[3];
sx q[3];
rz(0.26086807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.83539) q[2];
sx q[2];
rz(-1.0058879) q[2];
sx q[2];
rz(-0.13892697) q[2];
rz(2.1991918) q[3];
sx q[3];
rz(-1.5709632) q[3];
sx q[3];
rz(-0.55148235) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5979364) q[0];
sx q[0];
rz(-2.5718226) q[0];
sx q[0];
rz(-2.561835) q[0];
rz(-3.014091) q[1];
sx q[1];
rz(-1.189905) q[1];
sx q[1];
rz(1.5396083) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.471305) q[0];
sx q[0];
rz(-1.3677214) q[0];
sx q[0];
rz(-2.2905486) q[0];
rz(-pi) q[1];
rz(-0.33467218) q[2];
sx q[2];
rz(-0.13813189) q[2];
sx q[2];
rz(1.3484671) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0582038) q[1];
sx q[1];
rz(-1.1440047) q[1];
sx q[1];
rz(0.31174387) q[1];
x q[2];
rz(-1.2195915) q[3];
sx q[3];
rz(-0.71577365) q[3];
sx q[3];
rz(2.2367665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5876028) q[2];
sx q[2];
rz(-0.25045276) q[2];
sx q[2];
rz(0.26947752) q[2];
rz(-0.23412165) q[3];
sx q[3];
rz(-0.52474371) q[3];
sx q[3];
rz(0.060119303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6948029) q[0];
sx q[0];
rz(-2.5233874) q[0];
sx q[0];
rz(0.68429464) q[0];
rz(-0.11958312) q[1];
sx q[1];
rz(-1.2922829) q[1];
sx q[1];
rz(2.6228242) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9414026) q[0];
sx q[0];
rz(-1.7046283) q[0];
sx q[0];
rz(0.83394136) q[0];
rz(2.9154645) q[2];
sx q[2];
rz(-2.3580708) q[2];
sx q[2];
rz(2.4353611) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8956633) q[1];
sx q[1];
rz(-0.83918011) q[1];
sx q[1];
rz(-0.64114665) q[1];
rz(0.026168907) q[3];
sx q[3];
rz(-1.1702288) q[3];
sx q[3];
rz(0.27163423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8873022) q[2];
sx q[2];
rz(-1.7742949) q[2];
sx q[2];
rz(0.56345144) q[2];
rz(-3.0900132) q[3];
sx q[3];
rz(-2.0023465) q[3];
sx q[3];
rz(-0.95190597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1812487) q[0];
sx q[0];
rz(-0.5287756) q[0];
sx q[0];
rz(-1.3990336) q[0];
rz(0.78701204) q[1];
sx q[1];
rz(-2.0128638) q[1];
sx q[1];
rz(-0.74434892) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23119152) q[0];
sx q[0];
rz(-1.7185128) q[0];
sx q[0];
rz(1.5157248) q[0];
x q[1];
rz(2.0853945) q[2];
sx q[2];
rz(-1.7040164) q[2];
sx q[2];
rz(2.104987) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5163755) q[1];
sx q[1];
rz(-0.55621925) q[1];
sx q[1];
rz(0.43362995) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0934422) q[3];
sx q[3];
rz(-0.45504323) q[3];
sx q[3];
rz(-0.20850785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7156334) q[2];
sx q[2];
rz(-1.3954433) q[2];
sx q[2];
rz(-2.5320833) q[2];
rz(-0.65731796) q[3];
sx q[3];
rz(-2.4980563) q[3];
sx q[3];
rz(2.8801584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9534) q[0];
sx q[0];
rz(-0.094390079) q[0];
sx q[0];
rz(-1.5040065) q[0];
rz(1.2414744) q[1];
sx q[1];
rz(-1.991661) q[1];
sx q[1];
rz(2.3666568) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0477714) q[0];
sx q[0];
rz(-1.1822961) q[0];
sx q[0];
rz(0.85431487) q[0];
x q[1];
rz(2.9853285) q[2];
sx q[2];
rz(-2.0093577) q[2];
sx q[2];
rz(-1.4807448) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1800268) q[1];
sx q[1];
rz(-1.119009) q[1];
sx q[1];
rz(-1.3614484) q[1];
x q[2];
rz(-1.7899412) q[3];
sx q[3];
rz(-1.431576) q[3];
sx q[3];
rz(-2.5185891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.187414) q[2];
sx q[2];
rz(-0.22970197) q[2];
sx q[2];
rz(-2.9837218) q[2];
rz(1.9291417) q[3];
sx q[3];
rz(-1.0586497) q[3];
sx q[3];
rz(-1.2780179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.071844) q[0];
sx q[0];
rz(-0.97244111) q[0];
sx q[0];
rz(0.21433314) q[0];
rz(2.4841323) q[1];
sx q[1];
rz(-2.9174556) q[1];
sx q[1];
rz(2.0956031) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8612254) q[0];
sx q[0];
rz(-1.4645637) q[0];
sx q[0];
rz(2.7809814) q[0];
rz(-3.0984512) q[2];
sx q[2];
rz(-2.5501745) q[2];
sx q[2];
rz(-0.45515781) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6194832) q[1];
sx q[1];
rz(-1.2088747) q[1];
sx q[1];
rz(2.1543909) q[1];
rz(-pi) q[2];
rz(0.06185992) q[3];
sx q[3];
rz(-2.7354771) q[3];
sx q[3];
rz(-1.5554242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7252698) q[2];
sx q[2];
rz(-0.060083397) q[2];
sx q[2];
rz(2.0521169) q[2];
rz(-1.5754835) q[3];
sx q[3];
rz(-1.8982866) q[3];
sx q[3];
rz(2.4889448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2789223) q[0];
sx q[0];
rz(-2.537732) q[0];
sx q[0];
rz(-2.296007) q[0];
rz(1.5325585) q[1];
sx q[1];
rz(-1.6747723) q[1];
sx q[1];
rz(2.0369045) q[1];
rz(0.63411843) q[2];
sx q[2];
rz(-0.49370439) q[2];
sx q[2];
rz(1.6903071) q[2];
rz(2.5915626) q[3];
sx q[3];
rz(-1.6246272) q[3];
sx q[3];
rz(0.18679242) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
