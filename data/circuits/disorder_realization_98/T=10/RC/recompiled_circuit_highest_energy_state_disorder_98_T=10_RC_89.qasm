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
rz(0.75818169) q[0];
sx q[0];
rz(-2.7739006) q[0];
sx q[0];
rz(-2.0453069) q[0];
rz(0.81049377) q[1];
sx q[1];
rz(-0.23263045) q[1];
sx q[1];
rz(-2.0356324) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6152412) q[0];
sx q[0];
rz(-2.5509074) q[0];
sx q[0];
rz(-0.75449852) q[0];
rz(-pi) q[1];
rz(2.8739086) q[2];
sx q[2];
rz(-1.6639198) q[2];
sx q[2];
rz(-0.47506079) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3333862) q[1];
sx q[1];
rz(-1.5667586) q[1];
sx q[1];
rz(1.6081078) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5777784) q[3];
sx q[3];
rz(-1.2864224) q[3];
sx q[3];
rz(-2.6740536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7652863) q[2];
sx q[2];
rz(-0.47069612) q[2];
sx q[2];
rz(0.37962309) q[2];
rz(-1.6342573) q[3];
sx q[3];
rz(-2.0416656) q[3];
sx q[3];
rz(2.1716993) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64604243) q[0];
sx q[0];
rz(-0.33691418) q[0];
sx q[0];
rz(-0.90091339) q[0];
rz(3.0100929) q[1];
sx q[1];
rz(-1.200518) q[1];
sx q[1];
rz(2.6995755) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7666047) q[0];
sx q[0];
rz(-0.33549136) q[0];
sx q[0];
rz(0.25095244) q[0];
x q[1];
rz(-0.20698787) q[2];
sx q[2];
rz(-0.0047193165) q[2];
sx q[2];
rz(-2.8568792) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0626559) q[1];
sx q[1];
rz(-2.1265944) q[1];
sx q[1];
rz(-2.1318667) q[1];
x q[2];
rz(2.7412299) q[3];
sx q[3];
rz(-1.4678848) q[3];
sx q[3];
rz(-0.77983749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3758731) q[2];
sx q[2];
rz(-1.7123545) q[2];
sx q[2];
rz(2.5326552) q[2];
rz(0.40306148) q[3];
sx q[3];
rz(-1.9654704) q[3];
sx q[3];
rz(0.5016996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32560638) q[0];
sx q[0];
rz(-0.41929647) q[0];
sx q[0];
rz(-2.0035279) q[0];
rz(-1.5886547) q[1];
sx q[1];
rz(-2.373003) q[1];
sx q[1];
rz(0.80702153) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1083384) q[0];
sx q[0];
rz(-1.6699808) q[0];
sx q[0];
rz(1.5981111) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4078232) q[2];
sx q[2];
rz(-2.0103243) q[2];
sx q[2];
rz(-1.3541612) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.096232) q[1];
sx q[1];
rz(-1.7590471) q[1];
sx q[1];
rz(-1.6922616) q[1];
x q[2];
rz(-0.93385796) q[3];
sx q[3];
rz(-2.6818706) q[3];
sx q[3];
rz(0.62593725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1300065) q[2];
sx q[2];
rz(-0.61379543) q[2];
sx q[2];
rz(1.2137132) q[2];
rz(-0.064149292) q[3];
sx q[3];
rz(-1.6545273) q[3];
sx q[3];
rz(-0.00042644342) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5478058) q[0];
sx q[0];
rz(-1.8812027) q[0];
sx q[0];
rz(0.88687819) q[0];
rz(2.9828494) q[1];
sx q[1];
rz(-2.6491149) q[1];
sx q[1];
rz(-1.278272) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87065164) q[0];
sx q[0];
rz(-2.2566911) q[0];
sx q[0];
rz(-2.565388) q[0];
rz(-0.69235574) q[2];
sx q[2];
rz(-2.468363) q[2];
sx q[2];
rz(-3.0557951) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6311196) q[1];
sx q[1];
rz(-2.2620631) q[1];
sx q[1];
rz(3.0635628) q[1];
rz(-pi) q[2];
rz(-1.3298755) q[3];
sx q[3];
rz(-0.56801012) q[3];
sx q[3];
rz(2.3532075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.13901237) q[2];
sx q[2];
rz(-0.85550344) q[2];
sx q[2];
rz(-0.95376897) q[2];
rz(1.194713) q[3];
sx q[3];
rz(-2.298893) q[3];
sx q[3];
rz(-0.4755303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3139528) q[0];
sx q[0];
rz(-0.71641818) q[0];
sx q[0];
rz(-2.5415976) q[0];
rz(0.68296877) q[1];
sx q[1];
rz(-0.78790793) q[1];
sx q[1];
rz(0.33040985) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4709188) q[0];
sx q[0];
rz(-0.74711159) q[0];
sx q[0];
rz(-1.8404191) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9783713) q[2];
sx q[2];
rz(-1.7931869) q[2];
sx q[2];
rz(1.6878273) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.323954) q[1];
sx q[1];
rz(-1.3778701) q[1];
sx q[1];
rz(2.0617538) q[1];
rz(1.0903301) q[3];
sx q[3];
rz(-0.76101979) q[3];
sx q[3];
rz(0.42945592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.43634513) q[2];
sx q[2];
rz(-1.0785495) q[2];
sx q[2];
rz(0.62502965) q[2];
rz(-2.446512) q[3];
sx q[3];
rz(-1.2806226) q[3];
sx q[3];
rz(-0.052791031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-0.92796749) q[0];
sx q[0];
rz(-1.9897505) q[0];
sx q[0];
rz(1.7684162) q[0];
rz(-1.3881418) q[1];
sx q[1];
rz(-2.2699247) q[1];
sx q[1];
rz(-1.6067827) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6944748) q[0];
sx q[0];
rz(-0.15167323) q[0];
sx q[0];
rz(1.8063074) q[0];
rz(-pi) q[1];
rz(-0.18267554) q[2];
sx q[2];
rz(-1.8790882) q[2];
sx q[2];
rz(-1.3253554) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.88253792) q[1];
sx q[1];
rz(-1.1657622) q[1];
sx q[1];
rz(-0.71230131) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0032773) q[3];
sx q[3];
rz(-1.9683016) q[3];
sx q[3];
rz(-2.60633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9390949) q[2];
sx q[2];
rz(-1.9715344) q[2];
sx q[2];
rz(1.5411752) q[2];
rz(-1.0209068) q[3];
sx q[3];
rz(-1.3991791) q[3];
sx q[3];
rz(-2.8405564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0114667) q[0];
sx q[0];
rz(-3.0045894) q[0];
sx q[0];
rz(0.89114183) q[0];
rz(0.64542422) q[1];
sx q[1];
rz(-1.3963457) q[1];
sx q[1];
rz(2.3088764) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.010983) q[0];
sx q[0];
rz(-2.0314706) q[0];
sx q[0];
rz(1.4185262) q[0];
rz(-3.0751257) q[2];
sx q[2];
rz(-2.1938213) q[2];
sx q[2];
rz(2.6084171) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.11800217) q[1];
sx q[1];
rz(-1.0806135) q[1];
sx q[1];
rz(1.3497737) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54033684) q[3];
sx q[3];
rz(-1.9129921) q[3];
sx q[3];
rz(-2.1006753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.58925313) q[2];
sx q[2];
rz(-2.4746042) q[2];
sx q[2];
rz(1.5935295) q[2];
rz(-0.42406905) q[3];
sx q[3];
rz(-1.4196906) q[3];
sx q[3];
rz(0.29485318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1770723) q[0];
sx q[0];
rz(-1.1966713) q[0];
sx q[0];
rz(-0.82343423) q[0];
rz(1.9442762) q[1];
sx q[1];
rz(-1.4875393) q[1];
sx q[1];
rz(-1.4607325) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7410008) q[0];
sx q[0];
rz(-1.680458) q[0];
sx q[0];
rz(-1.2506668) q[0];
x q[1];
rz(2.1637971) q[2];
sx q[2];
rz(-1.7109181) q[2];
sx q[2];
rz(1.8210379) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.15732652) q[1];
sx q[1];
rz(-1.6132024) q[1];
sx q[1];
rz(2.4799816) q[1];
rz(-pi) q[2];
rz(-0.36007215) q[3];
sx q[3];
rz(-1.8723893) q[3];
sx q[3];
rz(2.2693279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0494277) q[2];
sx q[2];
rz(-2.4804513) q[2];
sx q[2];
rz(-0.1235505) q[2];
rz(0.83856797) q[3];
sx q[3];
rz(-1.1127915) q[3];
sx q[3];
rz(0.53028321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0274444) q[0];
sx q[0];
rz(-2.1338978) q[0];
sx q[0];
rz(1.7281519) q[0];
rz(-0.89883262) q[1];
sx q[1];
rz(-2.2347968) q[1];
sx q[1];
rz(-2.5487505) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6956447) q[0];
sx q[0];
rz(-0.88291016) q[0];
sx q[0];
rz(-0.77151973) q[0];
x q[1];
rz(-0.18211629) q[2];
sx q[2];
rz(-1.6556634) q[2];
sx q[2];
rz(-0.1402459) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.457691) q[1];
sx q[1];
rz(-1.8967034) q[1];
sx q[1];
rz(-2.0052471) q[1];
x q[2];
rz(-2.4580703) q[3];
sx q[3];
rz(-2.0647979) q[3];
sx q[3];
rz(2.0919111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2841407) q[2];
sx q[2];
rz(-2.393674) q[2];
sx q[2];
rz(0.22135529) q[2];
rz(2.0981233) q[3];
sx q[3];
rz(-0.9404434) q[3];
sx q[3];
rz(-2.7531457) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2986044) q[0];
sx q[0];
rz(-2.8589111) q[0];
sx q[0];
rz(-0.84841949) q[0];
rz(0.09659718) q[1];
sx q[1];
rz(-1.074147) q[1];
sx q[1];
rz(2.3883147) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71614186) q[0];
sx q[0];
rz(-1.4708859) q[0];
sx q[0];
rz(1.7990756) q[0];
rz(-pi) q[1];
rz(1.5346555) q[2];
sx q[2];
rz(-2.3996764) q[2];
sx q[2];
rz(-0.20797543) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.14695653) q[1];
sx q[1];
rz(-0.6377694) q[1];
sx q[1];
rz(2.2398148) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51414689) q[3];
sx q[3];
rz(-0.78350583) q[3];
sx q[3];
rz(-1.0422106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0640556) q[2];
sx q[2];
rz(-2.3040743) q[2];
sx q[2];
rz(-2.688664) q[2];
rz(2.5405267) q[3];
sx q[3];
rz(-1.6744924) q[3];
sx q[3];
rz(2.1452904) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8375028) q[0];
sx q[0];
rz(-1.6927728) q[0];
sx q[0];
rz(0.46020831) q[0];
rz(2.9651463) q[1];
sx q[1];
rz(-2.9063168) q[1];
sx q[1];
rz(1.6520687) q[1];
rz(0.34863451) q[2];
sx q[2];
rz(-1.5713816) q[2];
sx q[2];
rz(-1.1272507) q[2];
rz(-2.413977) q[3];
sx q[3];
rz(-2.03491) q[3];
sx q[3];
rz(2.2274957) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
