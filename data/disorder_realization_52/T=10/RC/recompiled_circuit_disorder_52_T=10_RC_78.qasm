OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.62587005) q[0];
sx q[0];
rz(6.8318879) q[0];
sx q[0];
rz(5.3988342) q[0];
rz(1.4305152) q[1];
sx q[1];
rz(-2.1880452) q[1];
sx q[1];
rz(1.5024827) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8774672) q[0];
sx q[0];
rz(-1.929951) q[0];
sx q[0];
rz(2.9496664) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0083562334) q[2];
sx q[2];
rz(-0.5651606) q[2];
sx q[2];
rz(1.9985808) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1391746) q[1];
sx q[1];
rz(-2.8867509) q[1];
sx q[1];
rz(2.8789218) q[1];
x q[2];
rz(-1.9804269) q[3];
sx q[3];
rz(-0.75244609) q[3];
sx q[3];
rz(0.15256552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3551336) q[2];
sx q[2];
rz(-0.81374514) q[2];
sx q[2];
rz(2.4856429) q[2];
rz(1.9338699) q[3];
sx q[3];
rz(-1.175712) q[3];
sx q[3];
rz(2.1470127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8939963) q[0];
sx q[0];
rz(-0.42034724) q[0];
sx q[0];
rz(-0.43352747) q[0];
rz(2.9128089) q[1];
sx q[1];
rz(-0.42963916) q[1];
sx q[1];
rz(-0.0072335009) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1996961) q[0];
sx q[0];
rz(-1.4775839) q[0];
sx q[0];
rz(1.9641563) q[0];
rz(-pi) q[1];
rz(-0.8231926) q[2];
sx q[2];
rz(-0.41831145) q[2];
sx q[2];
rz(0.38564607) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.65456395) q[1];
sx q[1];
rz(-1.9332758) q[1];
sx q[1];
rz(1.3351721) q[1];
x q[2];
rz(2.238027) q[3];
sx q[3];
rz(-1.9575319) q[3];
sx q[3];
rz(-0.29153338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6146415) q[2];
sx q[2];
rz(-0.80792892) q[2];
sx q[2];
rz(-0.69765222) q[2];
rz(0.12156045) q[3];
sx q[3];
rz(-1.9024885) q[3];
sx q[3];
rz(2.8377623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-1.6737297) q[0];
sx q[0];
rz(-2.2255852) q[0];
sx q[0];
rz(-1.3695705) q[0];
rz(1.9000152) q[1];
sx q[1];
rz(-1.4135655) q[1];
sx q[1];
rz(0.26161584) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81944377) q[0];
sx q[0];
rz(-1.5804287) q[0];
sx q[0];
rz(1.5887512) q[0];
rz(3.0171266) q[2];
sx q[2];
rz(-2.3111768) q[2];
sx q[2];
rz(-0.24701961) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2178104) q[1];
sx q[1];
rz(-0.79499704) q[1];
sx q[1];
rz(1.7407106) q[1];
rz(-2.6286078) q[3];
sx q[3];
rz(-1.5981711) q[3];
sx q[3];
rz(0.56752906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6802784) q[2];
sx q[2];
rz(-1.6363982) q[2];
sx q[2];
rz(0.20351163) q[2];
rz(-0.92173785) q[3];
sx q[3];
rz(-1.8739871) q[3];
sx q[3];
rz(0.27954277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1638284) q[0];
sx q[0];
rz(-1.5777359) q[0];
sx q[0];
rz(1.571636) q[0];
rz(1.0034026) q[1];
sx q[1];
rz(-1.827821) q[1];
sx q[1];
rz(-1.8932231) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.473339) q[0];
sx q[0];
rz(-3.0248397) q[0];
sx q[0];
rz(0.97572414) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3893045) q[2];
sx q[2];
rz(-1.3647807) q[2];
sx q[2];
rz(0.69603053) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6971671) q[1];
sx q[1];
rz(-0.66272347) q[1];
sx q[1];
rz(-0.18130937) q[1];
rz(-pi) q[2];
rz(-1.8321886) q[3];
sx q[3];
rz(-2.2769615) q[3];
sx q[3];
rz(-0.9725001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1127597) q[2];
sx q[2];
rz(-1.3112105) q[2];
sx q[2];
rz(2.3045585) q[2];
rz(-1.2083496) q[3];
sx q[3];
rz(-1.874606) q[3];
sx q[3];
rz(-2.3560431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0060624881) q[0];
sx q[0];
rz(-2.0879789) q[0];
sx q[0];
rz(2.3663882) q[0];
rz(-2.7397621) q[1];
sx q[1];
rz(-0.95087516) q[1];
sx q[1];
rz(-2.2391589) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31061253) q[0];
sx q[0];
rz(-1.6164301) q[0];
sx q[0];
rz(-2.5020585) q[0];
rz(-pi) q[1];
rz(2.6685153) q[2];
sx q[2];
rz(-0.50191754) q[2];
sx q[2];
rz(-2.9830473) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.202995) q[1];
sx q[1];
rz(-1.0293759) q[1];
sx q[1];
rz(0.90403892) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0028354672) q[3];
sx q[3];
rz(-1.2540725) q[3];
sx q[3];
rz(-1.0032652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.19501413) q[2];
sx q[2];
rz(-0.48823753) q[2];
sx q[2];
rz(1.1966594) q[2];
rz(1.442391) q[3];
sx q[3];
rz(-1.2714352) q[3];
sx q[3];
rz(1.4590013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.8114132) q[0];
sx q[0];
rz(-2.9646962) q[0];
sx q[0];
rz(2.6384171) q[0];
rz(1.4563837) q[1];
sx q[1];
rz(-2.0676985) q[1];
sx q[1];
rz(0.20176372) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70226442) q[0];
sx q[0];
rz(-1.6730671) q[0];
sx q[0];
rz(-0.066145397) q[0];
rz(-pi) q[1];
rz(2.6642338) q[2];
sx q[2];
rz(-1.9743894) q[2];
sx q[2];
rz(-1.3178283) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.7355431) q[1];
sx q[1];
rz(-0.65945259) q[1];
sx q[1];
rz(0.72592782) q[1];
x q[2];
rz(2.6334409) q[3];
sx q[3];
rz(-2.3414632) q[3];
sx q[3];
rz(-1.8964163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9266944) q[2];
sx q[2];
rz(-1.639067) q[2];
sx q[2];
rz(-2.8743437) q[2];
rz(-0.8231419) q[3];
sx q[3];
rz(-3.0624793) q[3];
sx q[3];
rz(-2.2657623) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8478407) q[0];
sx q[0];
rz(-0.1593312) q[0];
sx q[0];
rz(0.057549495) q[0];
rz(1.6607704) q[1];
sx q[1];
rz(-1.7819504) q[1];
sx q[1];
rz(2.1988791) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1377624) q[0];
sx q[0];
rz(-0.98235213) q[0];
sx q[0];
rz(-1.0914735) q[0];
x q[1];
rz(-0.16072388) q[2];
sx q[2];
rz(-1.3458999) q[2];
sx q[2];
rz(-0.42663923) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.47989935) q[1];
sx q[1];
rz(-2.6091895) q[1];
sx q[1];
rz(0.9049306) q[1];
x q[2];
rz(-1.058217) q[3];
sx q[3];
rz(-1.3430542) q[3];
sx q[3];
rz(-2.3678569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0223579) q[2];
sx q[2];
rz(-1.0417577) q[2];
sx q[2];
rz(0.38267246) q[2];
rz(-1.0391957) q[3];
sx q[3];
rz(-1.305205) q[3];
sx q[3];
rz(0.97810811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.4246178) q[0];
sx q[0];
rz(-3.0452403) q[0];
sx q[0];
rz(2.8714645) q[0];
rz(0.62942901) q[1];
sx q[1];
rz(-2.4286178) q[1];
sx q[1];
rz(2.8576635) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89284183) q[0];
sx q[0];
rz(-2.1584956) q[0];
sx q[0];
rz(-2.6177004) q[0];
x q[1];
rz(0.59992744) q[2];
sx q[2];
rz(-1.4495965) q[2];
sx q[2];
rz(1.3231414) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8995081) q[1];
sx q[1];
rz(-2.429932) q[1];
sx q[1];
rz(-2.3942024) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5408526) q[3];
sx q[3];
rz(-0.91120126) q[3];
sx q[3];
rz(-2.3823882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5876864) q[2];
sx q[2];
rz(-0.17051414) q[2];
sx q[2];
rz(-1.930687) q[2];
rz(-2.8816913) q[3];
sx q[3];
rz(-0.62429684) q[3];
sx q[3];
rz(-0.62817812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0652086) q[0];
sx q[0];
rz(-0.56088352) q[0];
sx q[0];
rz(-0.22924766) q[0];
rz(0.30300888) q[1];
sx q[1];
rz(-1.3906994) q[1];
sx q[1];
rz(1.4607666) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54366771) q[0];
sx q[0];
rz(-0.28408465) q[0];
sx q[0];
rz(-3.0790867) q[0];
rz(-pi) q[1];
rz(-0.42713366) q[2];
sx q[2];
rz(-0.41296994) q[2];
sx q[2];
rz(-2.4954748) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8742121) q[1];
sx q[1];
rz(-1.0460209) q[1];
sx q[1];
rz(1.5449779) q[1];
x q[2];
rz(1.1838412) q[3];
sx q[3];
rz(-2.317252) q[3];
sx q[3];
rz(1.4617621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71172697) q[2];
sx q[2];
rz(-1.9248328) q[2];
sx q[2];
rz(2.4460068) q[2];
rz(-0.43186489) q[3];
sx q[3];
rz(-2.6769107) q[3];
sx q[3];
rz(0.7152043) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5678976) q[0];
sx q[0];
rz(-1.2117813) q[0];
sx q[0];
rz(0.33690548) q[0];
rz(-0.20740549) q[1];
sx q[1];
rz(-1.0284871) q[1];
sx q[1];
rz(-0.38063231) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5779553) q[0];
sx q[0];
rz(-1.5008238) q[0];
sx q[0];
rz(1.7726583) q[0];
rz(-pi) q[1];
rz(0.74110465) q[2];
sx q[2];
rz(-2.0217102) q[2];
sx q[2];
rz(1.8795183) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8457348) q[1];
sx q[1];
rz(-2.0669792) q[1];
sx q[1];
rz(-0.22175281) q[1];
x q[2];
rz(-0.59488876) q[3];
sx q[3];
rz(-0.66685646) q[3];
sx q[3];
rz(0.24645933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.0011065817) q[2];
sx q[2];
rz(-1.1662741) q[2];
sx q[2];
rz(-1.1364737) q[2];
rz(3.100637) q[3];
sx q[3];
rz(-0.80871964) q[3];
sx q[3];
rz(1.827318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8052335) q[0];
sx q[0];
rz(-2.2362066) q[0];
sx q[0];
rz(2.8295828) q[0];
rz(-2.1144755) q[1];
sx q[1];
rz(-1.849091) q[1];
sx q[1];
rz(-1.0277933) q[1];
rz(0.74577352) q[2];
sx q[2];
rz(-0.33200982) q[2];
sx q[2];
rz(0.40340323) q[2];
rz(1.9789226) q[3];
sx q[3];
rz(-0.47331953) q[3];
sx q[3];
rz(2.4389653) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
