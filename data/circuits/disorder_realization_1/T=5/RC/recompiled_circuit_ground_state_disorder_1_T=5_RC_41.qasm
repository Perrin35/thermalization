OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.27788568) q[0];
sx q[0];
rz(3.570896) q[0];
sx q[0];
rz(12.255393) q[0];
rz(-1.3485981) q[1];
sx q[1];
rz(-2.0452979) q[1];
sx q[1];
rz(-0.57408339) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.047217) q[0];
sx q[0];
rz(-1.939491) q[0];
sx q[0];
rz(-2.019702) q[0];
x q[1];
rz(3.1068748) q[2];
sx q[2];
rz(-0.98811921) q[2];
sx q[2];
rz(2.5947844) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.61389) q[1];
sx q[1];
rz(-2.7655294) q[1];
sx q[1];
rz(1.3251873) q[1];
x q[2];
rz(-0.063954924) q[3];
sx q[3];
rz(-1.1950781) q[3];
sx q[3];
rz(-0.8448173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.2744039) q[2];
sx q[2];
rz(-1.8827266) q[2];
sx q[2];
rz(-1.3192419) q[2];
rz(0.073171767) q[3];
sx q[3];
rz(-1.3061482) q[3];
sx q[3];
rz(-0.81899548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0935593) q[0];
sx q[0];
rz(-1.9608542) q[0];
sx q[0];
rz(0.67980415) q[0];
rz(2.3381084) q[1];
sx q[1];
rz(-1.7796703) q[1];
sx q[1];
rz(0.63978535) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073994324) q[0];
sx q[0];
rz(-2.4129281) q[0];
sx q[0];
rz(2.0557899) q[0];
rz(-pi) q[1];
rz(2.4117208) q[2];
sx q[2];
rz(-1.4966855) q[2];
sx q[2];
rz(0.48219901) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8697193) q[1];
sx q[1];
rz(-0.98733989) q[1];
sx q[1];
rz(-0.81923749) q[1];
rz(0.30840318) q[3];
sx q[3];
rz(-1.7843912) q[3];
sx q[3];
rz(1.9540389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5378319) q[2];
sx q[2];
rz(-2.64309) q[2];
sx q[2];
rz(2.4375088) q[2];
rz(3.0294561) q[3];
sx q[3];
rz(-1.4487368) q[3];
sx q[3];
rz(-1.0224379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0391487) q[0];
sx q[0];
rz(-2.0344489) q[0];
sx q[0];
rz(2.802134) q[0];
rz(-1.1184232) q[1];
sx q[1];
rz(-2.2148841) q[1];
sx q[1];
rz(-0.035645398) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4887071) q[0];
sx q[0];
rz(-1.6566601) q[0];
sx q[0];
rz(0.20767943) q[0];
x q[1];
rz(1.8730803) q[2];
sx q[2];
rz(-1.2460104) q[2];
sx q[2];
rz(-0.12476441) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.798637) q[1];
sx q[1];
rz(-1.3561212) q[1];
sx q[1];
rz(-0.66185419) q[1];
rz(-pi) q[2];
rz(0.6142637) q[3];
sx q[3];
rz(-1.8395367) q[3];
sx q[3];
rz(0.72878557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.79225916) q[2];
sx q[2];
rz(-0.68816853) q[2];
sx q[2];
rz(-0.86000693) q[2];
rz(-2.3490014) q[3];
sx q[3];
rz(-0.20321295) q[3];
sx q[3];
rz(2.121076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0122796) q[0];
sx q[0];
rz(-1.3849994) q[0];
sx q[0];
rz(0.63283515) q[0];
rz(1.1526147) q[1];
sx q[1];
rz(-0.8747789) q[1];
sx q[1];
rz(-1.6713743) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33575422) q[0];
sx q[0];
rz(-2.4921761) q[0];
sx q[0];
rz(-1.5494256) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7196413) q[2];
sx q[2];
rz(-0.80721569) q[2];
sx q[2];
rz(2.5850353) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8006261) q[1];
sx q[1];
rz(-1.3142085) q[1];
sx q[1];
rz(-0.8534732) q[1];
rz(-pi) q[2];
rz(-0.74030627) q[3];
sx q[3];
rz(-2.3362219) q[3];
sx q[3];
rz(2.3399618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4319438) q[2];
sx q[2];
rz(-0.96082965) q[2];
sx q[2];
rz(-2.7107837) q[2];
rz(0.68239799) q[3];
sx q[3];
rz(-1.4902481) q[3];
sx q[3];
rz(-2.1916913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7948941) q[0];
sx q[0];
rz(-1.2656724) q[0];
sx q[0];
rz(-1.9516113) q[0];
rz(-0.31040141) q[1];
sx q[1];
rz(-1.0810532) q[1];
sx q[1];
rz(-0.50500542) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22689817) q[0];
sx q[0];
rz(-2.3932308) q[0];
sx q[0];
rz(-0.90004653) q[0];
x q[1];
rz(-1.1558542) q[2];
sx q[2];
rz(-1.4507074) q[2];
sx q[2];
rz(1.1121554) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0934) q[1];
sx q[1];
rz(-1.4835156) q[1];
sx q[1];
rz(2.643226) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.096432627) q[3];
sx q[3];
rz(-0.72800205) q[3];
sx q[3];
rz(-1.6681287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.40671047) q[2];
sx q[2];
rz(-2.4411185) q[2];
sx q[2];
rz(2.238838) q[2];
rz(1.0550176) q[3];
sx q[3];
rz(-2.2746634) q[3];
sx q[3];
rz(2.6796851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2257253) q[0];
sx q[0];
rz(-2.4287455) q[0];
sx q[0];
rz(-1.073904) q[0];
rz(1.3794927) q[1];
sx q[1];
rz(-2.4122489) q[1];
sx q[1];
rz(-3.0440547) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65501761) q[0];
sx q[0];
rz(-1.3897822) q[0];
sx q[0];
rz(-2.5979396) q[0];
rz(1.1647085) q[2];
sx q[2];
rz(-2.7345719) q[2];
sx q[2];
rz(2.2634552) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.40332023) q[1];
sx q[1];
rz(-3.0636859) q[1];
sx q[1];
rz(2.2408443) q[1];
rz(-pi) q[2];
rz(-1.0076526) q[3];
sx q[3];
rz(-1.1112257) q[3];
sx q[3];
rz(2.8569222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11233106) q[2];
sx q[2];
rz(-1.6882201) q[2];
sx q[2];
rz(0.41109273) q[2];
rz(1.7279651) q[3];
sx q[3];
rz(-2.3634383) q[3];
sx q[3];
rz(1.5467862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77808648) q[0];
sx q[0];
rz(-0.030487617) q[0];
sx q[0];
rz(-1.7331069) q[0];
rz(2.1259437) q[1];
sx q[1];
rz(-1.8135095) q[1];
sx q[1];
rz(1.4782864) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3625945) q[0];
sx q[0];
rz(-1.2752879) q[0];
sx q[0];
rz(-0.63531718) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4276349) q[2];
sx q[2];
rz(-0.78642008) q[2];
sx q[2];
rz(-2.8734549) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6802753) q[1];
sx q[1];
rz(-1.4948541) q[1];
sx q[1];
rz(2.0148177) q[1];
rz(-pi) q[2];
rz(1.0706581) q[3];
sx q[3];
rz(-1.4551509) q[3];
sx q[3];
rz(2.107419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1136855) q[2];
sx q[2];
rz(-2.5137641) q[2];
sx q[2];
rz(-0.17024635) q[2];
rz(1.2804735) q[3];
sx q[3];
rz(-1.6672641) q[3];
sx q[3];
rz(-1.9835651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4729507) q[0];
sx q[0];
rz(-0.96422115) q[0];
sx q[0];
rz(0.52841312) q[0];
rz(-0.93327418) q[1];
sx q[1];
rz(-2.2163138) q[1];
sx q[1];
rz(2.5859213) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4802769) q[0];
sx q[0];
rz(-0.6538891) q[0];
sx q[0];
rz(-1.2459027) q[0];
rz(-pi) q[1];
rz(3.0425983) q[2];
sx q[2];
rz(-2.3470391) q[2];
sx q[2];
rz(-2.0956479) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.96413999) q[1];
sx q[1];
rz(-2.429306) q[1];
sx q[1];
rz(-0.049055462) q[1];
x q[2];
rz(1.8818086) q[3];
sx q[3];
rz(-0.90299435) q[3];
sx q[3];
rz(-3.1329114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.8133424) q[2];
sx q[2];
rz(-2.3395061) q[2];
sx q[2];
rz(-2.8343406) q[2];
rz(0.79636374) q[3];
sx q[3];
rz(-0.73085228) q[3];
sx q[3];
rz(3.0439607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(2.5873213) q[0];
sx q[0];
rz(-2.0926496) q[0];
sx q[0];
rz(1.8121423) q[0];
rz(-2.9216271) q[1];
sx q[1];
rz(-2.797762) q[1];
sx q[1];
rz(-1.3185917) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47596915) q[0];
sx q[0];
rz(-2.7950056) q[0];
sx q[0];
rz(2.2524538) q[0];
rz(2.6418453) q[2];
sx q[2];
rz(-1.4683654) q[2];
sx q[2];
rz(2.5186755) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.73083774) q[1];
sx q[1];
rz(-1.2304092) q[1];
sx q[1];
rz(-2.0654668) q[1];
rz(-pi) q[2];
rz(2.6133282) q[3];
sx q[3];
rz(-0.79341989) q[3];
sx q[3];
rz(-1.8509282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4465176) q[2];
sx q[2];
rz(-1.7138285) q[2];
sx q[2];
rz(-1.2850777) q[2];
rz(2.255693) q[3];
sx q[3];
rz(-1.748184) q[3];
sx q[3];
rz(1.6133962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0415683) q[0];
sx q[0];
rz(-0.73606857) q[0];
sx q[0];
rz(-0.31684434) q[0];
rz(-0.44081229) q[1];
sx q[1];
rz(-2.3471954) q[1];
sx q[1];
rz(2.7712834) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0531385) q[0];
sx q[0];
rz(-1.1897108) q[0];
sx q[0];
rz(-3.1082112) q[0];
rz(0.86285186) q[2];
sx q[2];
rz(-0.86332317) q[2];
sx q[2];
rz(-1.488729) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3349975) q[1];
sx q[1];
rz(-1.8704526) q[1];
sx q[1];
rz(-0.2789168) q[1];
rz(-0.050906128) q[3];
sx q[3];
rz(-1.4671145) q[3];
sx q[3];
rz(3.0391676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4050498) q[2];
sx q[2];
rz(-2.6506347) q[2];
sx q[2];
rz(1.5092108) q[2];
rz(0.80471188) q[3];
sx q[3];
rz(-1.6377621) q[3];
sx q[3];
rz(0.10778431) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1152773) q[0];
sx q[0];
rz(-1.4125217) q[0];
sx q[0];
rz(-1.9052292) q[0];
rz(-2.9112877) q[1];
sx q[1];
rz(-0.31619148) q[1];
sx q[1];
rz(-2.2264623) q[1];
rz(-0.66441734) q[2];
sx q[2];
rz(-2.0521441) q[2];
sx q[2];
rz(-1.9569546) q[2];
rz(-0.32634278) q[3];
sx q[3];
rz(-2.1643066) q[3];
sx q[3];
rz(-1.7325919) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
