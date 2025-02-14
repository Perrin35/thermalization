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
rz(0.053057916) q[0];
sx q[0];
rz(-2.7795656) q[0];
sx q[0];
rz(1.9757353) q[0];
rz(-0.8968269) q[1];
sx q[1];
rz(-1.4520175) q[1];
sx q[1];
rz(1.4282164) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4025164) q[0];
sx q[0];
rz(-1.6310283) q[0];
sx q[0];
rz(-0.51243685) q[0];
rz(-pi) q[1];
rz(-2.8901454) q[2];
sx q[2];
rz(-1.5000952) q[2];
sx q[2];
rz(-1.9858907) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5184091) q[1];
sx q[1];
rz(-1.3819547) q[1];
sx q[1];
rz(-0.13407003) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0479508) q[3];
sx q[3];
rz(-0.53864281) q[3];
sx q[3];
rz(-0.25687309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4521788) q[2];
sx q[2];
rz(-0.016760085) q[2];
sx q[2];
rz(0.20081271) q[2];
rz(0.14874841) q[3];
sx q[3];
rz(-3.1368308) q[3];
sx q[3];
rz(-0.31140056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0018604) q[0];
sx q[0];
rz(-2.5480324) q[0];
sx q[0];
rz(-2.1019905) q[0];
rz(-0.014558583) q[1];
sx q[1];
rz(-1.9083551) q[1];
sx q[1];
rz(1.5537517) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4400875) q[0];
sx q[0];
rz(-2.0962414) q[0];
sx q[0];
rz(0.82869014) q[0];
rz(-1.7660242) q[2];
sx q[2];
rz(-3.0659817) q[2];
sx q[2];
rz(-1.4812673) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0224198) q[1];
sx q[1];
rz(-1.3117562) q[1];
sx q[1];
rz(-0.035373493) q[1];
rz(-pi) q[2];
x q[2];
rz(0.513345) q[3];
sx q[3];
rz(-0.37783315) q[3];
sx q[3];
rz(-2.6031818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5130834) q[2];
sx q[2];
rz(-1.6078948) q[2];
sx q[2];
rz(1.3879363) q[2];
rz(1.3758818) q[3];
sx q[3];
rz(-2.0949771) q[3];
sx q[3];
rz(0.29471135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7498216) q[0];
sx q[0];
rz(-2.9048558) q[0];
sx q[0];
rz(-0.6075052) q[0];
rz(-1.5982184) q[1];
sx q[1];
rz(-2.9608455) q[1];
sx q[1];
rz(0.9651331) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34244363) q[0];
sx q[0];
rz(-1.5913054) q[0];
sx q[0];
rz(-0.0040702013) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7982462) q[2];
sx q[2];
rz(-2.008524) q[2];
sx q[2];
rz(-0.9870607) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.17441347) q[1];
sx q[1];
rz(-1.426218) q[1];
sx q[1];
rz(0.036199613) q[1];
x q[2];
rz(1.4749281) q[3];
sx q[3];
rz(-1.4879259) q[3];
sx q[3];
rz(-1.3898943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8049916) q[2];
sx q[2];
rz(-0.6707297) q[2];
sx q[2];
rz(-0.86656183) q[2];
rz(-2.0349515) q[3];
sx q[3];
rz(-1.5525147) q[3];
sx q[3];
rz(-1.470587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.57257819) q[0];
sx q[0];
rz(-0.67936474) q[0];
sx q[0];
rz(1.6160075) q[0];
rz(-0.010604803) q[1];
sx q[1];
rz(-0.0037825982) q[1];
sx q[1];
rz(2.3912281) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5724154) q[0];
sx q[0];
rz(-0.67052013) q[0];
sx q[0];
rz(-2.9886093) q[0];
x q[1];
rz(-1.7979787) q[2];
sx q[2];
rz(-2.1641693) q[2];
sx q[2];
rz(-0.57791711) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.039783104) q[1];
sx q[1];
rz(-2.0763811) q[1];
sx q[1];
rz(-1.702448) q[1];
rz(-pi) q[2];
rz(-0.37639002) q[3];
sx q[3];
rz(-2.808066) q[3];
sx q[3];
rz(0.1025478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41287199) q[2];
sx q[2];
rz(-2.0527288) q[2];
sx q[2];
rz(1.2415761) q[2];
rz(-0.0056754644) q[3];
sx q[3];
rz(-2.3311876) q[3];
sx q[3];
rz(-0.84398213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4960957) q[0];
sx q[0];
rz(-0.050300751) q[0];
sx q[0];
rz(2.2182933) q[0];
rz(0.80054379) q[1];
sx q[1];
rz(-0.0034595483) q[1];
sx q[1];
rz(-0.18005767) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31949735) q[0];
sx q[0];
rz(-1.8094157) q[0];
sx q[0];
rz(-0.042763189) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6864941) q[2];
sx q[2];
rz(-3.0149934) q[2];
sx q[2];
rz(1.1243658) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6392526) q[1];
sx q[1];
rz(-2.709734) q[1];
sx q[1];
rz(0.60245241) q[1];
rz(-3.0308444) q[3];
sx q[3];
rz(-1.2071273) q[3];
sx q[3];
rz(1.8520385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8397612) q[2];
sx q[2];
rz(-1.3000891) q[2];
sx q[2];
rz(-1.4361471) q[2];
rz(1.885421) q[3];
sx q[3];
rz(-1.4830282) q[3];
sx q[3];
rz(3.0459611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65211463) q[0];
sx q[0];
rz(-0.57179946) q[0];
sx q[0];
rz(-0.44334626) q[0];
rz(1.9709142) q[1];
sx q[1];
rz(-3.1405293) q[1];
sx q[1];
rz(2.7098157) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74878377) q[0];
sx q[0];
rz(-1.8816299) q[0];
sx q[0];
rz(2.3883467) q[0];
rz(-pi) q[1];
rz(-1.121793) q[2];
sx q[2];
rz(-0.19059715) q[2];
sx q[2];
rz(-0.42241851) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1237632) q[1];
sx q[1];
rz(-2.4430954) q[1];
sx q[1];
rz(-2.5269954) q[1];
x q[2];
rz(-1.3408991) q[3];
sx q[3];
rz(-1.2480924) q[3];
sx q[3];
rz(1.6663807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.40054896) q[2];
sx q[2];
rz(-2.2811175) q[2];
sx q[2];
rz(1.384037) q[2];
rz(-2.4436229) q[3];
sx q[3];
rz(-2.3845086) q[3];
sx q[3];
rz(-2.82011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8292002) q[0];
sx q[0];
rz(-1.3425403) q[0];
sx q[0];
rz(0.37305748) q[0];
rz(2.864605) q[1];
sx q[1];
rz(-3.1412558) q[1];
sx q[1];
rz(0.75609797) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6269659) q[0];
sx q[0];
rz(-1.9494162) q[0];
sx q[0];
rz(-2.4832804) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2038846) q[2];
sx q[2];
rz(-1.1808625) q[2];
sx q[2];
rz(-0.45118794) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.71891975) q[1];
sx q[1];
rz(-0.70168272) q[1];
sx q[1];
rz(2.2223081) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1499665) q[3];
sx q[3];
rz(-1.2927353) q[3];
sx q[3];
rz(-2.7476877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.21195352) q[2];
sx q[2];
rz(-2.5888207) q[2];
sx q[2];
rz(-0.92145222) q[2];
rz(-0.067342162) q[3];
sx q[3];
rz(-1.208297) q[3];
sx q[3];
rz(-1.6578081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9901554) q[0];
sx q[0];
rz(-0.28011265) q[0];
sx q[0];
rz(0.13023278) q[0];
rz(2.3479334) q[1];
sx q[1];
rz(-0.001611324) q[1];
sx q[1];
rz(-2.8614955) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5406427) q[0];
sx q[0];
rz(-1.4919992) q[0];
sx q[0];
rz(3.0879091) q[0];
rz(-pi) q[1];
rz(-0.38298265) q[2];
sx q[2];
rz(-2.9071701) q[2];
sx q[2];
rz(1.2146666) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.51378) q[1];
sx q[1];
rz(-2.4901721) q[1];
sx q[1];
rz(2.2808204) q[1];
x q[2];
rz(0.71106588) q[3];
sx q[3];
rz(-1.1991522) q[3];
sx q[3];
rz(1.7106595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0555931) q[2];
sx q[2];
rz(-1.6722101) q[2];
sx q[2];
rz(-0.80297339) q[2];
rz(1.5351852) q[3];
sx q[3];
rz(-2.1740422) q[3];
sx q[3];
rz(-2.3101961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0292173) q[0];
sx q[0];
rz(-3.1387098) q[0];
sx q[0];
rz(3.032384) q[0];
rz(-2.7681007) q[1];
sx q[1];
rz(-1.9544574) q[1];
sx q[1];
rz(0.55140299) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63579381) q[0];
sx q[0];
rz(-0.83573816) q[0];
sx q[0];
rz(-0.73979017) q[0];
x q[1];
rz(-1.4457747) q[2];
sx q[2];
rz(-1.3879965) q[2];
sx q[2];
rz(-1.3318544) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.99433339) q[1];
sx q[1];
rz(-2.2698987) q[1];
sx q[1];
rz(-2.2963013) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53111303) q[3];
sx q[3];
rz(-0.45980849) q[3];
sx q[3];
rz(-2.8961033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6138844) q[2];
sx q[2];
rz(-1.3091427) q[2];
sx q[2];
rz(1.8180397) q[2];
rz(-1.8771133) q[3];
sx q[3];
rz(-1.8529961) q[3];
sx q[3];
rz(0.0059787353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-1.5962113) q[0];
sx q[0];
rz(-2.5061506) q[0];
sx q[0];
rz(-0.7363466) q[0];
rz(2.9296854) q[1];
sx q[1];
rz(-1.0286464) q[1];
sx q[1];
rz(1.5440936) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6074267) q[0];
sx q[0];
rz(-1.5887999) q[0];
sx q[0];
rz(-2.8624363) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0545636) q[2];
sx q[2];
rz(-0.030478796) q[2];
sx q[2];
rz(1.4914184) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5888583) q[1];
sx q[1];
rz(-0.97179669) q[1];
sx q[1];
rz(2.1373279) q[1];
rz(-pi) q[2];
rz(-2.5251901) q[3];
sx q[3];
rz(-0.52866259) q[3];
sx q[3];
rz(-1.5746436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.836901) q[2];
sx q[2];
rz(-2.3062314) q[2];
sx q[2];
rz(1.8677853) q[2];
rz(1.6972208) q[3];
sx q[3];
rz(-0.074051753) q[3];
sx q[3];
rz(1.4950289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.973751) q[0];
sx q[0];
rz(-1.583562) q[0];
sx q[0];
rz(-1.2927443) q[0];
rz(1.5372859) q[1];
sx q[1];
rz(-0.91265408) q[1];
sx q[1];
rz(0.18462054) q[1];
rz(0.22278788) q[2];
sx q[2];
rz(-3.1004799) q[2];
sx q[2];
rz(-2.0903175) q[2];
rz(1.7459695) q[3];
sx q[3];
rz(-0.61476954) q[3];
sx q[3];
rz(0.7634306) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
