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
rz(1.2443378) q[0];
sx q[0];
rz(-0.79255784) q[0];
sx q[0];
rz(-0.25294024) q[0];
rz(-0.77415544) q[1];
sx q[1];
rz(-0.3781265) q[1];
sx q[1];
rz(0.3869431) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3186424) q[0];
sx q[0];
rz(-1.3143263) q[0];
sx q[0];
rz(1.9049113) q[0];
rz(3.0091506) q[2];
sx q[2];
rz(-2.5597566) q[2];
sx q[2];
rz(1.0755607) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9305715) q[1];
sx q[1];
rz(-1.4912349) q[1];
sx q[1];
rz(0.08427517) q[1];
x q[2];
rz(0.59149318) q[3];
sx q[3];
rz(-1.9884571) q[3];
sx q[3];
rz(1.0491022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0953377) q[2];
sx q[2];
rz(-0.75901186) q[2];
sx q[2];
rz(-1.4026027) q[2];
rz(1.7272353) q[3];
sx q[3];
rz(-2.4284913) q[3];
sx q[3];
rz(-0.49899092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6801179) q[0];
sx q[0];
rz(-0.48826996) q[0];
sx q[0];
rz(2.4312191) q[0];
rz(1.0164227) q[1];
sx q[1];
rz(-1.2998394) q[1];
sx q[1];
rz(-2.3764835) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.349583) q[0];
sx q[0];
rz(-0.63711626) q[0];
sx q[0];
rz(-2.5883915) q[0];
rz(-0.76025195) q[2];
sx q[2];
rz(-2.2905802) q[2];
sx q[2];
rz(-0.30530294) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2110846) q[1];
sx q[1];
rz(-0.74568664) q[1];
sx q[1];
rz(-1.0812283) q[1];
rz(-2.9917172) q[3];
sx q[3];
rz(-1.2399763) q[3];
sx q[3];
rz(-0.56043032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0244828) q[2];
sx q[2];
rz(-2.4958002) q[2];
sx q[2];
rz(0.70466858) q[2];
rz(2.9674528) q[3];
sx q[3];
rz(-2.2691059) q[3];
sx q[3];
rz(-2.7599938) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0891377) q[0];
sx q[0];
rz(-1.8280886) q[0];
sx q[0];
rz(-2.254159) q[0];
rz(1.8916091) q[1];
sx q[1];
rz(-1.9307815) q[1];
sx q[1];
rz(-1.7807622) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1340536) q[0];
sx q[0];
rz(-2.0749395) q[0];
sx q[0];
rz(1.1372805) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.046692) q[2];
sx q[2];
rz(-2.9827318) q[2];
sx q[2];
rz(-2.9422974) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8831142) q[1];
sx q[1];
rz(-2.3732024) q[1];
sx q[1];
rz(1.7151296) q[1];
x q[2];
rz(3.0988337) q[3];
sx q[3];
rz(-1.1819289) q[3];
sx q[3];
rz(-2.1270909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7736194) q[2];
sx q[2];
rz(-2.7757288) q[2];
sx q[2];
rz(1.6486637) q[2];
rz(-0.44935539) q[3];
sx q[3];
rz(-1.8466693) q[3];
sx q[3];
rz(1.2202643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.087695) q[0];
sx q[0];
rz(-0.59006417) q[0];
sx q[0];
rz(-1.4087403) q[0];
rz(-2.159481) q[1];
sx q[1];
rz(-1.5487919) q[1];
sx q[1];
rz(1.7353479) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3228673) q[0];
sx q[0];
rz(-1.3565085) q[0];
sx q[0];
rz(-3.0279798) q[0];
rz(0.92806973) q[2];
sx q[2];
rz(-0.065970369) q[2];
sx q[2];
rz(-0.82974354) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6342369) q[1];
sx q[1];
rz(-1.0546396) q[1];
sx q[1];
rz(2.8576351) q[1];
x q[2];
rz(-0.26867433) q[3];
sx q[3];
rz(-2.2019535) q[3];
sx q[3];
rz(-0.97418601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1534319) q[2];
sx q[2];
rz(-1.016523) q[2];
sx q[2];
rz(0.15596685) q[2];
rz(0.65034741) q[3];
sx q[3];
rz(-1.1740843) q[3];
sx q[3];
rz(0.11421886) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0680189) q[0];
sx q[0];
rz(-1.2696215) q[0];
sx q[0];
rz(2.9408348) q[0];
rz(-1.1772032) q[1];
sx q[1];
rz(-2.2889844) q[1];
sx q[1];
rz(-1.1600561) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8667127) q[0];
sx q[0];
rz(-1.794941) q[0];
sx q[0];
rz(0.17205162) q[0];
rz(-pi) q[1];
rz(-1.6408553) q[2];
sx q[2];
rz(-1.4558917) q[2];
sx q[2];
rz(0.15945486) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1694035) q[1];
sx q[1];
rz(-2.6693925) q[1];
sx q[1];
rz(2.4257823) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9057285) q[3];
sx q[3];
rz(-0.55509243) q[3];
sx q[3];
rz(2.5344283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7572299) q[2];
sx q[2];
rz(-0.94567662) q[2];
sx q[2];
rz(-2.6833656) q[2];
rz(0.630817) q[3];
sx q[3];
rz(-1.9061371) q[3];
sx q[3];
rz(-2.6423776) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45796564) q[0];
sx q[0];
rz(-2.2891335) q[0];
sx q[0];
rz(-3.1347347) q[0];
rz(0.231617) q[1];
sx q[1];
rz(-1.7624785) q[1];
sx q[1];
rz(2.0793656) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7201014) q[0];
sx q[0];
rz(-1.694223) q[0];
sx q[0];
rz(-2.8272259) q[0];
rz(-pi) q[1];
rz(2.7659645) q[2];
sx q[2];
rz(-2.6115719) q[2];
sx q[2];
rz(-2.9935868) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4783027) q[1];
sx q[1];
rz(-1.8243102) q[1];
sx q[1];
rz(2.4624444) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2799479) q[3];
sx q[3];
rz(-1.6966736) q[3];
sx q[3];
rz(1.5365841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.14083938) q[2];
sx q[2];
rz(-1.2522298) q[2];
sx q[2];
rz(-1.8737277) q[2];
rz(2.2800692) q[3];
sx q[3];
rz(-1.0637161) q[3];
sx q[3];
rz(-1.7027732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2086901) q[0];
sx q[0];
rz(-2.8080495) q[0];
sx q[0];
rz(-0.21555899) q[0];
rz(-0.49628273) q[1];
sx q[1];
rz(-1.9535306) q[1];
sx q[1];
rz(-0.15942474) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27860363) q[0];
sx q[0];
rz(-1.5734324) q[0];
sx q[0];
rz(-0.9359261) q[0];
rz(-0.36783571) q[2];
sx q[2];
rz(-0.4994889) q[2];
sx q[2];
rz(-0.50423451) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.39495271) q[1];
sx q[1];
rz(-1.3356908) q[1];
sx q[1];
rz(3.0993153) q[1];
rz(2.6786118) q[3];
sx q[3];
rz(-2.0714875) q[3];
sx q[3];
rz(2.7241796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.40519199) q[2];
sx q[2];
rz(-2.0173732) q[2];
sx q[2];
rz(-2.6250725) q[2];
rz(1.2107595) q[3];
sx q[3];
rz(-1.6885933) q[3];
sx q[3];
rz(1.7621382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6571758) q[0];
sx q[0];
rz(-0.076198904) q[0];
sx q[0];
rz(0.54022378) q[0];
rz(-1.6172488) q[1];
sx q[1];
rz(-2.1397619) q[1];
sx q[1];
rz(-1.2670955) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62477797) q[0];
sx q[0];
rz(-2.5886726) q[0];
sx q[0];
rz(-2.7163158) q[0];
rz(-0.71155602) q[2];
sx q[2];
rz(-2.2268953) q[2];
sx q[2];
rz(-0.0023105769) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0896014) q[1];
sx q[1];
rz(-2.7447408) q[1];
sx q[1];
rz(2.7605935) q[1];
x q[2];
rz(-2.0319875) q[3];
sx q[3];
rz(-0.69503419) q[3];
sx q[3];
rz(-1.1170596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2544864) q[2];
sx q[2];
rz(-1.9820513) q[2];
sx q[2];
rz(2.9642588) q[2];
rz(2.0717428) q[3];
sx q[3];
rz(-2.0900574) q[3];
sx q[3];
rz(-2.9593318) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7808481) q[0];
sx q[0];
rz(-1.9875263) q[0];
sx q[0];
rz(-3.0408707) q[0];
rz(1.1489457) q[1];
sx q[1];
rz(-1.9457685) q[1];
sx q[1];
rz(-1.9675868) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74885313) q[0];
sx q[0];
rz(-2.1707188) q[0];
sx q[0];
rz(-1.3412745) q[0];
rz(-pi) q[1];
rz(1.6672544) q[2];
sx q[2];
rz(-2.3347179) q[2];
sx q[2];
rz(0.45713666) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2795978) q[1];
sx q[1];
rz(-1.5902385) q[1];
sx q[1];
rz(-1.6085546) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87971808) q[3];
sx q[3];
rz(-1.6501325) q[3];
sx q[3];
rz(-0.21896958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9775057) q[2];
sx q[2];
rz(-1.6567433) q[2];
sx q[2];
rz(-2.738415) q[2];
rz(2.4781503) q[3];
sx q[3];
rz(-0.57938975) q[3];
sx q[3];
rz(0.0042393953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73096257) q[0];
sx q[0];
rz(-2.0616489) q[0];
sx q[0];
rz(3.0626815) q[0];
rz(2.2321189) q[1];
sx q[1];
rz(-2.0957004) q[1];
sx q[1];
rz(0.39631072) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8412668) q[0];
sx q[0];
rz(-1.5864276) q[0];
sx q[0];
rz(-0.0090740694) q[0];
x q[1];
rz(-1.5780007) q[2];
sx q[2];
rz(-0.77709891) q[2];
sx q[2];
rz(-2.8030628) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4034319) q[1];
sx q[1];
rz(-2.6034174) q[1];
sx q[1];
rz(0.47026431) q[1];
rz(-pi) q[2];
rz(0.86045281) q[3];
sx q[3];
rz(-0.75569587) q[3];
sx q[3];
rz(-2.6967654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4497455) q[2];
sx q[2];
rz(-2.2578466) q[2];
sx q[2];
rz(0.79908243) q[2];
rz(-3.0633022) q[3];
sx q[3];
rz(-1.8872063) q[3];
sx q[3];
rz(0.6453132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4770724) q[0];
sx q[0];
rz(-1.399853) q[0];
sx q[0];
rz(0.54367263) q[0];
rz(-1.1935344) q[1];
sx q[1];
rz(-1.2763034) q[1];
sx q[1];
rz(-2.6313849) q[1];
rz(-1.2340679) q[2];
sx q[2];
rz(-0.48067817) q[2];
sx q[2];
rz(-0.59412063) q[2];
rz(-2.067749) q[3];
sx q[3];
rz(-0.90917773) q[3];
sx q[3];
rz(-2.3346693) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
