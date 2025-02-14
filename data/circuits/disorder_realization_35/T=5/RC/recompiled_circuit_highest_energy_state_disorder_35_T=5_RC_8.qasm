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
rz(-3.0776032) q[0];
sx q[0];
rz(-0.912323) q[0];
sx q[0];
rz(1.3349226) q[0];
rz(-1.1420684) q[1];
sx q[1];
rz(-0.51841441) q[1];
sx q[1];
rz(-1.4919182) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27150422) q[0];
sx q[0];
rz(-0.056500204) q[0];
sx q[0];
rz(1.9046049) q[0];
rz(0.9426192) q[2];
sx q[2];
rz(-0.97731263) q[2];
sx q[2];
rz(-1.3243937) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.838249) q[1];
sx q[1];
rz(-2.456899) q[1];
sx q[1];
rz(1.1176904) q[1];
rz(-pi) q[2];
rz(0.033871058) q[3];
sx q[3];
rz(-1.813214) q[3];
sx q[3];
rz(-2.9592379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3492744) q[2];
sx q[2];
rz(-2.8742542) q[2];
sx q[2];
rz(3.086669) q[2];
rz(-1.6867636) q[3];
sx q[3];
rz(-2.7456386) q[3];
sx q[3];
rz(-1.5168064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2199663) q[0];
sx q[0];
rz(-1.3082137) q[0];
sx q[0];
rz(-0.24800214) q[0];
rz(1.2451046) q[1];
sx q[1];
rz(-0.73687941) q[1];
sx q[1];
rz(-1.9128333) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1708536) q[0];
sx q[0];
rz(-1.6038994) q[0];
sx q[0];
rz(2.956361) q[0];
rz(2.6284559) q[2];
sx q[2];
rz(-2.0433132) q[2];
sx q[2];
rz(0.11829157) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.23763021) q[1];
sx q[1];
rz(-1.7824689) q[1];
sx q[1];
rz(-2.4304076) q[1];
x q[2];
rz(0.13350962) q[3];
sx q[3];
rz(-2.7179681) q[3];
sx q[3];
rz(2.6299919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2174786) q[2];
sx q[2];
rz(-2.3096297) q[2];
sx q[2];
rz(-0.58383101) q[2];
rz(-2.7122688) q[3];
sx q[3];
rz(-1.9177633) q[3];
sx q[3];
rz(1.0113641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9871224) q[0];
sx q[0];
rz(-0.38726375) q[0];
sx q[0];
rz(1.6744457) q[0];
rz(1.6636498) q[1];
sx q[1];
rz(-1.7324305) q[1];
sx q[1];
rz(2.0416226) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80187782) q[0];
sx q[0];
rz(-2.9168282) q[0];
sx q[0];
rz(-2.6532214) q[0];
rz(-pi) q[1];
rz(-0.87022484) q[2];
sx q[2];
rz(-0.43284518) q[2];
sx q[2];
rz(0.74898042) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4277122) q[1];
sx q[1];
rz(-1.1888479) q[1];
sx q[1];
rz(1.2101342) q[1];
rz(-pi) q[2];
rz(-2.1446225) q[3];
sx q[3];
rz(-1.1862635) q[3];
sx q[3];
rz(1.1763193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.786342) q[2];
sx q[2];
rz(-0.92266551) q[2];
sx q[2];
rz(-1.5477017) q[2];
rz(1.330438) q[3];
sx q[3];
rz(-1.3733613) q[3];
sx q[3];
rz(1.0734585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.380577) q[0];
sx q[0];
rz(-1.8208193) q[0];
sx q[0];
rz(-2.9837578) q[0];
rz(2.8997968) q[1];
sx q[1];
rz(-2.7561185) q[1];
sx q[1];
rz(0.48113021) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7173955) q[0];
sx q[0];
rz(-1.4823874) q[0];
sx q[0];
rz(0.7128977) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60637668) q[2];
sx q[2];
rz(-2.7131792) q[2];
sx q[2];
rz(1.0029175) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0193023) q[1];
sx q[1];
rz(-1.230352) q[1];
sx q[1];
rz(2.2346418) q[1];
rz(-pi) q[2];
rz(1.2838383) q[3];
sx q[3];
rz(-1.7718414) q[3];
sx q[3];
rz(1.505681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2699997) q[2];
sx q[2];
rz(-2.3303878) q[2];
sx q[2];
rz(-0.2743741) q[2];
rz(1.4659878) q[3];
sx q[3];
rz(-1.8612739) q[3];
sx q[3];
rz(-2.2082641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45193732) q[0];
sx q[0];
rz(-2.2640197) q[0];
sx q[0];
rz(1.0614606) q[0];
rz(0.44867107) q[1];
sx q[1];
rz(-0.45495382) q[1];
sx q[1];
rz(-1.4400858) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0827918) q[0];
sx q[0];
rz(-2.3661206) q[0];
sx q[0];
rz(2.916921) q[0];
x q[1];
rz(-1.2810122) q[2];
sx q[2];
rz(-2.7216879) q[2];
sx q[2];
rz(-0.29422255) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6430407) q[1];
sx q[1];
rz(-1.828755) q[1];
sx q[1];
rz(0.027679701) q[1];
rz(-0.96025709) q[3];
sx q[3];
rz(-1.6597851) q[3];
sx q[3];
rz(1.2582092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7777286) q[2];
sx q[2];
rz(-1.2960459) q[2];
sx q[2];
rz(-1.9484005) q[2];
rz(2.7423972) q[3];
sx q[3];
rz(-1.4123071) q[3];
sx q[3];
rz(-0.027755888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4611918) q[0];
sx q[0];
rz(-1.2874648) q[0];
sx q[0];
rz(-0.68049085) q[0];
rz(0.83241278) q[1];
sx q[1];
rz(-1.2544371) q[1];
sx q[1];
rz(0.15636538) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45711058) q[0];
sx q[0];
rz(-1.6231939) q[0];
sx q[0];
rz(2.8657317) q[0];
rz(0.83145468) q[2];
sx q[2];
rz(-2.105684) q[2];
sx q[2];
rz(-1.8597459) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.865327) q[1];
sx q[1];
rz(-0.4368096) q[1];
sx q[1];
rz(-2.275009) q[1];
rz(-pi) q[2];
x q[2];
rz(1.353305) q[3];
sx q[3];
rz(-1.6902349) q[3];
sx q[3];
rz(-0.34298204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.51123315) q[2];
sx q[2];
rz(-0.66143051) q[2];
sx q[2];
rz(2.2516001) q[2];
rz(-2.6651799) q[3];
sx q[3];
rz(-2.2178631) q[3];
sx q[3];
rz(-2.711003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7156859) q[0];
sx q[0];
rz(-1.3193193) q[0];
sx q[0];
rz(2.529378) q[0];
rz(-1.3587492) q[1];
sx q[1];
rz(-2.3814059) q[1];
sx q[1];
rz(2.8403958) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8071547) q[0];
sx q[0];
rz(-0.2858735) q[0];
sx q[0];
rz(-1.2561594) q[0];
rz(-pi) q[1];
rz(0.91353784) q[2];
sx q[2];
rz(-1.871167) q[2];
sx q[2];
rz(2.0180839) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0665022) q[1];
sx q[1];
rz(-1.9478223) q[1];
sx q[1];
rz(0.45555631) q[1];
x q[2];
rz(-2.4047818) q[3];
sx q[3];
rz(-0.34837803) q[3];
sx q[3];
rz(-0.37955561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.907454) q[2];
sx q[2];
rz(-0.80628866) q[2];
sx q[2];
rz(-1.0478728) q[2];
rz(-0.35510865) q[3];
sx q[3];
rz(-0.57098782) q[3];
sx q[3];
rz(0.6855489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2316786) q[0];
sx q[0];
rz(-0.34220085) q[0];
sx q[0];
rz(-0.32383305) q[0];
rz(0.58770761) q[1];
sx q[1];
rz(-0.67981845) q[1];
sx q[1];
rz(0.30276611) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97070951) q[0];
sx q[0];
rz(-1.3178749) q[0];
sx q[0];
rz(-0.33527261) q[0];
x q[1];
rz(0.23719533) q[2];
sx q[2];
rz(-1.6182119) q[2];
sx q[2];
rz(0.42979017) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1582022) q[1];
sx q[1];
rz(-1.2081971) q[1];
sx q[1];
rz(-0.98539488) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1015424) q[3];
sx q[3];
rz(-2.663739) q[3];
sx q[3];
rz(-0.32020204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4964464) q[2];
sx q[2];
rz(-2.1108997) q[2];
sx q[2];
rz(2.5808064) q[2];
rz(2.136611) q[3];
sx q[3];
rz(-1.8533665) q[3];
sx q[3];
rz(-2.0463478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1110558) q[0];
sx q[0];
rz(-0.66937864) q[0];
sx q[0];
rz(2.3916767) q[0];
rz(0.38052446) q[1];
sx q[1];
rz(-0.98697248) q[1];
sx q[1];
rz(-2.1662625) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048063979) q[0];
sx q[0];
rz(-0.87665999) q[0];
sx q[0];
rz(-2.3175236) q[0];
rz(-1.8504233) q[2];
sx q[2];
rz(-1.5131009) q[2];
sx q[2];
rz(-1.1448432) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5228032) q[1];
sx q[1];
rz(-1.0741884) q[1];
sx q[1];
rz(1.9267531) q[1];
rz(-pi) q[2];
rz(0.22244979) q[3];
sx q[3];
rz(-1.8651903) q[3];
sx q[3];
rz(0.45757145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.19899496) q[2];
sx q[2];
rz(-0.90513217) q[2];
sx q[2];
rz(-2.060176) q[2];
rz(-0.069247581) q[3];
sx q[3];
rz(-1.4360177) q[3];
sx q[3];
rz(-0.92250219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0104495) q[0];
sx q[0];
rz(-0.32587019) q[0];
sx q[0];
rz(-2.8358054) q[0];
rz(0.30793134) q[1];
sx q[1];
rz(-1.6653857) q[1];
sx q[1];
rz(1.9889132) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4534396) q[0];
sx q[0];
rz(-1.3092293) q[0];
sx q[0];
rz(-0.75193172) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2158578) q[2];
sx q[2];
rz(-1.4042428) q[2];
sx q[2];
rz(-1.8018617) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.252723) q[1];
sx q[1];
rz(-2.6111341) q[1];
sx q[1];
rz(-2.2525851) q[1];
rz(-pi) q[2];
rz(1.5074499) q[3];
sx q[3];
rz(-1.5858486) q[3];
sx q[3];
rz(1.9346332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5311188) q[2];
sx q[2];
rz(-1.0627154) q[2];
sx q[2];
rz(-1.6884241) q[2];
rz(-0.48062634) q[3];
sx q[3];
rz(-0.56699816) q[3];
sx q[3];
rz(-1.3948729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4881445) q[0];
sx q[0];
rz(-1.5339889) q[0];
sx q[0];
rz(1.4657159) q[0];
rz(-1.0427955) q[1];
sx q[1];
rz(-1.7386309) q[1];
sx q[1];
rz(2.244619) q[1];
rz(-2.7077012) q[2];
sx q[2];
rz(-0.57882751) q[2];
sx q[2];
rz(-0.71633518) q[2];
rz(-0.35814169) q[3];
sx q[3];
rz(-1.609057) q[3];
sx q[3];
rz(-1.9032794) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
