OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.96149093) q[0];
sx q[0];
rz(3.7082727) q[0];
sx q[0];
rz(10.215496) q[0];
rz(2.2136731) q[1];
sx q[1];
rz(3.1766422) q[1];
sx q[1];
rz(9.1215134) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1643675) q[0];
sx q[0];
rz(-0.70500716) q[0];
sx q[0];
rz(0.35568072) q[0];
rz(0.35776414) q[2];
sx q[2];
rz(-1.2105889) q[2];
sx q[2];
rz(-1.8829164) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3584663) q[1];
sx q[1];
rz(-2.3975961) q[1];
sx q[1];
rz(2.3143155) q[1];
x q[2];
rz(1.6536775) q[3];
sx q[3];
rz(-2.304416) q[3];
sx q[3];
rz(2.8745911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7683893) q[2];
sx q[2];
rz(-0.36036569) q[2];
sx q[2];
rz(-2.4599794) q[2];
rz(2.5810589) q[3];
sx q[3];
rz(-0.78442854) q[3];
sx q[3];
rz(0.64422977) q[3];
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
rz(0.19175567) q[0];
sx q[0];
rz(-0.90832174) q[0];
sx q[0];
rz(2.7857696) q[0];
rz(-0.25335723) q[1];
sx q[1];
rz(-2.2919877) q[1];
sx q[1];
rz(-1.8960948) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1278763) q[0];
sx q[0];
rz(-2.6111013) q[0];
sx q[0];
rz(-0.45376221) q[0];
rz(-2.5077057) q[2];
sx q[2];
rz(-2.0075463) q[2];
sx q[2];
rz(2.2417807) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.49719062) q[1];
sx q[1];
rz(-0.70919398) q[1];
sx q[1];
rz(2.4552114) q[1];
rz(-pi) q[2];
rz(-1.0222517) q[3];
sx q[3];
rz(-1.1303076) q[3];
sx q[3];
rz(-2.2817274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3972828) q[2];
sx q[2];
rz(-2.3645568) q[2];
sx q[2];
rz(-0.12766078) q[2];
rz(-2.4658261) q[3];
sx q[3];
rz(-0.45857576) q[3];
sx q[3];
rz(0.42135409) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0965939) q[0];
sx q[0];
rz(-2.8969722) q[0];
sx q[0];
rz(1.1340207) q[0];
rz(2.2484312) q[1];
sx q[1];
rz(-2.2375219) q[1];
sx q[1];
rz(-1.2718511) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.853493) q[0];
sx q[0];
rz(-1.5971759) q[0];
sx q[0];
rz(2.1335925) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8804466) q[2];
sx q[2];
rz(-1.4952588) q[2];
sx q[2];
rz(0.47020082) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0643029) q[1];
sx q[1];
rz(-1.7756808) q[1];
sx q[1];
rz(-2.2868392) q[1];
rz(-pi) q[2];
rz(2.2168733) q[3];
sx q[3];
rz(-2.8689403) q[3];
sx q[3];
rz(-0.70875185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9280055) q[2];
sx q[2];
rz(-0.86557937) q[2];
sx q[2];
rz(1.8013087) q[2];
rz(1.853893) q[3];
sx q[3];
rz(-1.6940593) q[3];
sx q[3];
rz(-2.652216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5260148) q[0];
sx q[0];
rz(-0.57681334) q[0];
sx q[0];
rz(2.4082129) q[0];
rz(0.79542696) q[1];
sx q[1];
rz(-0.19618244) q[1];
sx q[1];
rz(1.1682074) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2640894) q[0];
sx q[0];
rz(-1.4651708) q[0];
sx q[0];
rz(1.9838404) q[0];
rz(-pi) q[1];
rz(2.204037) q[2];
sx q[2];
rz(-1.0934693) q[2];
sx q[2];
rz(-0.027743159) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8454947) q[1];
sx q[1];
rz(-1.6137527) q[1];
sx q[1];
rz(-1.6011493) q[1];
x q[2];
rz(-0.9549035) q[3];
sx q[3];
rz(-2.1400937) q[3];
sx q[3];
rz(-2.3950441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1802133) q[2];
sx q[2];
rz(-0.81575477) q[2];
sx q[2];
rz(0.95480314) q[2];
rz(1.0500326) q[3];
sx q[3];
rz(-2.3539216) q[3];
sx q[3];
rz(0.55154705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0410205) q[0];
sx q[0];
rz(-0.53871483) q[0];
sx q[0];
rz(0.69514298) q[0];
rz(-1.9255385) q[1];
sx q[1];
rz(-1.0168409) q[1];
sx q[1];
rz(3.1207808) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1095525) q[0];
sx q[0];
rz(-1.4281245) q[0];
sx q[0];
rz(2.232216) q[0];
rz(0.32692744) q[2];
sx q[2];
rz(-0.86092868) q[2];
sx q[2];
rz(-0.29151379) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6823618) q[1];
sx q[1];
rz(-1.8098847) q[1];
sx q[1];
rz(3.1296631) q[1];
rz(-pi) q[2];
rz(1.5091875) q[3];
sx q[3];
rz(-1.4780668) q[3];
sx q[3];
rz(-1.217886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1702599) q[2];
sx q[2];
rz(-2.8016165) q[2];
sx q[2];
rz(0.60410947) q[2];
rz(2.4938834) q[3];
sx q[3];
rz(-0.52102399) q[3];
sx q[3];
rz(0.46085301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.0172417) q[0];
sx q[0];
rz(-1.0931953) q[0];
sx q[0];
rz(0.35853115) q[0];
rz(-0.56309807) q[1];
sx q[1];
rz(-0.9022572) q[1];
sx q[1];
rz(0.30555746) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9193076) q[0];
sx q[0];
rz(-0.74757517) q[0];
sx q[0];
rz(-0.54493746) q[0];
rz(-pi) q[1];
rz(2.9592909) q[2];
sx q[2];
rz(-1.4023702) q[2];
sx q[2];
rz(1.4252848) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2787152) q[1];
sx q[1];
rz(-0.76932615) q[1];
sx q[1];
rz(-2.3946986) q[1];
x q[2];
rz(0.73718585) q[3];
sx q[3];
rz(-0.10431001) q[3];
sx q[3];
rz(0.3462458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2629338) q[2];
sx q[2];
rz(-0.96076751) q[2];
sx q[2];
rz(0.052138694) q[2];
rz(0.42929286) q[3];
sx q[3];
rz(-1.41058) q[3];
sx q[3];
rz(0.2521635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4572064) q[0];
sx q[0];
rz(-0.52427137) q[0];
sx q[0];
rz(2.9065409) q[0];
rz(2.5382407) q[1];
sx q[1];
rz(-2.3124606) q[1];
sx q[1];
rz(-0.57650173) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3424732) q[0];
sx q[0];
rz(-2.2469889) q[0];
sx q[0];
rz(-0.66838092) q[0];
rz(-pi) q[1];
rz(-1.9095678) q[2];
sx q[2];
rz(-0.82451754) q[2];
sx q[2];
rz(-0.19245779) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.602142) q[1];
sx q[1];
rz(-1.2509402) q[1];
sx q[1];
rz(-1.1509239) q[1];
rz(-2.0527538) q[3];
sx q[3];
rz(-1.4432943) q[3];
sx q[3];
rz(-1.8021999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0779886) q[2];
sx q[2];
rz(-2.3006907) q[2];
sx q[2];
rz(1.7951175) q[2];
rz(2.6320362) q[3];
sx q[3];
rz(-1.2600803) q[3];
sx q[3];
rz(-2.8349561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6154196) q[0];
sx q[0];
rz(-1.0831447) q[0];
sx q[0];
rz(-2.7857067) q[0];
rz(-0.7016167) q[1];
sx q[1];
rz(-0.18467782) q[1];
sx q[1];
rz(-2.1772749) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6414757) q[0];
sx q[0];
rz(-1.4371514) q[0];
sx q[0];
rz(-1.1125314) q[0];
x q[1];
rz(-1.6379328) q[2];
sx q[2];
rz(-0.29342857) q[2];
sx q[2];
rz(-2.1591883) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.4533952) q[1];
sx q[1];
rz(-1.5904237) q[1];
sx q[1];
rz(2.6466796) q[1];
x q[2];
rz(-0.24941949) q[3];
sx q[3];
rz(-2.8207327) q[3];
sx q[3];
rz(-1.3414931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4926766) q[2];
sx q[2];
rz(-0.70768386) q[2];
sx q[2];
rz(-0.81563449) q[2];
rz(-2.2260769) q[3];
sx q[3];
rz(-2.2736277) q[3];
sx q[3];
rz(0.22667949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.8409519) q[0];
sx q[0];
rz(-2.9673567) q[0];
sx q[0];
rz(1.2409644) q[0];
rz(-3.0366483) q[1];
sx q[1];
rz(-2.8009156) q[1];
sx q[1];
rz(2.6822283) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7442087) q[0];
sx q[0];
rz(-2.7046015) q[0];
sx q[0];
rz(0.97399141) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89519897) q[2];
sx q[2];
rz(-0.4154993) q[2];
sx q[2];
rz(-2.7875827) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8018559) q[1];
sx q[1];
rz(-1.1366399) q[1];
sx q[1];
rz(-1.804002) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4308484) q[3];
sx q[3];
rz(-2.140518) q[3];
sx q[3];
rz(2.8960814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0910054) q[2];
sx q[2];
rz(-2.4531187) q[2];
sx q[2];
rz(-3.0509994) q[2];
rz(-0.76720864) q[3];
sx q[3];
rz(-1.4827261) q[3];
sx q[3];
rz(-1.4513133) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7200658) q[0];
sx q[0];
rz(-0.92362112) q[0];
sx q[0];
rz(1.6044755) q[0];
rz(0.9556669) q[1];
sx q[1];
rz(-1.4612528) q[1];
sx q[1];
rz(2.59424) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0091396686) q[0];
sx q[0];
rz(-1.4910772) q[0];
sx q[0];
rz(-0.0059156079) q[0];
x q[1];
rz(2.5432406) q[2];
sx q[2];
rz(-0.6171591) q[2];
sx q[2];
rz(2.0218818) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.20968425) q[1];
sx q[1];
rz(-0.80291498) q[1];
sx q[1];
rz(-1.241339) q[1];
x q[2];
rz(-0.56741107) q[3];
sx q[3];
rz(-2.3443522) q[3];
sx q[3];
rz(2.286943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0305816) q[2];
sx q[2];
rz(-0.17685282) q[2];
sx q[2];
rz(-0.64876968) q[2];
rz(-2.9099921) q[3];
sx q[3];
rz(-0.88445556) q[3];
sx q[3];
rz(-0.086656682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7975174) q[0];
sx q[0];
rz(-1.2114914) q[0];
sx q[0];
rz(-1.399566) q[0];
rz(-2.4823785) q[1];
sx q[1];
rz(-1.2792239) q[1];
sx q[1];
rz(-1.4469133) q[1];
rz(1.7805889) q[2];
sx q[2];
rz(-2.9083512) q[2];
sx q[2];
rz(-0.19610263) q[2];
rz(-1.4249887) q[3];
sx q[3];
rz(-1.5167699) q[3];
sx q[3];
rz(-1.9595893) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
