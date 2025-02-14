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
rz(2.6707668) q[0];
sx q[0];
rz(-0.45007053) q[0];
sx q[0];
rz(-0.39029628) q[0];
rz(0.14416873) q[1];
sx q[1];
rz(-1.498797) q[1];
sx q[1];
rz(2.1013451) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9594355) q[0];
sx q[0];
rz(-2.6819026) q[0];
sx q[0];
rz(-0.88031405) q[0];
rz(-2.7513946) q[2];
sx q[2];
rz(-0.39053655) q[2];
sx q[2];
rz(2.2817176) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2147602) q[1];
sx q[1];
rz(-2.4364987) q[1];
sx q[1];
rz(1.8929204) q[1];
x q[2];
rz(-0.066180996) q[3];
sx q[3];
rz(-0.77941637) q[3];
sx q[3];
rz(0.73252892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1990004) q[2];
sx q[2];
rz(-0.60574836) q[2];
sx q[2];
rz(-1.5462297) q[2];
rz(0.47510251) q[3];
sx q[3];
rz(-2.4964156) q[3];
sx q[3];
rz(-1.1554385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5763181) q[0];
sx q[0];
rz(-2.1529038) q[0];
sx q[0];
rz(-0.43462547) q[0];
rz(1.0644396) q[1];
sx q[1];
rz(-2.7480405) q[1];
sx q[1];
rz(2.2501066) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1284244) q[0];
sx q[0];
rz(-1.9697609) q[0];
sx q[0];
rz(1.9600887) q[0];
rz(-pi) q[1];
rz(-0.72618809) q[2];
sx q[2];
rz(-1.9627769) q[2];
sx q[2];
rz(1.3996982) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.75900092) q[1];
sx q[1];
rz(-1.0487952) q[1];
sx q[1];
rz(-1.0434138) q[1];
x q[2];
rz(0.91153146) q[3];
sx q[3];
rz(-1.5943822) q[3];
sx q[3];
rz(-0.12480907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6330304) q[2];
sx q[2];
rz(-0.56266251) q[2];
sx q[2];
rz(2.0821345) q[2];
rz(1.2262454) q[3];
sx q[3];
rz(-1.6112593) q[3];
sx q[3];
rz(0.1196158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0934963) q[0];
sx q[0];
rz(-2.7432848) q[0];
sx q[0];
rz(-0.68921971) q[0];
rz(2.7293909) q[1];
sx q[1];
rz(-2.0239794) q[1];
sx q[1];
rz(1.9788007) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9927558) q[0];
sx q[0];
rz(-0.34275469) q[0];
sx q[0];
rz(1.1008016) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5291205) q[2];
sx q[2];
rz(-2.2650913) q[2];
sx q[2];
rz(1.5907254) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38123576) q[1];
sx q[1];
rz(-2.3489526) q[1];
sx q[1];
rz(-2.4355678) q[1];
rz(-0.27966313) q[3];
sx q[3];
rz(-1.8496017) q[3];
sx q[3];
rz(1.5115084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.73843655) q[2];
sx q[2];
rz(-1.4102035) q[2];
sx q[2];
rz(2.9031244) q[2];
rz(2.9078935) q[3];
sx q[3];
rz(-2.3668079) q[3];
sx q[3];
rz(2.9888747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9018263) q[0];
sx q[0];
rz(-0.16401839) q[0];
sx q[0];
rz(-0.30583403) q[0];
rz(2.8726874) q[1];
sx q[1];
rz(-2.3482359) q[1];
sx q[1];
rz(2.7893524) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0687843) q[0];
sx q[0];
rz(-0.13958344) q[0];
sx q[0];
rz(-1.8603345) q[0];
rz(1.8369434) q[2];
sx q[2];
rz(-2.4983642) q[2];
sx q[2];
rz(0.51189724) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8246987) q[1];
sx q[1];
rz(-2.3815063) q[1];
sx q[1];
rz(2.555067) q[1];
x q[2];
rz(-2.8958578) q[3];
sx q[3];
rz(-2.338109) q[3];
sx q[3];
rz(-1.5218228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.29117808) q[2];
sx q[2];
rz(-1.520227) q[2];
sx q[2];
rz(-0.16695437) q[2];
rz(3.0252365) q[3];
sx q[3];
rz(-2.6226624) q[3];
sx q[3];
rz(-1.2701344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51036924) q[0];
sx q[0];
rz(-2.8041757) q[0];
sx q[0];
rz(1.1153197) q[0];
rz(-1.2315617) q[1];
sx q[1];
rz(-1.4304588) q[1];
sx q[1];
rz(-0.11428741) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5047805) q[0];
sx q[0];
rz(-2.0858602) q[0];
sx q[0];
rz(-2.5367141) q[0];
rz(-1.8744993) q[2];
sx q[2];
rz(-1.8371008) q[2];
sx q[2];
rz(1.6144621) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6553919) q[1];
sx q[1];
rz(-2.1874912) q[1];
sx q[1];
rz(1.060703) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.435569) q[3];
sx q[3];
rz(-1.5665652) q[3];
sx q[3];
rz(2.8805681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.085122434) q[2];
sx q[2];
rz(-1.9279927) q[2];
sx q[2];
rz(-1.6950133) q[2];
rz(-2.1794686) q[3];
sx q[3];
rz(-2.646793) q[3];
sx q[3];
rz(-1.3303293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98244786) q[0];
sx q[0];
rz(-1.5788989) q[0];
sx q[0];
rz(0.53318095) q[0];
rz(1.606733) q[1];
sx q[1];
rz(-0.77203647) q[1];
sx q[1];
rz(2.3981222) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0395775) q[0];
sx q[0];
rz(-1.8626889) q[0];
sx q[0];
rz(1.4881698) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23225611) q[2];
sx q[2];
rz(-2.1035367) q[2];
sx q[2];
rz(-1.7552623) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.096315) q[1];
sx q[1];
rz(-1.320134) q[1];
sx q[1];
rz(2.4947007) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1164769) q[3];
sx q[3];
rz(-0.66963306) q[3];
sx q[3];
rz(2.7463934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6712436) q[2];
sx q[2];
rz(-0.84731421) q[2];
sx q[2];
rz(0.49873763) q[2];
rz(2.4399452) q[3];
sx q[3];
rz(-2.4528153) q[3];
sx q[3];
rz(-2.4832895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43211234) q[0];
sx q[0];
rz(-1.226959) q[0];
sx q[0];
rz(3.0416601) q[0];
rz(1.2808895) q[1];
sx q[1];
rz(-2.62968) q[1];
sx q[1];
rz(2.4023712) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3456536) q[0];
sx q[0];
rz(-1.1689742) q[0];
sx q[0];
rz(1.9964735) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1110837) q[2];
sx q[2];
rz(-2.0000474) q[2];
sx q[2];
rz(2.3547518) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2120253) q[1];
sx q[1];
rz(-1.9627357) q[1];
sx q[1];
rz(2.2076616) q[1];
rz(-pi) q[2];
rz(-2.337226) q[3];
sx q[3];
rz(-1.5782246) q[3];
sx q[3];
rz(0.69350375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7382536) q[2];
sx q[2];
rz(-0.6311987) q[2];
sx q[2];
rz(-1.5367907) q[2];
rz(-1.4295476) q[3];
sx q[3];
rz(-1.6481383) q[3];
sx q[3];
rz(2.3454989) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1256063) q[0];
sx q[0];
rz(-2.7637389) q[0];
sx q[0];
rz(2.0269537) q[0];
rz(2.5573225) q[1];
sx q[1];
rz(-2.22157) q[1];
sx q[1];
rz(-3.027473) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18543359) q[0];
sx q[0];
rz(-1.2853649) q[0];
sx q[0];
rz(2.5117842) q[0];
x q[1];
rz(0.80107208) q[2];
sx q[2];
rz(-1.984904) q[2];
sx q[2];
rz(-2.0887449) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8128359) q[1];
sx q[1];
rz(-0.90531534) q[1];
sx q[1];
rz(-2.6129641) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1646284) q[3];
sx q[3];
rz(-0.28985786) q[3];
sx q[3];
rz(-1.690762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4383661) q[2];
sx q[2];
rz(-0.72202903) q[2];
sx q[2];
rz(2.2528193) q[2];
rz(1.599954) q[3];
sx q[3];
rz(-1.2069353) q[3];
sx q[3];
rz(0.82109872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.090076598) q[0];
sx q[0];
rz(-2.8335644) q[0];
sx q[0];
rz(-0.28710452) q[0];
rz(1.9376532) q[1];
sx q[1];
rz(-2.2724889) q[1];
sx q[1];
rz(-1.513011) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4925032) q[0];
sx q[0];
rz(-1.2375298) q[0];
sx q[0];
rz(-1.6514342) q[0];
x q[1];
rz(0.44782717) q[2];
sx q[2];
rz(-0.97140233) q[2];
sx q[2];
rz(-0.87132711) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4040428) q[1];
sx q[1];
rz(-2.3047631) q[1];
sx q[1];
rz(1.6307733) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.00079454409) q[3];
sx q[3];
rz(-1.5677139) q[3];
sx q[3];
rz(-2.2828988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5055351) q[2];
sx q[2];
rz(-1.9745741) q[2];
sx q[2];
rz(2.8324845) q[2];
rz(1.2067893) q[3];
sx q[3];
rz(-0.38201067) q[3];
sx q[3];
rz(-2.8822854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1768271) q[0];
sx q[0];
rz(-2.4137156) q[0];
sx q[0];
rz(0.65560174) q[0];
rz(2.5128095) q[1];
sx q[1];
rz(-1.4097593) q[1];
sx q[1];
rz(-1.7579196) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7611279) q[0];
sx q[0];
rz(-1.5232067) q[0];
sx q[0];
rz(-1.6143198) q[0];
rz(-pi) q[1];
rz(1.0760154) q[2];
sx q[2];
rz(-2.3411334) q[2];
sx q[2];
rz(0.59409522) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.261499) q[1];
sx q[1];
rz(-1.4215018) q[1];
sx q[1];
rz(1.5890676) q[1];
rz(-pi) q[2];
rz(1.4301368) q[3];
sx q[3];
rz(-1.0753618) q[3];
sx q[3];
rz(0.87498935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.76929602) q[2];
sx q[2];
rz(-1.6326222) q[2];
sx q[2];
rz(2.8995635) q[2];
rz(-2.5567143) q[3];
sx q[3];
rz(-2.4681028) q[3];
sx q[3];
rz(2.6175595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.723421) q[0];
sx q[0];
rz(-0.62234288) q[0];
sx q[0];
rz(-0.55707669) q[0];
rz(2.2611025) q[1];
sx q[1];
rz(-1.4028032) q[1];
sx q[1];
rz(1.0407851) q[1];
rz(-1.5401079) q[2];
sx q[2];
rz(-2.9353113) q[2];
sx q[2];
rz(-2.7842709) q[2];
rz(0.09481341) q[3];
sx q[3];
rz(-2.1349813) q[3];
sx q[3];
rz(2.3053942) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
