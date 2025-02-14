OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.6736421) q[0];
sx q[0];
rz(7.5543348) q[0];
sx q[0];
rz(13.332097) q[0];
rz(2.7453121) q[1];
sx q[1];
rz(-3.0488465) q[1];
sx q[1];
rz(-0.5527817) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39550135) q[0];
sx q[0];
rz(-1.1679018) q[0];
sx q[0];
rz(-0.42568107) q[0];
rz(-pi) q[1];
rz(-2.2596006) q[2];
sx q[2];
rz(-2.246736) q[2];
sx q[2];
rz(-1.1774398) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0243135) q[1];
sx q[1];
rz(-1.2543884) q[1];
sx q[1];
rz(-0.98543075) q[1];
rz(-0.44646937) q[3];
sx q[3];
rz(-0.95799082) q[3];
sx q[3];
rz(1.0366576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7736241) q[2];
sx q[2];
rz(-1.5017193) q[2];
sx q[2];
rz(0.088851301) q[2];
rz(-2.7446274) q[3];
sx q[3];
rz(-2.1019955) q[3];
sx q[3];
rz(1.4575492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.6615768) q[0];
sx q[0];
rz(-2.6341697) q[0];
sx q[0];
rz(-0.85451025) q[0];
rz(-0.62141934) q[1];
sx q[1];
rz(-0.3365376) q[1];
sx q[1];
rz(1.052676) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9493457) q[0];
sx q[0];
rz(-1.7070012) q[0];
sx q[0];
rz(-1.6736567) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1243049) q[2];
sx q[2];
rz(-1.948775) q[2];
sx q[2];
rz(-1.93731) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6748283) q[1];
sx q[1];
rz(-2.383814) q[1];
sx q[1];
rz(-1.8653052) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1747401) q[3];
sx q[3];
rz(-1.569442) q[3];
sx q[3];
rz(2.7905066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4738327) q[2];
sx q[2];
rz(-2.2809873) q[2];
sx q[2];
rz(0.94998002) q[2];
rz(-2.1330323) q[3];
sx q[3];
rz(-2.5086094) q[3];
sx q[3];
rz(-2.8620201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94473332) q[0];
sx q[0];
rz(-1.2161398) q[0];
sx q[0];
rz(-1.3470294) q[0];
rz(-2.0836209) q[1];
sx q[1];
rz(-2.8949013) q[1];
sx q[1];
rz(0.11014858) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5025217) q[0];
sx q[0];
rz(-0.57146954) q[0];
sx q[0];
rz(2.1854464) q[0];
rz(-3.0579849) q[2];
sx q[2];
rz(-1.7183398) q[2];
sx q[2];
rz(-1.7718441) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.55787841) q[1];
sx q[1];
rz(-1.4643351) q[1];
sx q[1];
rz(-0.6006247) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3055656) q[3];
sx q[3];
rz(-2.5213833) q[3];
sx q[3];
rz(2.7171752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9614253) q[2];
sx q[2];
rz(-1.578293) q[2];
sx q[2];
rz(-0.041672826) q[2];
rz(-2.9336119) q[3];
sx q[3];
rz(-0.22170034) q[3];
sx q[3];
rz(-0.25377932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5623077) q[0];
sx q[0];
rz(-1.8203745) q[0];
sx q[0];
rz(-0.98096171) q[0];
rz(-0.43308577) q[1];
sx q[1];
rz(-1.4739477) q[1];
sx q[1];
rz(-1.513419) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1103605) q[0];
sx q[0];
rz(-1.9680395) q[0];
sx q[0];
rz(1.375694) q[0];
rz(-pi) q[1];
rz(0.89940091) q[2];
sx q[2];
rz(-1.5601306) q[2];
sx q[2];
rz(2.991302) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.42964307) q[1];
sx q[1];
rz(-0.50772655) q[1];
sx q[1];
rz(1.3583378) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1918729) q[3];
sx q[3];
rz(-2.4746568) q[3];
sx q[3];
rz(-2.2569424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.607434) q[2];
sx q[2];
rz(-2.6223493) q[2];
sx q[2];
rz(-0.81238166) q[2];
rz(-1.8937998) q[3];
sx q[3];
rz(-1.8915853) q[3];
sx q[3];
rz(-1.8947424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8060551) q[0];
sx q[0];
rz(-2.1195109) q[0];
sx q[0];
rz(1.442765) q[0];
rz(0.45817786) q[1];
sx q[1];
rz(-1.5135601) q[1];
sx q[1];
rz(-0.75622574) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0201705) q[0];
sx q[0];
rz(-0.64682942) q[0];
sx q[0];
rz(1.8298762) q[0];
rz(2.4632023) q[2];
sx q[2];
rz(-2.3783461) q[2];
sx q[2];
rz(0.89707546) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9994558) q[1];
sx q[1];
rz(-2.9805866) q[1];
sx q[1];
rz(-0.57438897) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6709354) q[3];
sx q[3];
rz(-0.34557811) q[3];
sx q[3];
rz(0.80039961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4592287) q[2];
sx q[2];
rz(-1.7533147) q[2];
sx q[2];
rz(-0.96561042) q[2];
rz(-2.5718578) q[3];
sx q[3];
rz(-1.1252879) q[3];
sx q[3];
rz(-1.6857356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.4023034) q[0];
sx q[0];
rz(-1.9434513) q[0];
sx q[0];
rz(-1.7560316) q[0];
rz(-0.7631453) q[1];
sx q[1];
rz(-1.4475854) q[1];
sx q[1];
rz(0.10173434) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3107233) q[0];
sx q[0];
rz(-2.1889157) q[0];
sx q[0];
rz(-3.028426) q[0];
rz(-0.38282443) q[2];
sx q[2];
rz(-1.6144132) q[2];
sx q[2];
rz(2.5137465) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8489972) q[1];
sx q[1];
rz(-1.5096055) q[1];
sx q[1];
rz(-1.2996177) q[1];
rz(-2.2564628) q[3];
sx q[3];
rz(-0.99381522) q[3];
sx q[3];
rz(-0.5761522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.10852854) q[2];
sx q[2];
rz(-1.4427002) q[2];
sx q[2];
rz(-2.4326883) q[2];
rz(0.76062834) q[3];
sx q[3];
rz(-1.9899188) q[3];
sx q[3];
rz(2.2061548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075542299) q[0];
sx q[0];
rz(-1.9882555) q[0];
sx q[0];
rz(-2.4205038) q[0];
rz(2.1636294) q[1];
sx q[1];
rz(-2.5826192) q[1];
sx q[1];
rz(-1.5616547) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8472196) q[0];
sx q[0];
rz(-3.0356963) q[0];
sx q[0];
rz(-2.3787593) q[0];
rz(-pi) q[1];
rz(1.220854) q[2];
sx q[2];
rz(-2.5169543) q[2];
sx q[2];
rz(0.24764316) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.60266337) q[1];
sx q[1];
rz(-0.70405761) q[1];
sx q[1];
rz(-2.5119147) q[1];
x q[2];
rz(-2.7851289) q[3];
sx q[3];
rz(-0.93178669) q[3];
sx q[3];
rz(1.4872516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.28479031) q[2];
sx q[2];
rz(-0.52758104) q[2];
sx q[2];
rz(-1.52012) q[2];
rz(-0.4367477) q[3];
sx q[3];
rz(-0.89478651) q[3];
sx q[3];
rz(-2.6020218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0590416) q[0];
sx q[0];
rz(-2.3024547) q[0];
sx q[0];
rz(-2.3178597) q[0];
rz(1.1388904) q[1];
sx q[1];
rz(-1.3528119) q[1];
sx q[1];
rz(-1.1762071) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0498031) q[0];
sx q[0];
rz(-0.50423008) q[0];
sx q[0];
rz(-1.4780028) q[0];
rz(-2.5912639) q[2];
sx q[2];
rz(-0.68745733) q[2];
sx q[2];
rz(-0.94151173) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2580793) q[1];
sx q[1];
rz(-2.9032859) q[1];
sx q[1];
rz(-0.041447354) q[1];
x q[2];
rz(-1.5309991) q[3];
sx q[3];
rz(-2.1615513) q[3];
sx q[3];
rz(-1.9692957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1463683) q[2];
sx q[2];
rz(-0.58688846) q[2];
sx q[2];
rz(0.5874908) q[2];
rz(-0.38582173) q[3];
sx q[3];
rz(-2.2930175) q[3];
sx q[3];
rz(-0.90252701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0226444) q[0];
sx q[0];
rz(-0.97624874) q[0];
sx q[0];
rz(2.5318085) q[0];
rz(-1.3439641) q[1];
sx q[1];
rz(-1.1577497) q[1];
sx q[1];
rz(-2.9885898) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1399779) q[0];
sx q[0];
rz(-1.5904332) q[0];
sx q[0];
rz(0.46592142) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.581089) q[2];
sx q[2];
rz(-0.77559815) q[2];
sx q[2];
rz(-2.6960109) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.24125946) q[1];
sx q[1];
rz(-1.5521084) q[1];
sx q[1];
rz(0.92872031) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47346799) q[3];
sx q[3];
rz(-0.93515474) q[3];
sx q[3];
rz(-1.1421695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1028221) q[2];
sx q[2];
rz(-1.5615347) q[2];
sx q[2];
rz(0.62225303) q[2];
rz(-1.9410939) q[3];
sx q[3];
rz(-1.7222722) q[3];
sx q[3];
rz(-0.57435575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1822561) q[0];
sx q[0];
rz(-0.065956235) q[0];
sx q[0];
rz(2.7863853) q[0];
rz(1.1025053) q[1];
sx q[1];
rz(-3.1246154) q[1];
sx q[1];
rz(-0.65974081) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0416464) q[0];
sx q[0];
rz(-0.92586854) q[0];
sx q[0];
rz(1.3553331) q[0];
x q[1];
rz(1.2869455) q[2];
sx q[2];
rz(-2.7128993) q[2];
sx q[2];
rz(-2.5108166) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5641355) q[1];
sx q[1];
rz(-2.681985) q[1];
sx q[1];
rz(0.80776249) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0889782) q[3];
sx q[3];
rz(-0.89375118) q[3];
sx q[3];
rz(0.074124215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3416662) q[2];
sx q[2];
rz(-0.081192668) q[2];
sx q[2];
rz(2.9426835) q[2];
rz(-2.343446) q[3];
sx q[3];
rz(-1.1856368) q[3];
sx q[3];
rz(-1.2822436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2220919) q[0];
sx q[0];
rz(-1.6035447) q[0];
sx q[0];
rz(-1.0731687) q[0];
rz(-0.78492289) q[1];
sx q[1];
rz(-2.2833318) q[1];
sx q[1];
rz(-2.2699184) q[1];
rz(-1.4474519) q[2];
sx q[2];
rz(-1.6179177) q[2];
sx q[2];
rz(-0.23455591) q[2];
rz(-1.296923) q[3];
sx q[3];
rz(-0.92007888) q[3];
sx q[3];
rz(-2.32213) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
