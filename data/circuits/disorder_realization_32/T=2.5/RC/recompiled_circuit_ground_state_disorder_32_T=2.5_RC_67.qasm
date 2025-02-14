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
rz(-1.8704432) q[0];
sx q[0];
rz(-0.76572642) q[0];
rz(5.8869047) q[1];
sx q[1];
rz(6.1904391) q[1];
sx q[1];
rz(9.9775597) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46229773) q[0];
sx q[0];
rz(-0.57749417) q[0];
sx q[0];
rz(2.3403843) q[0];
rz(-2.4713995) q[2];
sx q[2];
rz(-2.2171221) q[2];
sx q[2];
rz(1.0431511) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2491984) q[1];
sx q[1];
rz(-0.65649872) q[1];
sx q[1];
rz(1.0358443) q[1];
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
rz(-pi) q[1];
rz(1.7736241) q[2];
sx q[2];
rz(-1.6398733) q[2];
sx q[2];
rz(-0.088851301) q[2];
rz(-2.7446274) q[3];
sx q[3];
rz(-1.0395972) q[3];
sx q[3];
rz(1.6840434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48001584) q[0];
sx q[0];
rz(-0.50742298) q[0];
sx q[0];
rz(0.85451025) q[0];
rz(0.62141934) q[1];
sx q[1];
rz(-0.3365376) q[1];
sx q[1];
rz(-1.052676) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1922469) q[0];
sx q[0];
rz(-1.4345915) q[0];
sx q[0];
rz(-1.4679359) q[0];
rz(-pi) q[1];
rz(1.1243049) q[2];
sx q[2];
rz(-1.948775) q[2];
sx q[2];
rz(-1.2042827) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6748283) q[1];
sx q[1];
rz(-2.383814) q[1];
sx q[1];
rz(-1.2762875) q[1];
rz(1.5684116) q[3];
sx q[3];
rz(-2.5376476) q[3];
sx q[3];
rz(1.9238453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4738327) q[2];
sx q[2];
rz(-0.86060539) q[2];
sx q[2];
rz(0.94998002) q[2];
rz(2.1330323) q[3];
sx q[3];
rz(-2.5086094) q[3];
sx q[3];
rz(2.8620201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94473332) q[0];
sx q[0];
rz(-1.9254528) q[0];
sx q[0];
rz(-1.3470294) q[0];
rz(1.0579717) q[1];
sx q[1];
rz(-0.24669138) q[1];
sx q[1];
rz(-0.11014858) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6740883) q[0];
sx q[0];
rz(-1.2535998) q[0];
sx q[0];
rz(1.0870743) q[0];
rz(-pi) q[1];
rz(3.0579849) q[2];
sx q[2];
rz(-1.4232529) q[2];
sx q[2];
rz(-1.7718441) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.55787841) q[1];
sx q[1];
rz(-1.6772575) q[1];
sx q[1];
rz(2.540968) q[1];
rz(-pi) q[2];
rz(-2.0581117) q[3];
sx q[3];
rz(-1.9710473) q[3];
sx q[3];
rz(-2.6292173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.18016732) q[2];
sx q[2];
rz(-1.5632997) q[2];
sx q[2];
rz(3.0999198) q[2];
rz(-2.9336119) q[3];
sx q[3];
rz(-0.22170034) q[3];
sx q[3];
rz(2.8878133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5623077) q[0];
sx q[0];
rz(-1.3212181) q[0];
sx q[0];
rz(2.1606309) q[0];
rz(-0.43308577) q[1];
sx q[1];
rz(-1.4739477) q[1];
sx q[1];
rz(-1.513419) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63811231) q[0];
sx q[0];
rz(-0.44027087) q[0];
sx q[0];
rz(2.7087337) q[0];
x q[1];
rz(1.5879405) q[2];
sx q[2];
rz(-2.4701256) q[2];
sx q[2];
rz(1.4070828) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.6716587) q[1];
sx q[1];
rz(-2.0660559) q[1];
sx q[1];
rz(-3.0248066) q[1];
rz(-pi) q[2];
rz(2.8582005) q[3];
sx q[3];
rz(-2.1830354) q[3];
sx q[3];
rz(1.3536842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.607434) q[2];
sx q[2];
rz(-0.51924339) q[2];
sx q[2];
rz(2.329211) q[2];
rz(-1.8937998) q[3];
sx q[3];
rz(-1.8915853) q[3];
sx q[3];
rz(-1.8947424) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8060551) q[0];
sx q[0];
rz(-2.1195109) q[0];
sx q[0];
rz(1.442765) q[0];
rz(-0.45817786) q[1];
sx q[1];
rz(-1.6280326) q[1];
sx q[1];
rz(-0.75622574) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9423507) q[0];
sx q[0];
rz(-2.1926542) q[0];
sx q[0];
rz(-2.950475) q[0];
rz(0.67839038) q[2];
sx q[2];
rz(-2.3783461) q[2];
sx q[2];
rz(2.2445172) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.56179433) q[1];
sx q[1];
rz(-1.7057912) q[1];
sx q[1];
rz(1.4827824) q[1];
rz(-pi) q[2];
rz(-1.2268158) q[3];
sx q[3];
rz(-1.6046673) q[3];
sx q[3];
rz(0.86465166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.68236399) q[2];
sx q[2];
rz(-1.3882779) q[2];
sx q[2];
rz(2.1759822) q[2];
rz(-0.56973488) q[3];
sx q[3];
rz(-2.0163048) q[3];
sx q[3];
rz(1.4558571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4023034) q[0];
sx q[0];
rz(-1.9434513) q[0];
sx q[0];
rz(1.385561) q[0];
rz(-0.7631453) q[1];
sx q[1];
rz(-1.4475854) q[1];
sx q[1];
rz(0.10173434) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67416009) q[0];
sx q[0];
rz(-1.4786353) q[0];
sx q[0];
rz(-0.94964315) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11631233) q[2];
sx q[2];
rz(-2.7564133) q[2];
sx q[2];
rz(1.0508089) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.49470529) q[1];
sx q[1];
rz(-0.27783074) q[1];
sx q[1];
rz(1.3459413) q[1];
rz(-pi) q[2];
rz(2.4423994) q[3];
sx q[3];
rz(-1.0113888) q[3];
sx q[3];
rz(1.727211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.10852854) q[2];
sx q[2];
rz(-1.6988924) q[2];
sx q[2];
rz(-2.4326883) q[2];
rz(-0.76062834) q[3];
sx q[3];
rz(-1.1516738) q[3];
sx q[3];
rz(2.2061548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075542299) q[0];
sx q[0];
rz(-1.1533371) q[0];
sx q[0];
rz(0.72108889) q[0];
rz(0.9779633) q[1];
sx q[1];
rz(-0.55897346) q[1];
sx q[1];
rz(1.579938) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8472196) q[0];
sx q[0];
rz(-3.0356963) q[0];
sx q[0];
rz(0.76283331) q[0];
rz(-1.9207387) q[2];
sx q[2];
rz(-0.62463838) q[2];
sx q[2];
rz(2.8939495) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9813919) q[1];
sx q[1];
rz(-1.0202279) q[1];
sx q[1];
rz(-1.1070613) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35646371) q[3];
sx q[3];
rz(-0.93178669) q[3];
sx q[3];
rz(-1.4872516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.28479031) q[2];
sx q[2];
rz(-2.6140116) q[2];
sx q[2];
rz(1.52012) q[2];
rz(0.4367477) q[3];
sx q[3];
rz(-2.2468061) q[3];
sx q[3];
rz(0.53957087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0590416) q[0];
sx q[0];
rz(-2.3024547) q[0];
sx q[0];
rz(2.3178597) q[0];
rz(-1.1388904) q[1];
sx q[1];
rz(-1.7887807) q[1];
sx q[1];
rz(-1.1762071) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7438854) q[0];
sx q[0];
rz(-1.6155786) q[0];
sx q[0];
rz(-1.0683879) q[0];
rz(-pi) q[1];
rz(0.55032879) q[2];
sx q[2];
rz(-0.68745733) q[2];
sx q[2];
rz(-0.94151173) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9261647) q[1];
sx q[1];
rz(-1.3326982) q[1];
sx q[1];
rz(-1.5808616) q[1];
rz(-pi) q[2];
rz(-0.59112143) q[3];
sx q[3];
rz(-1.6038461) q[3];
sx q[3];
rz(2.765268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1463683) q[2];
sx q[2];
rz(-2.5547042) q[2];
sx q[2];
rz(-2.5541019) q[2];
rz(-0.38582173) q[3];
sx q[3];
rz(-0.84857517) q[3];
sx q[3];
rz(-2.2390656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0226444) q[0];
sx q[0];
rz(-0.97624874) q[0];
sx q[0];
rz(0.60978419) q[0];
rz(-1.7976286) q[1];
sx q[1];
rz(-1.1577497) q[1];
sx q[1];
rz(2.9885898) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5333819) q[0];
sx q[0];
rz(-2.675288) q[0];
sx q[0];
rz(-3.0979041) q[0];
x q[1];
rz(2.346368) q[2];
sx q[2];
rz(-1.5635901) q[2];
sx q[2];
rz(-2.0237271) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9003332) q[1];
sx q[1];
rz(-1.5894842) q[1];
sx q[1];
rz(2.2128723) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47346799) q[3];
sx q[3];
rz(-2.2064379) q[3];
sx q[3];
rz(1.9994232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1028221) q[2];
sx q[2];
rz(-1.5800579) q[2];
sx q[2];
rz(0.62225303) q[2];
rz(1.2004987) q[3];
sx q[3];
rz(-1.4193204) q[3];
sx q[3];
rz(-2.5672369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.95933652) q[0];
sx q[0];
rz(-0.065956235) q[0];
sx q[0];
rz(-0.35520735) q[0];
rz(2.0390873) q[1];
sx q[1];
rz(-3.1246154) q[1];
sx q[1];
rz(-2.4818518) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3400354) q[0];
sx q[0];
rz(-1.3990972) q[0];
sx q[0];
rz(-0.6562018) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1573359) q[2];
sx q[2];
rz(-1.454118) q[2];
sx q[2];
rz(1.199374) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5641355) q[1];
sx q[1];
rz(-2.681985) q[1];
sx q[1];
rz(0.80776249) q[1];
rz(-pi) q[2];
rz(-3.0889782) q[3];
sx q[3];
rz(-0.89375118) q[3];
sx q[3];
rz(3.0674684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79992646) q[2];
sx q[2];
rz(-3.0604) q[2];
sx q[2];
rz(-0.1989092) q[2];
rz(-0.79814664) q[3];
sx q[3];
rz(-1.9559559) q[3];
sx q[3];
rz(-1.2822436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2220919) q[0];
sx q[0];
rz(-1.6035447) q[0];
sx q[0];
rz(-1.0731687) q[0];
rz(-2.3566698) q[1];
sx q[1];
rz(-0.85826086) q[1];
sx q[1];
rz(0.8716743) q[1];
rz(-3.0941111) q[2];
sx q[2];
rz(-1.4475895) q[2];
sx q[2];
rz(-1.7995126) q[2];
rz(-1.8446696) q[3];
sx q[3];
rz(-2.2215138) q[3];
sx q[3];
rz(0.81946269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
