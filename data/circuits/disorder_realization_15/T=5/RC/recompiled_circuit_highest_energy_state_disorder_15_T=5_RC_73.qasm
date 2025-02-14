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
rz(-1.4544571) q[0];
sx q[0];
rz(-1.9782826) q[0];
sx q[0];
rz(-1.8736725) q[0];
rz(-1.4866225) q[1];
sx q[1];
rz(-1.2868737) q[1];
sx q[1];
rz(-0.16965228) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1592703) q[0];
sx q[0];
rz(-1.2125612) q[0];
sx q[0];
rz(1.0591749) q[0];
rz(-1.9965788) q[2];
sx q[2];
rz(-0.55934956) q[2];
sx q[2];
rz(-2.1970791) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4625068) q[1];
sx q[1];
rz(-0.39606491) q[1];
sx q[1];
rz(1.6258662) q[1];
rz(-pi) q[2];
rz(0.98244169) q[3];
sx q[3];
rz(-0.37068493) q[3];
sx q[3];
rz(1.6904168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.22592813) q[2];
sx q[2];
rz(-1.516284) q[2];
sx q[2];
rz(-0.020326745) q[2];
rz(0.23975553) q[3];
sx q[3];
rz(-0.25529796) q[3];
sx q[3];
rz(-2.3026626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0073485) q[0];
sx q[0];
rz(-2.5237995) q[0];
sx q[0];
rz(-2.0452621) q[0];
rz(-1.4758551) q[1];
sx q[1];
rz(-1.219039) q[1];
sx q[1];
rz(0.79442564) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3045469) q[0];
sx q[0];
rz(-1.908806) q[0];
sx q[0];
rz(-1.90453) q[0];
x q[1];
rz(0.6068318) q[2];
sx q[2];
rz(-0.85222679) q[2];
sx q[2];
rz(2.5963714) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.55358821) q[1];
sx q[1];
rz(-0.4376227) q[1];
sx q[1];
rz(-2.8904084) q[1];
x q[2];
rz(1.3737455) q[3];
sx q[3];
rz(-2.3757977) q[3];
sx q[3];
rz(1.1546419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.085792556) q[2];
sx q[2];
rz(-1.5689359) q[2];
sx q[2];
rz(0.22573486) q[2];
rz(-2.8977532) q[3];
sx q[3];
rz(-0.95832458) q[3];
sx q[3];
rz(-0.17914151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8826411) q[0];
sx q[0];
rz(-1.8283586) q[0];
sx q[0];
rz(-2.7146345) q[0];
rz(-0.34307617) q[1];
sx q[1];
rz(-1.1767358) q[1];
sx q[1];
rz(-0.25161904) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.204827) q[0];
sx q[0];
rz(-2.511095) q[0];
sx q[0];
rz(-0.38602324) q[0];
x q[1];
rz(-1.3400836) q[2];
sx q[2];
rz(-2.3614411) q[2];
sx q[2];
rz(-3.0456269) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.15089825) q[1];
sx q[1];
rz(-0.37726918) q[1];
sx q[1];
rz(-0.47028415) q[1];
x q[2];
rz(-3.0915458) q[3];
sx q[3];
rz(-2.280683) q[3];
sx q[3];
rz(-1.0753701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.130927) q[2];
sx q[2];
rz(-1.5074573) q[2];
sx q[2];
rz(2.6962386) q[2];
rz(1.2238067) q[3];
sx q[3];
rz(-1.8755951) q[3];
sx q[3];
rz(-2.7538917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-2.7447516) q[0];
sx q[0];
rz(-2.3362609) q[0];
sx q[0];
rz(-1.9824363) q[0];
rz(-1.9388439) q[1];
sx q[1];
rz(-1.1801327) q[1];
sx q[1];
rz(-0.73701) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4241959) q[0];
sx q[0];
rz(-1.6336461) q[0];
sx q[0];
rz(-0.24817384) q[0];
rz(-pi) q[1];
rz(-2.2463655) q[2];
sx q[2];
rz(-1.8631808) q[2];
sx q[2];
rz(-1.043373) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.23589489) q[1];
sx q[1];
rz(-1.3970084) q[1];
sx q[1];
rz(0.64069616) q[1];
x q[2];
rz(-1.6287732) q[3];
sx q[3];
rz(-1.7372469) q[3];
sx q[3];
rz(0.49515192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0143339) q[2];
sx q[2];
rz(-2.2517683) q[2];
sx q[2];
rz(0.8117525) q[2];
rz(-2.3938866) q[3];
sx q[3];
rz(-0.85993189) q[3];
sx q[3];
rz(-3.0809793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49374813) q[0];
sx q[0];
rz(-2.0525377) q[0];
sx q[0];
rz(1.7895948) q[0];
rz(-2.3494675) q[1];
sx q[1];
rz(-0.89901662) q[1];
sx q[1];
rz(2.1314714) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57841741) q[0];
sx q[0];
rz(-2.0978598) q[0];
sx q[0];
rz(2.7521247) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.13212684) q[2];
sx q[2];
rz(-0.39574049) q[2];
sx q[2];
rz(-0.49116116) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0466442) q[1];
sx q[1];
rz(-1.6502336) q[1];
sx q[1];
rz(0.55752505) q[1];
rz(-pi) q[2];
x q[2];
rz(1.406448) q[3];
sx q[3];
rz(-1.4557299) q[3];
sx q[3];
rz(-3.1361406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1048364) q[2];
sx q[2];
rz(-2.3827621) q[2];
sx q[2];
rz(1.3834312) q[2];
rz(-3.1116327) q[3];
sx q[3];
rz(-0.99830097) q[3];
sx q[3];
rz(1.8647319) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5112011) q[0];
sx q[0];
rz(-2.7250405) q[0];
sx q[0];
rz(3.0101486) q[0];
rz(2.1427872) q[1];
sx q[1];
rz(-1.1591594) q[1];
sx q[1];
rz(-1.3712032) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5509439) q[0];
sx q[0];
rz(-1.9134132) q[0];
sx q[0];
rz(0.37775535) q[0];
x q[1];
rz(0.68626257) q[2];
sx q[2];
rz(-1.4126083) q[2];
sx q[2];
rz(-0.37294086) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.35624337) q[1];
sx q[1];
rz(-1.8438854) q[1];
sx q[1];
rz(-2.9778984) q[1];
rz(-2.430702) q[3];
sx q[3];
rz(-1.5681453) q[3];
sx q[3];
rz(-2.7972935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3580631) q[2];
sx q[2];
rz(-0.72856599) q[2];
sx q[2];
rz(1.0895458) q[2];
rz(-2.774636) q[3];
sx q[3];
rz(-0.4042545) q[3];
sx q[3];
rz(0.75016108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.628196) q[0];
sx q[0];
rz(-2.3371526) q[0];
sx q[0];
rz(2.2688769) q[0];
rz(-1.5645507) q[1];
sx q[1];
rz(-1.6460452) q[1];
sx q[1];
rz(-2.4094792) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66595378) q[0];
sx q[0];
rz(-1.7717965) q[0];
sx q[0];
rz(0.71626407) q[0];
rz(-pi) q[1];
rz(-2.4560628) q[2];
sx q[2];
rz(-1.6674124) q[2];
sx q[2];
rz(-1.8308507) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8120576) q[1];
sx q[1];
rz(-1.1033711) q[1];
sx q[1];
rz(-2.9653564) q[1];
x q[2];
rz(0.096166178) q[3];
sx q[3];
rz(-1.3447666) q[3];
sx q[3];
rz(1.1975675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6096036) q[2];
sx q[2];
rz(-2.1405818) q[2];
sx q[2];
rz(-0.58094376) q[2];
rz(-2.0161435) q[3];
sx q[3];
rz(-1.1939129) q[3];
sx q[3];
rz(1.9862991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4693212) q[0];
sx q[0];
rz(-0.10675616) q[0];
sx q[0];
rz(-2.971055) q[0];
rz(1.8960309) q[1];
sx q[1];
rz(-1.5252557) q[1];
sx q[1];
rz(-2.4571498) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4218773) q[0];
sx q[0];
rz(-1.9027896) q[0];
sx q[0];
rz(-0.50194711) q[0];
x q[1];
rz(-0.23875321) q[2];
sx q[2];
rz(-0.49477067) q[2];
sx q[2];
rz(-0.35428167) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0741006) q[1];
sx q[1];
rz(-1.3252186) q[1];
sx q[1];
rz(2.8771993) q[1];
rz(-0.6489469) q[3];
sx q[3];
rz(-1.5228378) q[3];
sx q[3];
rz(-0.063604442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30274063) q[2];
sx q[2];
rz(-1.3417696) q[2];
sx q[2];
rz(1.1203241) q[2];
rz(-2.4773347) q[3];
sx q[3];
rz(-2.7393326) q[3];
sx q[3];
rz(-0.50814381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6054194) q[0];
sx q[0];
rz(-2.7689731) q[0];
sx q[0];
rz(2.5033409) q[0];
rz(-0.0087180184) q[1];
sx q[1];
rz(-2.0254841) q[1];
sx q[1];
rz(0.4062103) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8883977) q[0];
sx q[0];
rz(-2.8450845) q[0];
sx q[0];
rz(-1.1738846) q[0];
rz(-pi) q[1];
rz(-1.4582602) q[2];
sx q[2];
rz(-2.7055158) q[2];
sx q[2];
rz(-1.2051518) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8779341) q[1];
sx q[1];
rz(-2.3535427) q[1];
sx q[1];
rz(-1.7321943) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33274098) q[3];
sx q[3];
rz(-0.077692835) q[3];
sx q[3];
rz(1.2487703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0476039) q[2];
sx q[2];
rz(-1.1270019) q[2];
sx q[2];
rz(0.36063933) q[2];
rz(-0.7274729) q[3];
sx q[3];
rz(-2.45939) q[3];
sx q[3];
rz(0.31807652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8135081) q[0];
sx q[0];
rz(-0.94539517) q[0];
sx q[0];
rz(2.9720921) q[0];
rz(-2.4318579) q[1];
sx q[1];
rz(-1.4163481) q[1];
sx q[1];
rz(1.08606) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5616578) q[0];
sx q[0];
rz(-0.57889639) q[0];
sx q[0];
rz(-2.7282342) q[0];
rz(-pi) q[1];
rz(2.1745944) q[2];
sx q[2];
rz(-2.6259086) q[2];
sx q[2];
rz(2.5623851) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5063613) q[1];
sx q[1];
rz(-1.6690134) q[1];
sx q[1];
rz(1.0088831) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38137718) q[3];
sx q[3];
rz(-1.3998195) q[3];
sx q[3];
rz(1.387763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3544932) q[2];
sx q[2];
rz(-0.084986173) q[2];
sx q[2];
rz(0.3698012) q[2];
rz(-1.1981111) q[3];
sx q[3];
rz(-1.5912278) q[3];
sx q[3];
rz(2.2280391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-3.0030466) q[0];
sx q[0];
rz(-1.9771165) q[0];
sx q[0];
rz(1.8856915) q[0];
rz(2.1239602) q[1];
sx q[1];
rz(-1.7435278) q[1];
sx q[1];
rz(1.7494038) q[1];
rz(1.6412777) q[2];
sx q[2];
rz(-1.6432091) q[2];
sx q[2];
rz(2.5753449) q[2];
rz(0.16416488) q[3];
sx q[3];
rz(-1.4420684) q[3];
sx q[3];
rz(1.3185929) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
