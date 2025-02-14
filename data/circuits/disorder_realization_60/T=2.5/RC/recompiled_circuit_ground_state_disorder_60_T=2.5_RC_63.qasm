OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.085091703) q[0];
sx q[0];
rz(-0.29274517) q[0];
sx q[0];
rz(2.2226287) q[0];
rz(-0.38209823) q[1];
sx q[1];
rz(-2.9736019) q[1];
sx q[1];
rz(1.1408495) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1699693) q[0];
sx q[0];
rz(-0.21339082) q[0];
sx q[0];
rz(-0.9728699) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4130526) q[2];
sx q[2];
rz(-2.72914) q[2];
sx q[2];
rz(-2.0714687) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.39841043) q[1];
sx q[1];
rz(-2.7355621) q[1];
sx q[1];
rz(0.34161411) q[1];
rz(-pi) q[2];
x q[2];
rz(1.857418) q[3];
sx q[3];
rz(-2.686073) q[3];
sx q[3];
rz(-0.16715967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4484278) q[2];
sx q[2];
rz(-1.0591256) q[2];
sx q[2];
rz(1.9654174) q[2];
rz(-3.0991992) q[3];
sx q[3];
rz(-1.3375125) q[3];
sx q[3];
rz(0.5347518) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6723044) q[0];
sx q[0];
rz(-0.35457087) q[0];
sx q[0];
rz(2.9571423) q[0];
rz(-0.09672673) q[1];
sx q[1];
rz(-1.6177142) q[1];
sx q[1];
rz(-1.9116037) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7589129) q[0];
sx q[0];
rz(-1.935674) q[0];
sx q[0];
rz(1.7880102) q[0];
x q[1];
rz(2.582091) q[2];
sx q[2];
rz(-1.6116665) q[2];
sx q[2];
rz(0.80383435) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0945541) q[1];
sx q[1];
rz(-2.8287005) q[1];
sx q[1];
rz(-2.439586) q[1];
x q[2];
rz(-1.9490864) q[3];
sx q[3];
rz(-1.6069901) q[3];
sx q[3];
rz(0.32603273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35473287) q[2];
sx q[2];
rz(-0.41149461) q[2];
sx q[2];
rz(0.011431781) q[2];
rz(-1.8276021) q[3];
sx q[3];
rz(-1.7766137) q[3];
sx q[3];
rz(-2.438681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74362022) q[0];
sx q[0];
rz(-3.008606) q[0];
sx q[0];
rz(-0.61070329) q[0];
rz(0.01344219) q[1];
sx q[1];
rz(-1.2862658) q[1];
sx q[1];
rz(-0.59649831) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72815182) q[0];
sx q[0];
rz(-2.1755009) q[0];
sx q[0];
rz(2.938943) q[0];
rz(-pi) q[1];
rz(1.1223511) q[2];
sx q[2];
rz(-2.5150013) q[2];
sx q[2];
rz(0.5465051) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.92745078) q[1];
sx q[1];
rz(-0.79793707) q[1];
sx q[1];
rz(0.44254704) q[1];
x q[2];
rz(1.0168996) q[3];
sx q[3];
rz(-1.5427329) q[3];
sx q[3];
rz(-1.156776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4639123) q[2];
sx q[2];
rz(-1.8296506) q[2];
sx q[2];
rz(3.1413191) q[2];
rz(-1.6655946) q[3];
sx q[3];
rz(-0.6690343) q[3];
sx q[3];
rz(3.0345548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45604712) q[0];
sx q[0];
rz(-2.9750329) q[0];
sx q[0];
rz(1.1628994) q[0];
rz(-2.1641425) q[1];
sx q[1];
rz(-1.5403427) q[1];
sx q[1];
rz(-1.5553442) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.453131) q[0];
sx q[0];
rz(-1.5630088) q[0];
sx q[0];
rz(1.5611737) q[0];
rz(-pi) q[1];
rz(-0.78537099) q[2];
sx q[2];
rz(-1.5391239) q[2];
sx q[2];
rz(-2.9065064) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4443497) q[1];
sx q[1];
rz(-2.067579) q[1];
sx q[1];
rz(1.5543943) q[1];
rz(-2.8951696) q[3];
sx q[3];
rz(-2.0520254) q[3];
sx q[3];
rz(0.53477609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.9419452) q[2];
sx q[2];
rz(-0.5117828) q[2];
sx q[2];
rz(-0.23435782) q[2];
rz(2.6394081) q[3];
sx q[3];
rz(-2.1000704) q[3];
sx q[3];
rz(-0.65619367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8239044) q[0];
sx q[0];
rz(-2.3264139) q[0];
sx q[0];
rz(1.748388) q[0];
rz(2.4488917) q[1];
sx q[1];
rz(-2.0557978) q[1];
sx q[1];
rz(-2.0511973) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4596953) q[0];
sx q[0];
rz(-2.5715264) q[0];
sx q[0];
rz(1.1718114) q[0];
rz(2.1887921) q[2];
sx q[2];
rz(-1.1296461) q[2];
sx q[2];
rz(1.6447826) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0648374) q[1];
sx q[1];
rz(-1.4204111) q[1];
sx q[1];
rz(1.1527747) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1107858) q[3];
sx q[3];
rz(-1.705369) q[3];
sx q[3];
rz(2.7887044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.56110567) q[2];
sx q[2];
rz(-1.9876391) q[2];
sx q[2];
rz(2.6055276) q[2];
rz(2.8481893) q[3];
sx q[3];
rz(-1.9304201) q[3];
sx q[3];
rz(-0.48040473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90477657) q[0];
sx q[0];
rz(-2.1892956) q[0];
sx q[0];
rz(-3.0425518) q[0];
rz(2.1095443) q[1];
sx q[1];
rz(-2.2910304) q[1];
sx q[1];
rz(1.7031857) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8598762) q[0];
sx q[0];
rz(-1.3163573) q[0];
sx q[0];
rz(-2.2734725) q[0];
x q[1];
rz(1.4969669) q[2];
sx q[2];
rz(-1.2880688) q[2];
sx q[2];
rz(1.7746995) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9248396) q[1];
sx q[1];
rz(-2.3509563) q[1];
sx q[1];
rz(3.1086712) q[1];
rz(-pi) q[2];
rz(0.29797192) q[3];
sx q[3];
rz(-1.4557585) q[3];
sx q[3];
rz(-1.486766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9859887) q[2];
sx q[2];
rz(-2.6376548) q[2];
sx q[2];
rz(-0.8209374) q[2];
rz(1.4486897) q[3];
sx q[3];
rz(-0.71805787) q[3];
sx q[3];
rz(1.6528355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89707017) q[0];
sx q[0];
rz(-0.91461602) q[0];
sx q[0];
rz(0.50265092) q[0];
rz(-1.9082759) q[1];
sx q[1];
rz(-0.93524593) q[1];
sx q[1];
rz(2.1615692) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5222447) q[0];
sx q[0];
rz(-2.7947593) q[0];
sx q[0];
rz(2.3330487) q[0];
rz(-pi) q[1];
rz(1.2792311) q[2];
sx q[2];
rz(-0.81613651) q[2];
sx q[2];
rz(1.7486435) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7777694) q[1];
sx q[1];
rz(-2.2968493) q[1];
sx q[1];
rz(-1.2593624) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3408889) q[3];
sx q[3];
rz(-0.74384825) q[3];
sx q[3];
rz(0.4901674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.22214733) q[2];
sx q[2];
rz(-1.28747) q[2];
sx q[2];
rz(-1.3745098) q[2];
rz(1.3162656) q[3];
sx q[3];
rz(-2.2856183) q[3];
sx q[3];
rz(-2.159806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027243622) q[0];
sx q[0];
rz(-0.39334941) q[0];
sx q[0];
rz(0.91841665) q[0];
rz(-2.9367327) q[1];
sx q[1];
rz(-1.058895) q[1];
sx q[1];
rz(-1.4835666) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.525104) q[0];
sx q[0];
rz(-1.8601928) q[0];
sx q[0];
rz(2.9887236) q[0];
x q[1];
rz(-2.3307033) q[2];
sx q[2];
rz(-1.4810908) q[2];
sx q[2];
rz(0.16505884) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0369954) q[1];
sx q[1];
rz(-0.87338398) q[1];
sx q[1];
rz(-3.0265977) q[1];
rz(-pi) q[2];
rz(-2.4224255) q[3];
sx q[3];
rz(-1.6079818) q[3];
sx q[3];
rz(-2.3553762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1398937) q[2];
sx q[2];
rz(-2.6504982) q[2];
sx q[2];
rz(1.3341382) q[2];
rz(0.88207465) q[3];
sx q[3];
rz(-0.69552723) q[3];
sx q[3];
rz(-1.2923366) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0598711) q[0];
sx q[0];
rz(-0.43161714) q[0];
sx q[0];
rz(0.01817848) q[0];
rz(2.7320618) q[1];
sx q[1];
rz(-1.7117701) q[1];
sx q[1];
rz(-3.0407564) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5699974) q[0];
sx q[0];
rz(-2.2048108) q[0];
sx q[0];
rz(1.6136531) q[0];
rz(-pi) q[1];
rz(-1.6214431) q[2];
sx q[2];
rz(-1.1505923) q[2];
sx q[2];
rz(2.8977608) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3253357) q[1];
sx q[1];
rz(-0.44202572) q[1];
sx q[1];
rz(2.1620552) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0206817) q[3];
sx q[3];
rz(-1.3252581) q[3];
sx q[3];
rz(0.80805486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.33710256) q[2];
sx q[2];
rz(-3.0323995) q[2];
sx q[2];
rz(-1.8936554) q[2];
rz(1.9167871) q[3];
sx q[3];
rz(-2.2962511) q[3];
sx q[3];
rz(0.065751806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2357904) q[0];
sx q[0];
rz(-1.6658655) q[0];
sx q[0];
rz(0.32050785) q[0];
rz(-0.39101741) q[1];
sx q[1];
rz(-2.0000439) q[1];
sx q[1];
rz(1.549622) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81330106) q[0];
sx q[0];
rz(-1.3718954) q[0];
sx q[0];
rz(1.3764253) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7334601) q[2];
sx q[2];
rz(-2.0792336) q[2];
sx q[2];
rz(-0.66772738) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1913073) q[1];
sx q[1];
rz(-2.7627594) q[1];
sx q[1];
rz(-1.3326419) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5783089) q[3];
sx q[3];
rz(-1.9931317) q[3];
sx q[3];
rz(1.0653121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0512507) q[2];
sx q[2];
rz(-2.0917459) q[2];
sx q[2];
rz(2.7412097) q[2];
rz(-2.9054437) q[3];
sx q[3];
rz(-0.103424) q[3];
sx q[3];
rz(2.1307438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.878933) q[0];
sx q[0];
rz(-0.50983179) q[0];
sx q[0];
rz(-1.8806993) q[0];
rz(0.48514584) q[1];
sx q[1];
rz(-2.3427675) q[1];
sx q[1];
rz(2.6642703) q[1];
rz(1.0322124) q[2];
sx q[2];
rz(-0.24145024) q[2];
sx q[2];
rz(-0.35035124) q[2];
rz(-3.0780001) q[3];
sx q[3];
rz(-1.7983754) q[3];
sx q[3];
rz(1.6354431) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
