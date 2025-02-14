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
rz(1.3996742) q[0];
sx q[0];
rz(-1.6784486) q[0];
sx q[0];
rz(1.680467) q[0];
rz(-0.48179102) q[1];
sx q[1];
rz(-0.81463373) q[1];
sx q[1];
rz(0.93946594) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50892539) q[0];
sx q[0];
rz(-1.5804176) q[0];
sx q[0];
rz(-1.6144572) q[0];
rz(-pi) q[1];
rz(2.6109735) q[2];
sx q[2];
rz(-1.24959) q[2];
sx q[2];
rz(1.2752334) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.8535244) q[1];
sx q[1];
rz(-1.5306092) q[1];
sx q[1];
rz(-0.010577135) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.788673) q[3];
sx q[3];
rz(-1.1718201) q[3];
sx q[3];
rz(2.6348595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7751094) q[2];
sx q[2];
rz(-1.6747687) q[2];
sx q[2];
rz(-0.62466204) q[2];
rz(0.96056739) q[3];
sx q[3];
rz(-1.4550236) q[3];
sx q[3];
rz(-0.87489405) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7333154) q[0];
sx q[0];
rz(-2.4508219) q[0];
sx q[0];
rz(-2.9498003) q[0];
rz(0.57418144) q[1];
sx q[1];
rz(-1.5230813) q[1];
sx q[1];
rz(0.40723732) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8005368) q[0];
sx q[0];
rz(-1.5828907) q[0];
sx q[0];
rz(-1.3950134) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3941462) q[2];
sx q[2];
rz(-1.6345858) q[2];
sx q[2];
rz(2.0690567) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.25000962) q[1];
sx q[1];
rz(-0.56743562) q[1];
sx q[1];
rz(-2.7944348) q[1];
rz(-pi) q[2];
rz(1.0683162) q[3];
sx q[3];
rz(-1.4921616) q[3];
sx q[3];
rz(-2.7689611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.189397) q[2];
sx q[2];
rz(-1.1789097) q[2];
sx q[2];
rz(2.4883545) q[2];
rz(-2.2630528) q[3];
sx q[3];
rz(-2.0432751) q[3];
sx q[3];
rz(3.0285335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99881309) q[0];
sx q[0];
rz(-1.4343364) q[0];
sx q[0];
rz(2.0913731) q[0];
rz(-2.5255919) q[1];
sx q[1];
rz(-1.417825) q[1];
sx q[1];
rz(0.91316694) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8693059) q[0];
sx q[0];
rz(-1.7401631) q[0];
sx q[0];
rz(-0.10793751) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8575322) q[2];
sx q[2];
rz(-1.456919) q[2];
sx q[2];
rz(-1.8658474) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.84308594) q[1];
sx q[1];
rz(-1.8131327) q[1];
sx q[1];
rz(1.0964811) q[1];
rz(2.2596879) q[3];
sx q[3];
rz(-1.8671452) q[3];
sx q[3];
rz(2.0810082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9203165) q[2];
sx q[2];
rz(-0.90711275) q[2];
sx q[2];
rz(3.082412) q[2];
rz(0.83909285) q[3];
sx q[3];
rz(-1.9234761) q[3];
sx q[3];
rz(2.2742417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7729823) q[0];
sx q[0];
rz(-1.0091857) q[0];
sx q[0];
rz(-0.34977812) q[0];
rz(1.3444258) q[1];
sx q[1];
rz(-2.5549922) q[1];
sx q[1];
rz(0.1127359) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7737425) q[0];
sx q[0];
rz(-0.58184701) q[0];
sx q[0];
rz(0.655158) q[0];
x q[1];
rz(2.8543831) q[2];
sx q[2];
rz(-2.2104077) q[2];
sx q[2];
rz(-2.2869598) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5391721) q[1];
sx q[1];
rz(-2.2090171) q[1];
sx q[1];
rz(-0.37810564) q[1];
x q[2];
rz(0.90972439) q[3];
sx q[3];
rz(-0.99791699) q[3];
sx q[3];
rz(-2.6343144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.747644) q[2];
sx q[2];
rz(-1.5851721) q[2];
sx q[2];
rz(-3.0843688) q[2];
rz(-1.8852437) q[3];
sx q[3];
rz(-2.5689954) q[3];
sx q[3];
rz(0.73067874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8771862) q[0];
sx q[0];
rz(-1.1671966) q[0];
sx q[0];
rz(-2.4315244) q[0];
rz(2.8801628) q[1];
sx q[1];
rz(-1.556004) q[1];
sx q[1];
rz(-2.1905621) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9316021) q[0];
sx q[0];
rz(-1.2438124) q[0];
sx q[0];
rz(2.6698584) q[0];
rz(-2.0567953) q[2];
sx q[2];
rz(-2.1446781) q[2];
sx q[2];
rz(1.6465555) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3465865) q[1];
sx q[1];
rz(-1.7730496) q[1];
sx q[1];
rz(2.4319885) q[1];
rz(-pi) q[2];
x q[2];
rz(3.042438) q[3];
sx q[3];
rz(-0.94010779) q[3];
sx q[3];
rz(-2.658228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6422673) q[2];
sx q[2];
rz(-0.32411164) q[2];
sx q[2];
rz(2.8322241) q[2];
rz(2.046509) q[3];
sx q[3];
rz(-1.1284004) q[3];
sx q[3];
rz(0.40422082) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22769044) q[0];
sx q[0];
rz(-2.9588283) q[0];
sx q[0];
rz(-0.90232724) q[0];
rz(-0.062571049) q[1];
sx q[1];
rz(-1.8726655) q[1];
sx q[1];
rz(-1.8960457) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6645836) q[0];
sx q[0];
rz(-2.4467108) q[0];
sx q[0];
rz(0.24212153) q[0];
x q[1];
rz(1.0311216) q[2];
sx q[2];
rz(-1.6189685) q[2];
sx q[2];
rz(2.0210514) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3712103) q[1];
sx q[1];
rz(-0.66714811) q[1];
sx q[1];
rz(-0.70828153) q[1];
rz(0.67288237) q[3];
sx q[3];
rz(-1.9599293) q[3];
sx q[3];
rz(1.0997888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7237599) q[2];
sx q[2];
rz(-1.2699026) q[2];
sx q[2];
rz(-1.343824) q[2];
rz(1.8580681) q[3];
sx q[3];
rz(-0.44461861) q[3];
sx q[3];
rz(-0.72492391) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1062687) q[0];
sx q[0];
rz(-0.037287354) q[0];
sx q[0];
rz(-0.37621793) q[0];
rz(-0.47209921) q[1];
sx q[1];
rz(-0.68196982) q[1];
sx q[1];
rz(2.6825486) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3694341) q[0];
sx q[0];
rz(-0.87352814) q[0];
sx q[0];
rz(0.56241547) q[0];
x q[1];
rz(1.7455581) q[2];
sx q[2];
rz(-1.1272001) q[2];
sx q[2];
rz(0.52973739) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7524779) q[1];
sx q[1];
rz(-1.6387617) q[1];
sx q[1];
rz(0.1631921) q[1];
x q[2];
rz(3.1230733) q[3];
sx q[3];
rz(-2.5438622) q[3];
sx q[3];
rz(1.2883119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.072319357) q[2];
sx q[2];
rz(-0.71162144) q[2];
sx q[2];
rz(-2.7834328) q[2];
rz(0.88851309) q[3];
sx q[3];
rz(-1.0212967) q[3];
sx q[3];
rz(2.2108938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0446123) q[0];
sx q[0];
rz(-0.29878578) q[0];
sx q[0];
rz(-0.97691798) q[0];
rz(2.3878429) q[1];
sx q[1];
rz(-1.4434554) q[1];
sx q[1];
rz(1.5026198) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96732765) q[0];
sx q[0];
rz(-2.2696724) q[0];
sx q[0];
rz(1.4682659) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66960488) q[2];
sx q[2];
rz(-1.9960224) q[2];
sx q[2];
rz(-2.2086672) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5411631) q[1];
sx q[1];
rz(-0.90879295) q[1];
sx q[1];
rz(-1.1837436) q[1];
rz(2.1945403) q[3];
sx q[3];
rz(-1.3823518) q[3];
sx q[3];
rz(2.541813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1961907) q[2];
sx q[2];
rz(-1.6795029) q[2];
sx q[2];
rz(0.84575829) q[2];
rz(-0.78019199) q[3];
sx q[3];
rz(-1.3581685) q[3];
sx q[3];
rz(0.63854027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82275259) q[0];
sx q[0];
rz(-1.0718811) q[0];
sx q[0];
rz(2.7986797) q[0];
rz(1.0036537) q[1];
sx q[1];
rz(-1.9099078) q[1];
sx q[1];
rz(0.71663219) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4591689) q[0];
sx q[0];
rz(-2.3083516) q[0];
sx q[0];
rz(-2.6618945) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8877257) q[2];
sx q[2];
rz(-1.0980596) q[2];
sx q[2];
rz(-0.46893083) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.20880675) q[1];
sx q[1];
rz(-0.98742551) q[1];
sx q[1];
rz(1.8933184) q[1];
rz(2.6541084) q[3];
sx q[3];
rz(-2.4085452) q[3];
sx q[3];
rz(1.2300075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0922682) q[2];
sx q[2];
rz(-1.6246139) q[2];
sx q[2];
rz(-2.1410904) q[2];
rz(2.3246824) q[3];
sx q[3];
rz(-1.2714081) q[3];
sx q[3];
rz(0.36417714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7526806) q[0];
sx q[0];
rz(-2.2081544) q[0];
sx q[0];
rz(-2.3626784) q[0];
rz(2.3498416) q[1];
sx q[1];
rz(-1.2261032) q[1];
sx q[1];
rz(2.7955999) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2591922) q[0];
sx q[0];
rz(-1.1244052) q[0];
sx q[0];
rz(1.6631699) q[0];
rz(-pi) q[1];
rz(1.2674763) q[2];
sx q[2];
rz(-1.5823213) q[2];
sx q[2];
rz(-1.8755902) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.83904356) q[1];
sx q[1];
rz(-1.5356774) q[1];
sx q[1];
rz(-2.1167064) q[1];
x q[2];
rz(1.168604) q[3];
sx q[3];
rz(-2.3762083) q[3];
sx q[3];
rz(-2.6266498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9927696) q[2];
sx q[2];
rz(-1.537348) q[2];
sx q[2];
rz(-1.0619304) q[2];
rz(2.3289833) q[3];
sx q[3];
rz(-1.142623) q[3];
sx q[3];
rz(-2.699471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2963895) q[0];
sx q[0];
rz(-2.3136105) q[0];
sx q[0];
rz(-1.4985341) q[0];
rz(-0.2801815) q[1];
sx q[1];
rz(-1.2138841) q[1];
sx q[1];
rz(-0.79634204) q[1];
rz(-2.7325999) q[2];
sx q[2];
rz(-1.5341285) q[2];
sx q[2];
rz(-1.6899012) q[2];
rz(-2.0186043) q[3];
sx q[3];
rz(-2.347694) q[3];
sx q[3];
rz(-2.6705042) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
