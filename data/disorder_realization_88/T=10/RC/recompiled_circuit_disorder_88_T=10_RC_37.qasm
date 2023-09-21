OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8945583) q[0];
sx q[0];
rz(-1.2556827) q[0];
sx q[0];
rz(2.8136301) q[0];
rz(-0.2344996) q[1];
sx q[1];
rz(3.3426715) q[1];
sx q[1];
rz(9.3333416) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5782195) q[0];
sx q[0];
rz(-1.1735998) q[0];
sx q[0];
rz(2.5114775) q[0];
x q[1];
rz(2.882471) q[2];
sx q[2];
rz(-1.7012351) q[2];
sx q[2];
rz(0.63333095) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4788471) q[1];
sx q[1];
rz(-1.6531684) q[1];
sx q[1];
rz(1.3326416) q[1];
rz(-pi) q[2];
rz(-0.9760194) q[3];
sx q[3];
rz(-1.3888437) q[3];
sx q[3];
rz(2.6973157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.66449195) q[2];
sx q[2];
rz(-0.97350073) q[2];
sx q[2];
rz(-2.0155902) q[2];
rz(-2.8664355) q[3];
sx q[3];
rz(-0.61029172) q[3];
sx q[3];
rz(0.90308213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1317516) q[0];
sx q[0];
rz(-2.4480902) q[0];
sx q[0];
rz(-2.4480208) q[0];
rz(-1.0961078) q[1];
sx q[1];
rz(-2.1577436) q[1];
sx q[1];
rz(0.19031659) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41102558) q[0];
sx q[0];
rz(-1.1682296) q[0];
sx q[0];
rz(2.9257724) q[0];
rz(-pi) q[1];
rz(-1.9637738) q[2];
sx q[2];
rz(-1.0699546) q[2];
sx q[2];
rz(2.9289392) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.39365921) q[1];
sx q[1];
rz(-0.84888443) q[1];
sx q[1];
rz(3.1118803) q[1];
rz(-2.6211561) q[3];
sx q[3];
rz(-2.3428829) q[3];
sx q[3];
rz(-0.48898104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6341614) q[2];
sx q[2];
rz(-2.4298411) q[2];
sx q[2];
rz(-1.2197536) q[2];
rz(-2.9988585) q[3];
sx q[3];
rz(-2.1295363) q[3];
sx q[3];
rz(-0.90863168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.398657) q[0];
sx q[0];
rz(-2.0516899) q[0];
sx q[0];
rz(-2.3213342) q[0];
rz(-0.29207686) q[1];
sx q[1];
rz(-2.0675817) q[1];
sx q[1];
rz(-1.8935727) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1093729) q[0];
sx q[0];
rz(-1.4782227) q[0];
sx q[0];
rz(-2.5412482) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9343611) q[2];
sx q[2];
rz(-1.6336294) q[2];
sx q[2];
rz(2.732423) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1420715) q[1];
sx q[1];
rz(-1.6377249) q[1];
sx q[1];
rz(2.4740069) q[1];
rz(-pi) q[2];
rz(-2.8415801) q[3];
sx q[3];
rz(-2.1999151) q[3];
sx q[3];
rz(-0.56738561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32039207) q[2];
sx q[2];
rz(-2.0911066) q[2];
sx q[2];
rz(1.2134264) q[2];
rz(-0.16472566) q[3];
sx q[3];
rz(-0.89403331) q[3];
sx q[3];
rz(0.081610672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40959013) q[0];
sx q[0];
rz(-1.859917) q[0];
sx q[0];
rz(2.9555292) q[0];
rz(-0.20446725) q[1];
sx q[1];
rz(-0.47195131) q[1];
sx q[1];
rz(-1.8444555) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1870677) q[0];
sx q[0];
rz(-1.6141119) q[0];
sx q[0];
rz(1.4403507) q[0];
x q[1];
rz(-1.859971) q[2];
sx q[2];
rz(-1.7316069) q[2];
sx q[2];
rz(2.5922054) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.12249494) q[1];
sx q[1];
rz(-2.2593453) q[1];
sx q[1];
rz(0.81329878) q[1];
x q[2];
rz(-1.7503731) q[3];
sx q[3];
rz(-1.5844126) q[3];
sx q[3];
rz(-2.7333958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.90594784) q[2];
sx q[2];
rz(-0.79259688) q[2];
sx q[2];
rz(-0.91919351) q[2];
rz(2.8202608) q[3];
sx q[3];
rz(-1.0682169) q[3];
sx q[3];
rz(1.2478158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17386757) q[0];
sx q[0];
rz(-0.50650948) q[0];
sx q[0];
rz(2.2633973) q[0];
rz(-1.8163266) q[1];
sx q[1];
rz(-1.7293431) q[1];
sx q[1];
rz(-1.426288) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.722055) q[0];
sx q[0];
rz(-0.66156045) q[0];
sx q[0];
rz(-1.7046335) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23767383) q[2];
sx q[2];
rz(-1.875068) q[2];
sx q[2];
rz(-1.8397457) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1713472) q[1];
sx q[1];
rz(-2.4023224) q[1];
sx q[1];
rz(1.7319748) q[1];
x q[2];
rz(-1.8910847) q[3];
sx q[3];
rz(-0.8257782) q[3];
sx q[3];
rz(-0.68009963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8695801) q[2];
sx q[2];
rz(-2*pi/13) q[2];
sx q[2];
rz(3.0997979) q[2];
rz(-3.0801008) q[3];
sx q[3];
rz(-1.0890591) q[3];
sx q[3];
rz(2.7048236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7866661) q[0];
sx q[0];
rz(-0.085952856) q[0];
sx q[0];
rz(-2.1110995) q[0];
rz(-2.4018535) q[1];
sx q[1];
rz(-1.61295) q[1];
sx q[1];
rz(-2.5700263) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8067236) q[0];
sx q[0];
rz(-2.3248701) q[0];
sx q[0];
rz(-0.73210324) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7320485) q[2];
sx q[2];
rz(-2.469392) q[2];
sx q[2];
rz(-1.6598998) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.188365) q[1];
sx q[1];
rz(-1.1168915) q[1];
sx q[1];
rz(0.26785775) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81482859) q[3];
sx q[3];
rz(-1.0199254) q[3];
sx q[3];
rz(1.8967472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.968154) q[2];
sx q[2];
rz(-2.3996694) q[2];
sx q[2];
rz(-0.56419939) q[2];
rz(0.12600222) q[3];
sx q[3];
rz(-1.6832451) q[3];
sx q[3];
rz(-0.29461598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0903704) q[0];
sx q[0];
rz(-0.19989656) q[0];
sx q[0];
rz(-0.71227658) q[0];
rz(2.6157216) q[1];
sx q[1];
rz(-2.7253175) q[1];
sx q[1];
rz(-2.4760822) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3367046) q[0];
sx q[0];
rz(-2.8993653) q[0];
sx q[0];
rz(-1.5664943) q[0];
rz(-pi) q[1];
rz(-2.9930816) q[2];
sx q[2];
rz(-2.2627137) q[2];
sx q[2];
rz(1.6910545) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9741386) q[1];
sx q[1];
rz(-1.1145089) q[1];
sx q[1];
rz(0.2445226) q[1];
x q[2];
rz(-1.8129559) q[3];
sx q[3];
rz(-1.4302963) q[3];
sx q[3];
rz(-1.5797918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9843288) q[2];
sx q[2];
rz(-2.4815346) q[2];
sx q[2];
rz(1.7187913) q[2];
rz(-0.11519365) q[3];
sx q[3];
rz(-2.6264103) q[3];
sx q[3];
rz(-2.9490024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.049906235) q[0];
sx q[0];
rz(-1.1356857) q[0];
sx q[0];
rz(-2.9507622) q[0];
rz(2.514839) q[1];
sx q[1];
rz(-1.0160867) q[1];
sx q[1];
rz(0.33871067) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3074293) q[0];
sx q[0];
rz(-1.2273664) q[0];
sx q[0];
rz(0.15983454) q[0];
x q[1];
rz(1.8024826) q[2];
sx q[2];
rz(-1.1294239) q[2];
sx q[2];
rz(-2.1053134) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.056811) q[1];
sx q[1];
rz(-1.6817131) q[1];
sx q[1];
rz(-2.0463498) q[1];
rz(-2.268928) q[3];
sx q[3];
rz(-2.7250395) q[3];
sx q[3];
rz(2.0779028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.64615858) q[2];
sx q[2];
rz(-0.63445264) q[2];
sx q[2];
rz(-0.70043606) q[2];
rz(2.2436079) q[3];
sx q[3];
rz(-1.286641) q[3];
sx q[3];
rz(-2.9680796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53196466) q[0];
sx q[0];
rz(-0.48159972) q[0];
sx q[0];
rz(-2.4832446) q[0];
rz(-2.530653) q[1];
sx q[1];
rz(-1.8860203) q[1];
sx q[1];
rz(-3.0019965) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0848207) q[0];
sx q[0];
rz(-3.0773101) q[0];
sx q[0];
rz(1.8875185) q[0];
rz(-pi) q[1];
rz(1.6330958) q[2];
sx q[2];
rz(-0.75714105) q[2];
sx q[2];
rz(1.3296668) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.86233625) q[1];
sx q[1];
rz(-1.5039872) q[1];
sx q[1];
rz(-0.52211296) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92630444) q[3];
sx q[3];
rz(-1.6127805) q[3];
sx q[3];
rz(0.55461649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5247941) q[2];
sx q[2];
rz(-1.0107661) q[2];
sx q[2];
rz(2.810478) q[2];
rz(2.3838499) q[3];
sx q[3];
rz(-2.7527633) q[3];
sx q[3];
rz(-0.087879114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2963294) q[0];
sx q[0];
rz(-2.1305278) q[0];
sx q[0];
rz(-2.9560126) q[0];
rz(-2.045385) q[1];
sx q[1];
rz(-2.9269693) q[1];
sx q[1];
rz(-1.4846444) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.164924) q[0];
sx q[0];
rz(-2.6161368) q[0];
sx q[0];
rz(1.0111965) q[0];
rz(-2.9160203) q[2];
sx q[2];
rz(-1.301287) q[2];
sx q[2];
rz(-3.0849506) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6757322) q[1];
sx q[1];
rz(-1.670174) q[1];
sx q[1];
rz(-1.9468716) q[1];
x q[2];
rz(-1.5649892) q[3];
sx q[3];
rz(-2.4366637) q[3];
sx q[3];
rz(0.73202902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.93402702) q[2];
sx q[2];
rz(-2.2521314) q[2];
sx q[2];
rz(0.55220848) q[2];
rz(0.77783716) q[3];
sx q[3];
rz(-0.85770291) q[3];
sx q[3];
rz(-0.51789969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26375297) q[0];
sx q[0];
rz(-1.629566) q[0];
sx q[0];
rz(-1.4019479) q[0];
rz(-2.9539625) q[1];
sx q[1];
rz(-1.7788806) q[1];
sx q[1];
rz(2.3685041) q[1];
rz(0.4734584) q[2];
sx q[2];
rz(-0.31242328) q[2];
sx q[2];
rz(-1.7996126) q[2];
rz(-0.7836624) q[3];
sx q[3];
rz(-0.85859921) q[3];
sx q[3];
rz(2.8006299) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];