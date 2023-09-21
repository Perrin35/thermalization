OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7251627) q[0];
sx q[0];
rz(-3.0017612) q[0];
sx q[0];
rz(-0.60959417) q[0];
rz(0.66863376) q[1];
sx q[1];
rz(-2.2761087) q[1];
sx q[1];
rz(3.0545711) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4092769) q[0];
sx q[0];
rz(-1.9415932) q[0];
sx q[0];
rz(2.2229574) q[0];
x q[1];
rz(1.6992703) q[2];
sx q[2];
rz(-0.80768425) q[2];
sx q[2];
rz(1.3070004) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7056071) q[1];
sx q[1];
rz(-1.7090624) q[1];
sx q[1];
rz(2.0558753) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5243516) q[3];
sx q[3];
rz(-1.159045) q[3];
sx q[3];
rz(-0.1842894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0027851) q[2];
sx q[2];
rz(-1.3802718) q[2];
sx q[2];
rz(0.37386093) q[2];
rz(2.8047681) q[3];
sx q[3];
rz(-1.5954433) q[3];
sx q[3];
rz(-0.22836223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8841298) q[0];
sx q[0];
rz(-0.3733491) q[0];
sx q[0];
rz(1.9467547) q[0];
rz(-3.0589814) q[1];
sx q[1];
rz(-1.1673085) q[1];
sx q[1];
rz(-0.00037489051) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9449687) q[0];
sx q[0];
rz(-1.9163016) q[0];
sx q[0];
rz(-2.0161122) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20034321) q[2];
sx q[2];
rz(-1.8274954) q[2];
sx q[2];
rz(0.30470195) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2692711) q[1];
sx q[1];
rz(-1.3355796) q[1];
sx q[1];
rz(-2.5260731) q[1];
rz(3.0307816) q[3];
sx q[3];
rz(-2.6786945) q[3];
sx q[3];
rz(-0.09679951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.80883819) q[2];
sx q[2];
rz(-0.19503441) q[2];
sx q[2];
rz(-0.17671281) q[2];
rz(0.79408944) q[3];
sx q[3];
rz(-0.77787557) q[3];
sx q[3];
rz(-0.55364048) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2394543) q[0];
sx q[0];
rz(-2.4376526) q[0];
sx q[0];
rz(2.7368271) q[0];
rz(1.8602712) q[1];
sx q[1];
rz(-2.5679913) q[1];
sx q[1];
rz(-1.8331029) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3151911) q[0];
sx q[0];
rz(-2.0668525) q[0];
sx q[0];
rz(2.1000923) q[0];
rz(2.7576315) q[2];
sx q[2];
rz(-0.57120354) q[2];
sx q[2];
rz(2.2307894) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3209826) q[1];
sx q[1];
rz(-0.87249407) q[1];
sx q[1];
rz(-0.33336158) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43867302) q[3];
sx q[3];
rz(-1.6524501) q[3];
sx q[3];
rz(0.7338394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.219316) q[2];
sx q[2];
rz(-2.6185161) q[2];
sx q[2];
rz(3.1211839) q[2];
rz(1.071788) q[3];
sx q[3];
rz(-2.0139549) q[3];
sx q[3];
rz(0.66334692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.9999009) q[0];
sx q[0];
rz(-0.90943709) q[0];
sx q[0];
rz(-0.91598696) q[0];
rz(-2.6782716) q[1];
sx q[1];
rz(-2.0756192) q[1];
sx q[1];
rz(-1.0571009) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6592641) q[0];
sx q[0];
rz(-1.5963012) q[0];
sx q[0];
rz(-0.43735023) q[0];
rz(-pi) q[1];
rz(-0.96287231) q[2];
sx q[2];
rz(-2.0553556) q[2];
sx q[2];
rz(-1.0819266) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.108236) q[1];
sx q[1];
rz(-1.5844966) q[1];
sx q[1];
rz(1.800888) q[1];
rz(-pi) q[2];
rz(1.7980868) q[3];
sx q[3];
rz(-0.82633457) q[3];
sx q[3];
rz(-2.8518761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2477734) q[2];
sx q[2];
rz(-2.5017068) q[2];
sx q[2];
rz(2.7988953) q[2];
rz(-1.4438859) q[3];
sx q[3];
rz(-0.87564898) q[3];
sx q[3];
rz(-2.4263583) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7856359) q[0];
sx q[0];
rz(-2.5781093) q[0];
sx q[0];
rz(-0.086439565) q[0];
rz(1.7516288) q[1];
sx q[1];
rz(-2.2098863) q[1];
sx q[1];
rz(2.6729029) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2738004) q[0];
sx q[0];
rz(-2.3228354) q[0];
sx q[0];
rz(2.1208982) q[0];
rz(0.30829633) q[2];
sx q[2];
rz(-1.397965) q[2];
sx q[2];
rz(1.0002713) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3735818) q[1];
sx q[1];
rz(-1.7005159) q[1];
sx q[1];
rz(0.20357657) q[1];
rz(-pi) q[2];
rz(-0.080805578) q[3];
sx q[3];
rz(-2.3256133) q[3];
sx q[3];
rz(0.24085837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9995352) q[2];
sx q[2];
rz(-2.2613328) q[2];
sx q[2];
rz(2.2772677) q[2];
rz(0.21720973) q[3];
sx q[3];
rz(-0.61621284) q[3];
sx q[3];
rz(2.8488081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42751673) q[0];
sx q[0];
rz(-1.4070516) q[0];
sx q[0];
rz(-2.9445904) q[0];
rz(1.7794094) q[1];
sx q[1];
rz(-1.4809337) q[1];
sx q[1];
rz(-0.33624712) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67664458) q[0];
sx q[0];
rz(-2.0083545) q[0];
sx q[0];
rz(-0.69533555) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98724987) q[2];
sx q[2];
rz(-1.0580214) q[2];
sx q[2];
rz(0.057904569) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0834004) q[1];
sx q[1];
rz(-0.5914878) q[1];
sx q[1];
rz(-3.0565492) q[1];
rz(-pi) q[2];
rz(-1.398596) q[3];
sx q[3];
rz(-2.2599054) q[3];
sx q[3];
rz(0.25758753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6062935) q[2];
sx q[2];
rz(-1.4761816) q[2];
sx q[2];
rz(2.2952648) q[2];
rz(-1.9019295) q[3];
sx q[3];
rz(-1.2318434) q[3];
sx q[3];
rz(0.0011750778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3939312) q[0];
sx q[0];
rz(-1.0378391) q[0];
sx q[0];
rz(2.8826662) q[0];
rz(-1.3461643) q[1];
sx q[1];
rz(-1.3849473) q[1];
sx q[1];
rz(-1.0940201) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0897113) q[0];
sx q[0];
rz(-1.2298755) q[0];
sx q[0];
rz(-2.0345576) q[0];
x q[1];
rz(-1.6254243) q[2];
sx q[2];
rz(-2.9187181) q[2];
sx q[2];
rz(-2.5958027) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8158405) q[1];
sx q[1];
rz(-0.098512352) q[1];
sx q[1];
rz(-1.8580395) q[1];
x q[2];
rz(-3.0435211) q[3];
sx q[3];
rz(-1.5925928) q[3];
sx q[3];
rz(-1.9779713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.1743494) q[2];
sx q[2];
rz(-1.3568342) q[2];
sx q[2];
rz(-2.8209177) q[2];
rz(0.43618068) q[3];
sx q[3];
rz(-2.6685721) q[3];
sx q[3];
rz(-2.6935553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86673474) q[0];
sx q[0];
rz(-0.47645706) q[0];
sx q[0];
rz(1.0048237) q[0];
rz(0.59016219) q[1];
sx q[1];
rz(-2.2361123) q[1];
sx q[1];
rz(0.80387962) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.852254) q[0];
sx q[0];
rz(-0.58710557) q[0];
sx q[0];
rz(3.0703074) q[0];
rz(2.8477746) q[2];
sx q[2];
rz(-2.567798) q[2];
sx q[2];
rz(2.3125355) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1357321) q[1];
sx q[1];
rz(-1.5605968) q[1];
sx q[1];
rz(-1.6594396) q[1];
rz(-pi) q[2];
rz(0.11532468) q[3];
sx q[3];
rz(-0.66676312) q[3];
sx q[3];
rz(0.61570864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8089495) q[2];
sx q[2];
rz(-0.84471622) q[2];
sx q[2];
rz(2.1311029) q[2];
rz(0.60339749) q[3];
sx q[3];
rz(-1.6167275) q[3];
sx q[3];
rz(0.55019125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5378961) q[0];
sx q[0];
rz(-1.2768856) q[0];
sx q[0];
rz(2.7340775) q[0];
rz(-0.28911668) q[1];
sx q[1];
rz(-2.018785) q[1];
sx q[1];
rz(2.3908652) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84839612) q[0];
sx q[0];
rz(-1.8627394) q[0];
sx q[0];
rz(1.2743203) q[0];
x q[1];
rz(-2.8640792) q[2];
sx q[2];
rz(-0.29007426) q[2];
sx q[2];
rz(-1.8857423) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1322861) q[1];
sx q[1];
rz(-0.3711776) q[1];
sx q[1];
rz(-0.98307857) q[1];
rz(-pi) q[2];
rz(-1.5446072) q[3];
sx q[3];
rz(-1.2695754) q[3];
sx q[3];
rz(1.477004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0728545) q[2];
sx q[2];
rz(-0.72454238) q[2];
sx q[2];
rz(1.9753974) q[2];
rz(-1.7769622) q[3];
sx q[3];
rz(-2.3659673) q[3];
sx q[3];
rz(-0.021818074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5584548) q[0];
sx q[0];
rz(-2.3151509) q[0];
sx q[0];
rz(1.3903842) q[0];
rz(-2.8109) q[1];
sx q[1];
rz(-0.76534098) q[1];
sx q[1];
rz(-1.4601382) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6351785) q[0];
sx q[0];
rz(-1.9902475) q[0];
sx q[0];
rz(-1.0181396) q[0];
x q[1];
rz(-2.0748126) q[2];
sx q[2];
rz(-1.5504312) q[2];
sx q[2];
rz(0.45769826) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15522659) q[1];
sx q[1];
rz(-2.6990597) q[1];
sx q[1];
rz(2.7000484) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71698935) q[3];
sx q[3];
rz(-2.2319712) q[3];
sx q[3];
rz(-2.9644074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.34974393) q[2];
sx q[2];
rz(-1.8169553) q[2];
sx q[2];
rz(1.4617408) q[2];
rz(-2.0215624) q[3];
sx q[3];
rz(-0.62168613) q[3];
sx q[3];
rz(-2.6859443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5158952) q[0];
sx q[0];
rz(-1.6202171) q[0];
sx q[0];
rz(1.1465999) q[0];
rz(-1.760578) q[1];
sx q[1];
rz(-1.2881423) q[1];
sx q[1];
rz(-1.2013411) q[1];
rz(1.9258326) q[2];
sx q[2];
rz(-0.81510614) q[2];
sx q[2];
rz(0.19744273) q[2];
rz(-0.47386668) q[3];
sx q[3];
rz(-0.73054536) q[3];
sx q[3];
rz(1.2252145) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
