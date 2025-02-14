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
rz(1.1433831) q[0];
sx q[0];
rz(-1.5902061) q[0];
sx q[0];
rz(-2.8144612) q[0];
rz(1.7786572) q[1];
sx q[1];
rz(-2.5442446) q[1];
sx q[1];
rz(2.8391431) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0344447) q[0];
sx q[0];
rz(-1.4054789) q[0];
sx q[0];
rz(-1.3470696) q[0];
rz(2.601702) q[2];
sx q[2];
rz(-0.98289012) q[2];
sx q[2];
rz(-1.8251652) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.55284269) q[1];
sx q[1];
rz(-1.8859914) q[1];
sx q[1];
rz(-2.9957459) q[1];
rz(-pi) q[2];
rz(0.20154736) q[3];
sx q[3];
rz(-2.6601861) q[3];
sx q[3];
rz(2.4498482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7678307) q[2];
sx q[2];
rz(-0.094002873) q[2];
sx q[2];
rz(0.2730228) q[2];
rz(-0.36429575) q[3];
sx q[3];
rz(-1.7184869) q[3];
sx q[3];
rz(0.020615904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
rz(-1.6881707) q[0];
sx q[0];
rz(-1.0978798) q[0];
sx q[0];
rz(1.706644) q[0];
rz(0.20225987) q[1];
sx q[1];
rz(-0.63019284) q[1];
sx q[1];
rz(-0.30466255) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7633491) q[0];
sx q[0];
rz(-1.8190666) q[0];
sx q[0];
rz(-0.85064485) q[0];
rz(-0.29918798) q[2];
sx q[2];
rz(-1.9721037) q[2];
sx q[2];
rz(-2.4517259) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4604291) q[1];
sx q[1];
rz(-0.39084706) q[1];
sx q[1];
rz(-2.6397698) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84159236) q[3];
sx q[3];
rz(-2.9805956) q[3];
sx q[3];
rz(0.58486933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.27663407) q[2];
sx q[2];
rz(-0.43143299) q[2];
sx q[2];
rz(1.6032093) q[2];
rz(0.81613427) q[3];
sx q[3];
rz(-0.82681257) q[3];
sx q[3];
rz(-0.89474595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5012539) q[0];
sx q[0];
rz(-2.8479939) q[0];
sx q[0];
rz(1.5119934) q[0];
rz(1.8005796) q[1];
sx q[1];
rz(-2.2684542) q[1];
sx q[1];
rz(0.27892932) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19353774) q[0];
sx q[0];
rz(-2.9729261) q[0];
sx q[0];
rz(2.4381766) q[0];
rz(-pi) q[1];
rz(0.52742676) q[2];
sx q[2];
rz(-2.1373539) q[2];
sx q[2];
rz(1.1184831) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2532369) q[1];
sx q[1];
rz(-2.2742892) q[1];
sx q[1];
rz(0.40526595) q[1];
x q[2];
rz(0.0035905152) q[3];
sx q[3];
rz(-0.85455214) q[3];
sx q[3];
rz(0.2936756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5035847) q[2];
sx q[2];
rz(-1.6562485) q[2];
sx q[2];
rz(0.33835641) q[2];
rz(1.407297) q[3];
sx q[3];
rz(-2.0140078) q[3];
sx q[3];
rz(-2.1710763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19891837) q[0];
sx q[0];
rz(-0.022138683) q[0];
sx q[0];
rz(2.4778147) q[0];
rz(-0.55555073) q[1];
sx q[1];
rz(-1.5063565) q[1];
sx q[1];
rz(1.7864071) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2289425) q[0];
sx q[0];
rz(-1.602409) q[0];
sx q[0];
rz(-2.4949882) q[0];
rz(-pi) q[1];
rz(1.2122447) q[2];
sx q[2];
rz(-1.5927763) q[2];
sx q[2];
rz(0.41273261) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.042204) q[1];
sx q[1];
rz(-0.70800938) q[1];
sx q[1];
rz(-0.96570496) q[1];
rz(-1.438645) q[3];
sx q[3];
rz(-2.6717917) q[3];
sx q[3];
rz(2.2061004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3163471) q[2];
sx q[2];
rz(-0.97858846) q[2];
sx q[2];
rz(0.52581954) q[2];
rz(0.60872269) q[3];
sx q[3];
rz(-0.63568297) q[3];
sx q[3];
rz(2.4354602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71259251) q[0];
sx q[0];
rz(-1.9569995) q[0];
sx q[0];
rz(-0.96479601) q[0];
rz(2.7305799) q[1];
sx q[1];
rz(-0.95714584) q[1];
sx q[1];
rz(-1.1119276) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72090688) q[0];
sx q[0];
rz(-1.5474651) q[0];
sx q[0];
rz(-0.015714808) q[0];
rz(-pi) q[1];
rz(-2.0093727) q[2];
sx q[2];
rz(-2.2242332) q[2];
sx q[2];
rz(-2.6613622) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7605684) q[1];
sx q[1];
rz(-2.3443017) q[1];
sx q[1];
rz(-2.3098195) q[1];
rz(1.7181362) q[3];
sx q[3];
rz(-1.8354355) q[3];
sx q[3];
rz(1.6291813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3045706) q[2];
sx q[2];
rz(-1.3304173) q[2];
sx q[2];
rz(-2.0988317) q[2];
rz(1.1834831) q[3];
sx q[3];
rz(-0.8684929) q[3];
sx q[3];
rz(-2.1709501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17827621) q[0];
sx q[0];
rz(-1.9177508) q[0];
sx q[0];
rz(0.29603145) q[0];
rz(-3.1405247) q[1];
sx q[1];
rz(-1.3384534) q[1];
sx q[1];
rz(-1.0329049) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4460239) q[0];
sx q[0];
rz(-0.006651314) q[0];
sx q[0];
rz(-2.5435996) q[0];
rz(-pi) q[1];
rz(-1.8708816) q[2];
sx q[2];
rz(-0.3870766) q[2];
sx q[2];
rz(-0.29190266) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.57622785) q[1];
sx q[1];
rz(-1.8180483) q[1];
sx q[1];
rz(-1.6281716) q[1];
rz(2.37866) q[3];
sx q[3];
rz(-1.8682533) q[3];
sx q[3];
rz(-0.92253387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7615776) q[2];
sx q[2];
rz(-1.6753316) q[2];
sx q[2];
rz(-1.638691) q[2];
rz(2.5401529) q[3];
sx q[3];
rz(-0.99965874) q[3];
sx q[3];
rz(-2.7433266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10694207) q[0];
sx q[0];
rz(-1.5926188) q[0];
sx q[0];
rz(-0.73367992) q[0];
rz(-0.94379464) q[1];
sx q[1];
rz(-1.5422041) q[1];
sx q[1];
rz(-2.5234047) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3933936) q[0];
sx q[0];
rz(-1.7500779) q[0];
sx q[0];
rz(-0.050672942) q[0];
x q[1];
rz(-0.84039495) q[2];
sx q[2];
rz(-1.5638509) q[2];
sx q[2];
rz(0.23887979) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71431464) q[1];
sx q[1];
rz(-0.66525092) q[1];
sx q[1];
rz(-1.5284431) q[1];
rz(-1.9312956) q[3];
sx q[3];
rz(-2.1749718) q[3];
sx q[3];
rz(1.9856356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9435297) q[2];
sx q[2];
rz(-1.1897831) q[2];
sx q[2];
rz(2.5243536) q[2];
rz(-1.1715568) q[3];
sx q[3];
rz(-2.7660683) q[3];
sx q[3];
rz(0.61081162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.010855762) q[0];
sx q[0];
rz(-0.78748381) q[0];
sx q[0];
rz(2.6686344) q[0];
rz(1.6821776) q[1];
sx q[1];
rz(-1.727203) q[1];
sx q[1];
rz(0.032729538) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9851042) q[0];
sx q[0];
rz(-1.5384983) q[0];
sx q[0];
rz(-1.558139) q[0];
rz(-pi) q[1];
rz(0.47686994) q[2];
sx q[2];
rz(-0.82399659) q[2];
sx q[2];
rz(-1.0007953) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.02272005) q[1];
sx q[1];
rz(-1.1274844) q[1];
sx q[1];
rz(-2.9037649) q[1];
rz(0.030278997) q[3];
sx q[3];
rz(-1.6784199) q[3];
sx q[3];
rz(0.054747907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98256436) q[2];
sx q[2];
rz(-2.6391567) q[2];
sx q[2];
rz(2.0459797) q[2];
rz(0.76121965) q[3];
sx q[3];
rz(-1.2844362) q[3];
sx q[3];
rz(0.42263862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6765167) q[0];
sx q[0];
rz(-0.07971555) q[0];
sx q[0];
rz(-2.4327143) q[0];
rz(-0.31157663) q[1];
sx q[1];
rz(-1.0919002) q[1];
sx q[1];
rz(-0.30837217) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73357433) q[0];
sx q[0];
rz(-2.7112506) q[0];
sx q[0];
rz(0.4056613) q[0];
rz(-pi) q[1];
rz(-1.9498655) q[2];
sx q[2];
rz(-1.1806025) q[2];
sx q[2];
rz(-1.6782111) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2110097) q[1];
sx q[1];
rz(-0.66291729) q[1];
sx q[1];
rz(0.40056132) q[1];
rz(-pi) q[2];
rz(-2.3544054) q[3];
sx q[3];
rz(-1.3861674) q[3];
sx q[3];
rz(0.61239374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.10820216) q[2];
sx q[2];
rz(-1.4806925) q[2];
sx q[2];
rz(-2.9208276) q[2];
rz(2.0681785) q[3];
sx q[3];
rz(-0.99075166) q[3];
sx q[3];
rz(-0.78607917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.49755001) q[0];
sx q[0];
rz(-0.86174091) q[0];
sx q[0];
rz(-2.1047237) q[0];
rz(1.4700302) q[1];
sx q[1];
rz(-2.122888) q[1];
sx q[1];
rz(-1.1933614) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65914175) q[0];
sx q[0];
rz(-1.5406784) q[0];
sx q[0];
rz(-0.44216604) q[0];
x q[1];
rz(0.5035053) q[2];
sx q[2];
rz(-0.76447884) q[2];
sx q[2];
rz(1.6703005) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2630849) q[1];
sx q[1];
rz(-0.89205884) q[1];
sx q[1];
rz(1.8128859) q[1];
x q[2];
rz(2.0641238) q[3];
sx q[3];
rz(-0.26255373) q[3];
sx q[3];
rz(1.5535823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4736453) q[2];
sx q[2];
rz(-1.8569943) q[2];
sx q[2];
rz(-0.13772193) q[2];
rz(-1.5406476) q[3];
sx q[3];
rz(-2.9100304) q[3];
sx q[3];
rz(0.93129492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.4731428) q[0];
sx q[0];
rz(-1.4746329) q[0];
sx q[0];
rz(2.0869577) q[0];
rz(2.6558381) q[1];
sx q[1];
rz(-1.2243441) q[1];
sx q[1];
rz(-1.4996554) q[1];
rz(2.457224) q[2];
sx q[2];
rz(-1.3791313) q[2];
sx q[2];
rz(-2.287938) q[2];
rz(-0.53268643) q[3];
sx q[3];
rz(-2.3045427) q[3];
sx q[3];
rz(2.7994313) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
