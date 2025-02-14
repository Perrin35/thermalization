OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.35535204) q[0];
sx q[0];
rz(2.8347637) q[0];
sx q[0];
rz(9.1260202) q[0];
rz(2.4755251) q[1];
sx q[1];
rz(-2.5112285) q[1];
sx q[1];
rz(-1.4730374) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0399892) q[0];
sx q[0];
rz(-0.77058661) q[0];
sx q[0];
rz(1.6165074) q[0];
rz(-1.5130784) q[2];
sx q[2];
rz(-0.30122631) q[2];
sx q[2];
rz(2.0026375) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8572644) q[1];
sx q[1];
rz(-0.92580399) q[1];
sx q[1];
rz(-2.2086992) q[1];
x q[2];
rz(-2.8657718) q[3];
sx q[3];
rz(-1.6383871) q[3];
sx q[3];
rz(-1.7662314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1851958) q[2];
sx q[2];
rz(-0.38072017) q[2];
sx q[2];
rz(1.8308651) q[2];
rz(-3.0500566) q[3];
sx q[3];
rz(-0.59713489) q[3];
sx q[3];
rz(2.6105647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.071534261) q[0];
sx q[0];
rz(-1.0907084) q[0];
sx q[0];
rz(-2.6614905) q[0];
rz(-1.9942888) q[1];
sx q[1];
rz(-0.6310178) q[1];
sx q[1];
rz(-0.70153418) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6525878) q[0];
sx q[0];
rz(-2.2092144) q[0];
sx q[0];
rz(-0.63853635) q[0];
rz(-pi) q[1];
rz(2.7240997) q[2];
sx q[2];
rz(-1.1402297) q[2];
sx q[2];
rz(2.9232962) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.32075927) q[1];
sx q[1];
rz(-2.1693128) q[1];
sx q[1];
rz(-2.5952336) q[1];
x q[2];
rz(-0.40579943) q[3];
sx q[3];
rz(-2.4630952) q[3];
sx q[3];
rz(0.12756995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.39701617) q[2];
sx q[2];
rz(-1.1822367) q[2];
sx q[2];
rz(1.5632632) q[2];
rz(-0.27966106) q[3];
sx q[3];
rz(-2.4249488) q[3];
sx q[3];
rz(-0.40951148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3438943) q[0];
sx q[0];
rz(-0.41223031) q[0];
sx q[0];
rz(-1.718234) q[0];
rz(-3.0176945) q[1];
sx q[1];
rz(-2.2975477) q[1];
sx q[1];
rz(1.5843676) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14649571) q[0];
sx q[0];
rz(-0.87281681) q[0];
sx q[0];
rz(-0.3905889) q[0];
x q[1];
rz(-0.84752797) q[2];
sx q[2];
rz(-1.28994) q[2];
sx q[2];
rz(2.6658863) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6137984) q[1];
sx q[1];
rz(-1.8168194) q[1];
sx q[1];
rz(-0.08155827) q[1];
rz(-pi) q[2];
rz(-2.4067114) q[3];
sx q[3];
rz(-0.85414825) q[3];
sx q[3];
rz(1.7295966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.065217) q[2];
sx q[2];
rz(-0.6535483) q[2];
sx q[2];
rz(0.93488133) q[2];
rz(-0.30385083) q[3];
sx q[3];
rz(-0.6128208) q[3];
sx q[3];
rz(0.8479619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76871753) q[0];
sx q[0];
rz(-2.3621552) q[0];
sx q[0];
rz(2.8745162) q[0];
rz(0.13485394) q[1];
sx q[1];
rz(-2.5649773) q[1];
sx q[1];
rz(-2.3042302) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7202111) q[0];
sx q[0];
rz(-2.1108339) q[0];
sx q[0];
rz(-2.8296986) q[0];
x q[1];
rz(0.58394152) q[2];
sx q[2];
rz(-2.2211233) q[2];
sx q[2];
rz(-3.1027628) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3738651) q[1];
sx q[1];
rz(-0.30555913) q[1];
sx q[1];
rz(-0.69626804) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5092974) q[3];
sx q[3];
rz(-1.617518) q[3];
sx q[3];
rz(-0.72544569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5460633) q[2];
sx q[2];
rz(-2.9661621) q[2];
sx q[2];
rz(-2.9626633) q[2];
rz(-1.9723802) q[3];
sx q[3];
rz(-1.7763276) q[3];
sx q[3];
rz(2.9235212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021407481) q[0];
sx q[0];
rz(-0.51161259) q[0];
sx q[0];
rz(2.2688493) q[0];
rz(-1.1574289) q[1];
sx q[1];
rz(-0.87481421) q[1];
sx q[1];
rz(0.016955888) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.884753) q[0];
sx q[0];
rz(-1.5829979) q[0];
sx q[0];
rz(-0.73464616) q[0];
rz(-pi) q[1];
rz(2.5418607) q[2];
sx q[2];
rz(-1.5625192) q[2];
sx q[2];
rz(0.042009609) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1106564) q[1];
sx q[1];
rz(-2.3157694) q[1];
sx q[1];
rz(2.7361672) q[1];
rz(-2.8722829) q[3];
sx q[3];
rz(-1.8519173) q[3];
sx q[3];
rz(-2.4914329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0852647) q[2];
sx q[2];
rz(-1.7076098) q[2];
sx q[2];
rz(2.8223574) q[2];
rz(2.4625835) q[3];
sx q[3];
rz(-1.7947936) q[3];
sx q[3];
rz(-0.80176789) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3646669) q[0];
sx q[0];
rz(-0.32061446) q[0];
sx q[0];
rz(-2.4989682) q[0];
rz(1.369426) q[1];
sx q[1];
rz(-0.6232999) q[1];
sx q[1];
rz(-0.97698897) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8606023) q[0];
sx q[0];
rz(-0.15781584) q[0];
sx q[0];
rz(0.67276038) q[0];
x q[1];
rz(-1.0426635) q[2];
sx q[2];
rz(-1.1234049) q[2];
sx q[2];
rz(-1.7698947) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9357537) q[1];
sx q[1];
rz(-0.69200695) q[1];
sx q[1];
rz(1.5924953) q[1];
rz(-pi) q[2];
x q[2];
rz(0.053835458) q[3];
sx q[3];
rz(-0.98568688) q[3];
sx q[3];
rz(1.1424292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.52648181) q[2];
sx q[2];
rz(-1.7340163) q[2];
sx q[2];
rz(1.9132445) q[2];
rz(-0.86722106) q[3];
sx q[3];
rz(-2.7286178) q[3];
sx q[3];
rz(2.373608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41384554) q[0];
sx q[0];
rz(-0.21738805) q[0];
sx q[0];
rz(-2.9687498) q[0];
rz(-2.3364283) q[1];
sx q[1];
rz(-2.5925437) q[1];
sx q[1];
rz(-3.0056312) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2826709) q[0];
sx q[0];
rz(-1.6902807) q[0];
sx q[0];
rz(-0.8651328) q[0];
rz(3.0375541) q[2];
sx q[2];
rz(-0.82471961) q[2];
sx q[2];
rz(-0.64268836) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0522642) q[1];
sx q[1];
rz(-1.58984) q[1];
sx q[1];
rz(-0.11472265) q[1];
rz(2.8118531) q[3];
sx q[3];
rz(-1.6981594) q[3];
sx q[3];
rz(1.691526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2148296) q[2];
sx q[2];
rz(-1.5952933) q[2];
sx q[2];
rz(2.9570441) q[2];
rz(0.18937011) q[3];
sx q[3];
rz(-2.6034077) q[3];
sx q[3];
rz(2.2034933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9823343) q[0];
sx q[0];
rz(-0.33894798) q[0];
sx q[0];
rz(-0.92639297) q[0];
rz(-3.1023846) q[1];
sx q[1];
rz(-0.67752939) q[1];
sx q[1];
rz(-1.0367941) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77039546) q[0];
sx q[0];
rz(-1.7465034) q[0];
sx q[0];
rz(-1.2716588) q[0];
rz(-pi) q[1];
rz(-0.37929566) q[2];
sx q[2];
rz(-0.46078983) q[2];
sx q[2];
rz(2.5658105) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7651789) q[1];
sx q[1];
rz(-2.661099) q[1];
sx q[1];
rz(-0.54554598) q[1];
rz(-1.758427) q[3];
sx q[3];
rz(-1.1867282) q[3];
sx q[3];
rz(1.7373067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.89256531) q[2];
sx q[2];
rz(-2.0165063) q[2];
sx q[2];
rz(-0.79891515) q[2];
rz(0.51698452) q[3];
sx q[3];
rz(-2.7345782) q[3];
sx q[3];
rz(-2.5911205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7901881) q[0];
sx q[0];
rz(-0.0049954448) q[0];
sx q[0];
rz(-1.5193526) q[0];
rz(-1.5937357) q[1];
sx q[1];
rz(-0.82031119) q[1];
sx q[1];
rz(0.69040745) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.522377) q[0];
sx q[0];
rz(-1.3478312) q[0];
sx q[0];
rz(-0.35943835) q[0];
rz(-pi) q[1];
rz(-1.4682653) q[2];
sx q[2];
rz(-2.3674115) q[2];
sx q[2];
rz(0.098086327) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4286297) q[1];
sx q[1];
rz(-0.74364122) q[1];
sx q[1];
rz(-0.78674591) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3314186) q[3];
sx q[3];
rz(-0.55174151) q[3];
sx q[3];
rz(0.61659471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4385628) q[2];
sx q[2];
rz(-2.5327693) q[2];
sx q[2];
rz(-2.1717066) q[2];
rz(-0.25740933) q[3];
sx q[3];
rz(-0.22203797) q[3];
sx q[3];
rz(1.359587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.53512204) q[0];
sx q[0];
rz(-0.73000014) q[0];
sx q[0];
rz(3.0560793) q[0];
rz(-0.27587786) q[1];
sx q[1];
rz(-0.62092263) q[1];
sx q[1];
rz(1.1891018) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10201926) q[0];
sx q[0];
rz(-1.9763401) q[0];
sx q[0];
rz(-1.9692375) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.027935426) q[2];
sx q[2];
rz(-0.31012529) q[2];
sx q[2];
rz(-1.5512229) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4170749) q[1];
sx q[1];
rz(-2.7226884) q[1];
sx q[1];
rz(-1.275283) q[1];
rz(-0.25190763) q[3];
sx q[3];
rz(-0.35874507) q[3];
sx q[3];
rz(-0.22049604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6659866) q[2];
sx q[2];
rz(-2.2624272) q[2];
sx q[2];
rz(0.48412588) q[2];
rz(0.13949805) q[3];
sx q[3];
rz(-0.77553427) q[3];
sx q[3];
rz(-0.16105306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6763247) q[0];
sx q[0];
rz(-1.9857255) q[0];
sx q[0];
rz(2.3406512) q[0];
rz(-1.6839266) q[1];
sx q[1];
rz(-1.6438345) q[1];
sx q[1];
rz(1.8297292) q[1];
rz(-0.64907907) q[2];
sx q[2];
rz(-1.4076283) q[2];
sx q[2];
rz(0.84244737) q[2];
rz(1.5994208) q[3];
sx q[3];
rz(-2.327649) q[3];
sx q[3];
rz(-1.5690371) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
