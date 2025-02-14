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
rz(0.64463717) q[0];
sx q[0];
rz(3.7312464) q[0];
sx q[0];
rz(10.51242) q[0];
rz(3.1261858) q[1];
sx q[1];
rz(-0.3897804) q[1];
sx q[1];
rz(-1.9773693) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.312266) q[0];
sx q[0];
rz(-0.13008936) q[0];
sx q[0];
rz(1.6072558) q[0];
rz(2.7373123) q[2];
sx q[2];
rz(-2.4089185) q[2];
sx q[2];
rz(2.3892185) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1401745) q[1];
sx q[1];
rz(-1.3571897) q[1];
sx q[1];
rz(2.1013538) q[1];
rz(1.4527997) q[3];
sx q[3];
rz(-1.495365) q[3];
sx q[3];
rz(2.3153967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2414134) q[2];
sx q[2];
rz(-2.5799077) q[2];
sx q[2];
rz(-2.4955595) q[2];
rz(-2.8863886) q[3];
sx q[3];
rz(-0.45707688) q[3];
sx q[3];
rz(-1.7469762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0274886) q[0];
sx q[0];
rz(-2.8059967) q[0];
sx q[0];
rz(-2.3345729) q[0];
rz(1.457816) q[1];
sx q[1];
rz(-2.5648263) q[1];
sx q[1];
rz(-1.797765) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4167359) q[0];
sx q[0];
rz(-1.8586419) q[0];
sx q[0];
rz(1.4132981) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.730663) q[2];
sx q[2];
rz(-1.0419453) q[2];
sx q[2];
rz(-1.5134606) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8763153) q[1];
sx q[1];
rz(-2.1723299) q[1];
sx q[1];
rz(2.5391891) q[1];
rz(-pi) q[2];
rz(-0.84193167) q[3];
sx q[3];
rz(-1.537383) q[3];
sx q[3];
rz(-0.15671003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.9606216) q[2];
sx q[2];
rz(-2.0714859) q[2];
sx q[2];
rz(-0.19372678) q[2];
rz(-1.6173897) q[3];
sx q[3];
rz(-0.75810713) q[3];
sx q[3];
rz(2.962842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5209565) q[0];
sx q[0];
rz(-0.23796029) q[0];
sx q[0];
rz(2.1654907) q[0];
rz(2.2528265) q[1];
sx q[1];
rz(-2.3975394) q[1];
sx q[1];
rz(2.9685453) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9219396) q[0];
sx q[0];
rz(-1.573607) q[0];
sx q[0];
rz(-1.5332743) q[0];
rz(0.52061527) q[2];
sx q[2];
rz(-1.6797425) q[2];
sx q[2];
rz(1.1586939) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0806607) q[1];
sx q[1];
rz(-2.073604) q[1];
sx q[1];
rz(-1.3814371) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8260396) q[3];
sx q[3];
rz(-2.8036661) q[3];
sx q[3];
rz(2.7546143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.139107) q[2];
sx q[2];
rz(-1.8128914) q[2];
sx q[2];
rz(-2.4277182) q[2];
rz(-0.21126963) q[3];
sx q[3];
rz(-2.1043089) q[3];
sx q[3];
rz(2.5762288) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7738889) q[0];
sx q[0];
rz(-1.1834894) q[0];
sx q[0];
rz(-2.0088038) q[0];
rz(0.47138131) q[1];
sx q[1];
rz(-1.8476723) q[1];
sx q[1];
rz(0.43101355) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16551183) q[0];
sx q[0];
rz(-1.3595194) q[0];
sx q[0];
rz(1.6863281) q[0];
x q[1];
rz(-0.36194275) q[2];
sx q[2];
rz(-1.2133994) q[2];
sx q[2];
rz(-1.5533642) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1872594) q[1];
sx q[1];
rz(-2.1012839) q[1];
sx q[1];
rz(1.9217291) q[1];
rz(-pi) q[2];
rz(-1.2005476) q[3];
sx q[3];
rz(-2.4863829) q[3];
sx q[3];
rz(-2.6225775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5203349) q[2];
sx q[2];
rz(-2.4020577) q[2];
sx q[2];
rz(2.7353103) q[2];
rz(1.9351561) q[3];
sx q[3];
rz(-1.0932357) q[3];
sx q[3];
rz(2.5549197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63221145) q[0];
sx q[0];
rz(-0.87566942) q[0];
sx q[0];
rz(-0.32633728) q[0];
rz(0.95411602) q[1];
sx q[1];
rz(-0.64847821) q[1];
sx q[1];
rz(-0.023524806) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3259567) q[0];
sx q[0];
rz(-1.1679497) q[0];
sx q[0];
rz(2.7340552) q[0];
x q[1];
rz(-0.40074375) q[2];
sx q[2];
rz(-0.6983499) q[2];
sx q[2];
rz(1.9633479) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6885029) q[1];
sx q[1];
rz(-1.2283235) q[1];
sx q[1];
rz(-0.30330412) q[1];
x q[2];
rz(2.708528) q[3];
sx q[3];
rz(-0.32363656) q[3];
sx q[3];
rz(1.3201081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.24713369) q[2];
sx q[2];
rz(-2.6829312) q[2];
sx q[2];
rz(-2.4776754) q[2];
rz(-2.2513921) q[3];
sx q[3];
rz(-1.5492487) q[3];
sx q[3];
rz(0.012705407) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7909872) q[0];
sx q[0];
rz(-2.7101639) q[0];
sx q[0];
rz(0.33082333) q[0];
rz(-0.36258969) q[1];
sx q[1];
rz(-1.9622842) q[1];
sx q[1];
rz(-0.51574743) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88721758) q[0];
sx q[0];
rz(-1.3332813) q[0];
sx q[0];
rz(-2.919048) q[0];
x q[1];
rz(0.2771122) q[2];
sx q[2];
rz(-0.61219575) q[2];
sx q[2];
rz(-0.0099358048) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1860604) q[1];
sx q[1];
rz(-1.493425) q[1];
sx q[1];
rz(0.87826985) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3456818) q[3];
sx q[3];
rz(-1.4778959) q[3];
sx q[3];
rz(1.9167629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.35046878) q[2];
sx q[2];
rz(-1.6914657) q[2];
sx q[2];
rz(-0.80751944) q[2];
rz(-0.098585419) q[3];
sx q[3];
rz(-2.2435296) q[3];
sx q[3];
rz(1.0928104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7035141) q[0];
sx q[0];
rz(-2.5683537) q[0];
sx q[0];
rz(-0.44952965) q[0];
rz(-2.4564157) q[1];
sx q[1];
rz(-0.40760577) q[1];
sx q[1];
rz(0.61946851) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9307053) q[0];
sx q[0];
rz(-0.31085098) q[0];
sx q[0];
rz(-3.0489281) q[0];
x q[1];
rz(0.52823587) q[2];
sx q[2];
rz(-1.6124469) q[2];
sx q[2];
rz(1.3170751) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5743235) q[1];
sx q[1];
rz(-2.3593214) q[1];
sx q[1];
rz(2.6560654) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18467091) q[3];
sx q[3];
rz(-1.7725367) q[3];
sx q[3];
rz(-1.7137485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5770136) q[2];
sx q[2];
rz(-1.6926293) q[2];
sx q[2];
rz(-2.9435834) q[2];
rz(-1.3624582) q[3];
sx q[3];
rz(-0.30006108) q[3];
sx q[3];
rz(-0.70039606) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9971767) q[0];
sx q[0];
rz(-0.62189019) q[0];
sx q[0];
rz(1.8560334) q[0];
rz(0.27664912) q[1];
sx q[1];
rz(-1.3686907) q[1];
sx q[1];
rz(3.0329774) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40688694) q[0];
sx q[0];
rz(-0.50841516) q[0];
sx q[0];
rz(-1.7395354) q[0];
rz(-pi) q[1];
rz(-0.94679657) q[2];
sx q[2];
rz(-1.5063926) q[2];
sx q[2];
rz(1.4381806) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.120879) q[1];
sx q[1];
rz(-2.2710137) q[1];
sx q[1];
rz(2.8597699) q[1];
x q[2];
rz(-1.9572419) q[3];
sx q[3];
rz(-1.4542011) q[3];
sx q[3];
rz(-2.1045096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.30839977) q[2];
sx q[2];
rz(-1.9449642) q[2];
sx q[2];
rz(-1.8673372) q[2];
rz(1.3537004) q[3];
sx q[3];
rz(-2.7022868) q[3];
sx q[3];
rz(-0.86658365) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15247791) q[0];
sx q[0];
rz(-2.415933) q[0];
sx q[0];
rz(-0.12055483) q[0];
rz(2.4001135) q[1];
sx q[1];
rz(-1.2040141) q[1];
sx q[1];
rz(3.0659058) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3593589) q[0];
sx q[0];
rz(-1.7348716) q[0];
sx q[0];
rz(0.30996451) q[0];
rz(-pi) q[1];
rz(1.4160749) q[2];
sx q[2];
rz(-0.46636367) q[2];
sx q[2];
rz(-1.5049962) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.2774128) q[1];
sx q[1];
rz(-2.2267739) q[1];
sx q[1];
rz(-1.5389062) q[1];
x q[2];
rz(-2.5406557) q[3];
sx q[3];
rz(-1.7225725) q[3];
sx q[3];
rz(-0.58702089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.072448298) q[2];
sx q[2];
rz(-2.2587903) q[2];
sx q[2];
rz(-0.34633386) q[2];
rz(-2.196178) q[3];
sx q[3];
rz(-1.8190705) q[3];
sx q[3];
rz(-2.1944428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.866975) q[0];
sx q[0];
rz(-1.9434384) q[0];
sx q[0];
rz(-0.44788885) q[0];
rz(-2.2121494) q[1];
sx q[1];
rz(-1.399469) q[1];
sx q[1];
rz(-2.4670752) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0513368) q[0];
sx q[0];
rz(-2.3469527) q[0];
sx q[0];
rz(-0.28213536) q[0];
x q[1];
rz(-1.0638164) q[2];
sx q[2];
rz(-0.21179767) q[2];
sx q[2];
rz(-1.3174881) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.50215534) q[1];
sx q[1];
rz(-1.5719255) q[1];
sx q[1];
rz(-1.7441895) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9895248) q[3];
sx q[3];
rz(-1.5018641) q[3];
sx q[3];
rz(1.513371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.21005361) q[2];
sx q[2];
rz(-1.9320107) q[2];
sx q[2];
rz(-2.8741969) q[2];
rz(-2.0343272) q[3];
sx q[3];
rz(-0.78247726) q[3];
sx q[3];
rz(-0.84899181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4103107) q[0];
sx q[0];
rz(-1.622643) q[0];
sx q[0];
rz(-1.0868764) q[0];
rz(0.80401737) q[1];
sx q[1];
rz(-1.6537279) q[1];
sx q[1];
rz(-2.0860685) q[1];
rz(1.7179391) q[2];
sx q[2];
rz(-1.4999119) q[2];
sx q[2];
rz(-2.9905408) q[2];
rz(0.08392423) q[3];
sx q[3];
rz(-0.90463432) q[3];
sx q[3];
rz(-2.6107364) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
