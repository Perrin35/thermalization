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
rz(3.0930003) q[0];
sx q[0];
rz(-1.4253923) q[0];
sx q[0];
rz(0.26560321) q[0];
rz(2.6939997) q[1];
sx q[1];
rz(-1.5040553) q[1];
sx q[1];
rz(0.12463364) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4701047) q[0];
sx q[0];
rz(-0.82297814) q[0];
sx q[0];
rz(0.95938553) q[0];
rz(-pi) q[1];
x q[1];
rz(0.12829576) q[2];
sx q[2];
rz(-1.1991457) q[2];
sx q[2];
rz(-2.9662463) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9198299) q[1];
sx q[1];
rz(-1.3897093) q[1];
sx q[1];
rz(0.20125514) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3105495) q[3];
sx q[3];
rz(-2.9141015) q[3];
sx q[3];
rz(-2.8854407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0894185) q[2];
sx q[2];
rz(-0.82145059) q[2];
sx q[2];
rz(-0.87665147) q[2];
rz(-0.18757251) q[3];
sx q[3];
rz(-0.98637527) q[3];
sx q[3];
rz(3.0228289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12980421) q[0];
sx q[0];
rz(-3.0846444) q[0];
sx q[0];
rz(2.7563128) q[0];
rz(-2.724559) q[1];
sx q[1];
rz(-0.72154355) q[1];
sx q[1];
rz(2.1293652) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6906805) q[0];
sx q[0];
rz(-0.60784303) q[0];
sx q[0];
rz(-1.9588821) q[0];
x q[1];
rz(-2.5778722) q[2];
sx q[2];
rz(-2.4032058) q[2];
sx q[2];
rz(-0.2283048) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.23823638) q[1];
sx q[1];
rz(-2.2975249) q[1];
sx q[1];
rz(1.958117) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1249097) q[3];
sx q[3];
rz(-0.7601217) q[3];
sx q[3];
rz(2.7393642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9932844) q[2];
sx q[2];
rz(-0.40996429) q[2];
sx q[2];
rz(-0.5245463) q[2];
rz(2.2927393) q[3];
sx q[3];
rz(-0.69791228) q[3];
sx q[3];
rz(-2.2822288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9808905) q[0];
sx q[0];
rz(-2.7301259) q[0];
sx q[0];
rz(2.8149862) q[0];
rz(1.7713361) q[1];
sx q[1];
rz(-2.0854918) q[1];
sx q[1];
rz(-0.036570963) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9865788) q[0];
sx q[0];
rz(-2.554092) q[0];
sx q[0];
rz(2.4632238) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0033135) q[2];
sx q[2];
rz(-0.9126513) q[2];
sx q[2];
rz(-2.7213241) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.72552437) q[1];
sx q[1];
rz(-1.6139133) q[1];
sx q[1];
rz(-1.7579702) q[1];
rz(-pi) q[2];
x q[2];
rz(2.960409) q[3];
sx q[3];
rz(-1.8474839) q[3];
sx q[3];
rz(1.7152648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6959491) q[2];
sx q[2];
rz(-2.8358938) q[2];
sx q[2];
rz(1.1841527) q[2];
rz(0.46237692) q[3];
sx q[3];
rz(-2.0591044) q[3];
sx q[3];
rz(2.6760127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9584123) q[0];
sx q[0];
rz(-1.6540225) q[0];
sx q[0];
rz(2.2271449) q[0];
rz(2.27521) q[1];
sx q[1];
rz(-1.8439081) q[1];
sx q[1];
rz(2.8241209) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7144517) q[0];
sx q[0];
rz(-1.1981816) q[0];
sx q[0];
rz(-1.6880549) q[0];
x q[1];
rz(-0.16902618) q[2];
sx q[2];
rz(-1.9155972) q[2];
sx q[2];
rz(1.930869) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.13617789) q[1];
sx q[1];
rz(-1.4168315) q[1];
sx q[1];
rz(1.8421768) q[1];
rz(-pi) q[2];
rz(1.1906641) q[3];
sx q[3];
rz(-1.5695342) q[3];
sx q[3];
rz(2.900687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0471961) q[2];
sx q[2];
rz(-0.4717584) q[2];
sx q[2];
rz(2.6867234) q[2];
rz(1.1826078) q[3];
sx q[3];
rz(-2.9136361) q[3];
sx q[3];
rz(-0.14127775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32893786) q[0];
sx q[0];
rz(-1.8998572) q[0];
sx q[0];
rz(-3.090233) q[0];
rz(-0.32737577) q[1];
sx q[1];
rz(-0.86762571) q[1];
sx q[1];
rz(3.0816269) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7236008) q[0];
sx q[0];
rz(-1.9034042) q[0];
sx q[0];
rz(-0.88653112) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9554638) q[2];
sx q[2];
rz(-0.95343381) q[2];
sx q[2];
rz(2.3996446) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9717945) q[1];
sx q[1];
rz(-2.4996539) q[1];
sx q[1];
rz(-2.1084014) q[1];
rz(-pi) q[2];
rz(2.9465605) q[3];
sx q[3];
rz(-2.2761506) q[3];
sx q[3];
rz(2.5759199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3943693) q[2];
sx q[2];
rz(-0.9129492) q[2];
sx q[2];
rz(-0.24820122) q[2];
rz(0.40211755) q[3];
sx q[3];
rz(-0.44200236) q[3];
sx q[3];
rz(-2.2615652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1345271) q[0];
sx q[0];
rz(-1.8777254) q[0];
sx q[0];
rz(1.7279351) q[0];
rz(-1.2029348) q[1];
sx q[1];
rz(-0.92034942) q[1];
sx q[1];
rz(-0.10341067) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72194874) q[0];
sx q[0];
rz(-1.6323315) q[0];
sx q[0];
rz(3.0716672) q[0];
rz(-pi) q[1];
rz(2.6194917) q[2];
sx q[2];
rz(-2.7828597) q[2];
sx q[2];
rz(-1.9821253) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.23650552) q[1];
sx q[1];
rz(-3.0154722) q[1];
sx q[1];
rz(2.990475) q[1];
x q[2];
rz(2.1391368) q[3];
sx q[3];
rz(-2.6939658) q[3];
sx q[3];
rz(0.8608495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.47205135) q[2];
sx q[2];
rz(-0.59110385) q[2];
sx q[2];
rz(2.9284787) q[2];
rz(-1.537568) q[3];
sx q[3];
rz(-3.0209916) q[3];
sx q[3];
rz(-0.11046256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9977426) q[0];
sx q[0];
rz(-0.85875964) q[0];
sx q[0];
rz(2.6570038) q[0];
rz(0.46802014) q[1];
sx q[1];
rz(-1.3222398) q[1];
sx q[1];
rz(-3.0048634) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3829684) q[0];
sx q[0];
rz(-2.6892085) q[0];
sx q[0];
rz(-2.1470619) q[0];
x q[1];
rz(-2.7523106) q[2];
sx q[2];
rz(-1.8912669) q[2];
sx q[2];
rz(-0.58347246) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9715106) q[1];
sx q[1];
rz(-0.85935837) q[1];
sx q[1];
rz(0.17632874) q[1];
x q[2];
rz(-2.8806512) q[3];
sx q[3];
rz(-2.709124) q[3];
sx q[3];
rz(-2.3245426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.7318657) q[2];
sx q[2];
rz(-0.12696433) q[2];
sx q[2];
rz(0.48180386) q[2];
rz(-1.9110154) q[3];
sx q[3];
rz(-1.9989719) q[3];
sx q[3];
rz(-3.0060911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-2.9442673) q[0];
sx q[0];
rz(-0.33536401) q[0];
sx q[0];
rz(2.770597) q[0];
rz(0.15759298) q[1];
sx q[1];
rz(-1.7540365) q[1];
sx q[1];
rz(0.88034672) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94778316) q[0];
sx q[0];
rz(-1.0968465) q[0];
sx q[0];
rz(3.0633846) q[0];
rz(-2.4331941) q[2];
sx q[2];
rz(-1.4428291) q[2];
sx q[2];
rz(-2.5789946) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.48260716) q[1];
sx q[1];
rz(-0.12276608) q[1];
sx q[1];
rz(1.6020847) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2796466) q[3];
sx q[3];
rz(-1.2978122) q[3];
sx q[3];
rz(-3.027944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2788435) q[2];
sx q[2];
rz(-0.66484386) q[2];
sx q[2];
rz(2.1319907) q[2];
rz(0.75262117) q[3];
sx q[3];
rz(-1.8372583) q[3];
sx q[3];
rz(2.2032004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9143592) q[0];
sx q[0];
rz(-0.14150134) q[0];
sx q[0];
rz(0.046382647) q[0];
rz(-1.4917689) q[1];
sx q[1];
rz(-1.7581419) q[1];
sx q[1];
rz(2.6683064) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21264874) q[0];
sx q[0];
rz(-1.4601652) q[0];
sx q[0];
rz(-1.380019) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32765179) q[2];
sx q[2];
rz(-1.1197512) q[2];
sx q[2];
rz(2.4137677) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9732194) q[1];
sx q[1];
rz(-2.1811317) q[1];
sx q[1];
rz(-1.8944505) q[1];
rz(-0.00086176894) q[3];
sx q[3];
rz(-1.6026845) q[3];
sx q[3];
rz(-2.3621967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16372323) q[2];
sx q[2];
rz(-2.551584) q[2];
sx q[2];
rz(-0.02656492) q[2];
rz(-0.50655347) q[3];
sx q[3];
rz(-0.81304628) q[3];
sx q[3];
rz(-2.626239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9764421) q[0];
sx q[0];
rz(-0.9557752) q[0];
sx q[0];
rz(-0.15560786) q[0];
rz(-0.43442976) q[1];
sx q[1];
rz(-2.5948718) q[1];
sx q[1];
rz(-0.2880407) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1032938) q[0];
sx q[0];
rz(-1.0464051) q[0];
sx q[0];
rz(2.45837) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1278438) q[2];
sx q[2];
rz(-1.7239769) q[2];
sx q[2];
rz(1.2438229) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.23032863) q[1];
sx q[1];
rz(-2.2161178) q[1];
sx q[1];
rz(-0.55379587) q[1];
rz(-pi) q[2];
rz(2.1197917) q[3];
sx q[3];
rz(-0.70216252) q[3];
sx q[3];
rz(-2.8150866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7274999) q[2];
sx q[2];
rz(-0.45656559) q[2];
sx q[2];
rz(-0.045819316) q[2];
rz(-3.0231061) q[3];
sx q[3];
rz(-0.89209569) q[3];
sx q[3];
rz(0.74046016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4988149) q[0];
sx q[0];
rz(-1.5121664) q[0];
sx q[0];
rz(1.2335516) q[0];
rz(-1.9998101) q[1];
sx q[1];
rz(-0.70527609) q[1];
sx q[1];
rz(-1.3377778) q[1];
rz(-2.4579688) q[2];
sx q[2];
rz(-2.0487006) q[2];
sx q[2];
rz(-0.13681199) q[2];
rz(-3.0790515) q[3];
sx q[3];
rz(-1.5086255) q[3];
sx q[3];
rz(3.0695741) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
