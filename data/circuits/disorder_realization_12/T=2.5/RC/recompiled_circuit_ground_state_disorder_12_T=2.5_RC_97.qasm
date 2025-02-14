OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.96149093) q[0];
sx q[0];
rz(-2.5749126) q[0];
sx q[0];
rz(0.79071796) q[0];
rz(2.2136731) q[1];
sx q[1];
rz(3.1766422) q[1];
sx q[1];
rz(9.1215134) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1643675) q[0];
sx q[0];
rz(-2.4365855) q[0];
sx q[0];
rz(2.7859119) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7838285) q[2];
sx q[2];
rz(-1.2105889) q[2];
sx q[2];
rz(-1.2586762) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78312633) q[1];
sx q[1];
rz(-0.74399656) q[1];
sx q[1];
rz(-2.3143155) q[1];
rz(-1.4879151) q[3];
sx q[3];
rz(-0.83717665) q[3];
sx q[3];
rz(-2.8745911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.37320331) q[2];
sx q[2];
rz(-0.36036569) q[2];
sx q[2];
rz(2.4599794) q[2];
rz(0.56053376) q[3];
sx q[3];
rz(-0.78442854) q[3];
sx q[3];
rz(-0.64422977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.949837) q[0];
sx q[0];
rz(-2.2332709) q[0];
sx q[0];
rz(0.35582304) q[0];
rz(-2.8882354) q[1];
sx q[1];
rz(-0.8496049) q[1];
sx q[1];
rz(-1.8960948) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0137164) q[0];
sx q[0];
rz(-0.53049131) q[0];
sx q[0];
rz(2.6878304) q[0];
rz(-2.4741089) q[2];
sx q[2];
rz(-2.3892623) q[2];
sx q[2];
rz(-2.9927727) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51735866) q[1];
sx q[1];
rz(-1.9962203) q[1];
sx q[1];
rz(-2.5555771) q[1];
rz(-2.3058071) q[3];
sx q[3];
rz(-2.4525149) q[3];
sx q[3];
rz(-0.101735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3972828) q[2];
sx q[2];
rz(-2.3645568) q[2];
sx q[2];
rz(-0.12766078) q[2];
rz(-2.4658261) q[3];
sx q[3];
rz(-2.6830169) q[3];
sx q[3];
rz(-0.42135409) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0449988) q[0];
sx q[0];
rz(-2.8969722) q[0];
sx q[0];
rz(-2.007572) q[0];
rz(2.2484312) q[1];
sx q[1];
rz(-0.90407073) q[1];
sx q[1];
rz(1.2718511) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9006827) q[0];
sx q[0];
rz(-0.56334746) q[0];
sx q[0];
rz(-1.5213837) q[0];
x q[1];
rz(1.261146) q[2];
sx q[2];
rz(-1.6463338) q[2];
sx q[2];
rz(-0.47020082) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.33128502) q[1];
sx q[1];
rz(-0.87281094) q[1];
sx q[1];
rz(2.8728102) q[1];
rz(-pi) q[2];
rz(-1.7904537) q[3];
sx q[3];
rz(-1.4079511) q[3];
sx q[3];
rz(1.4901249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9280055) q[2];
sx q[2];
rz(-2.2760133) q[2];
sx q[2];
rz(1.8013087) q[2];
rz(-1.2876997) q[3];
sx q[3];
rz(-1.6940593) q[3];
sx q[3];
rz(-2.652216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-0.61557788) q[0];
sx q[0];
rz(-0.57681334) q[0];
sx q[0];
rz(-2.4082129) q[0];
rz(2.3461657) q[1];
sx q[1];
rz(-2.9454102) q[1];
sx q[1];
rz(-1.9733852) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35287898) q[0];
sx q[0];
rz(-1.1601935) q[0];
sx q[0];
rz(-3.0263508) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57046095) q[2];
sx q[2];
rz(-2.1242363) q[2];
sx q[2];
rz(-1.9238071) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8655905) q[1];
sx q[1];
rz(-1.6011213) q[1];
sx q[1];
rz(-0.042976168) q[1];
rz(-pi) q[2];
rz(-0.73428728) q[3];
sx q[3];
rz(-0.81273116) q[3];
sx q[3];
rz(-1.6664291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9613793) q[2];
sx q[2];
rz(-0.81575477) q[2];
sx q[2];
rz(-2.1867895) q[2];
rz(2.09156) q[3];
sx q[3];
rz(-0.78767109) q[3];
sx q[3];
rz(-2.5900456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10057218) q[0];
sx q[0];
rz(-2.6028778) q[0];
sx q[0];
rz(2.4464497) q[0];
rz(-1.9255385) q[1];
sx q[1];
rz(-2.1247517) q[1];
sx q[1];
rz(0.020811828) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6489844) q[0];
sx q[0];
rz(-0.91726117) q[0];
sx q[0];
rz(-0.18006353) q[0];
rz(-pi) q[1];
rz(-2.8146652) q[2];
sx q[2];
rz(-0.86092868) q[2];
sx q[2];
rz(-0.29151379) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0328524) q[1];
sx q[1];
rz(-1.5592062) q[1];
sx q[1];
rz(-1.3316915) q[1];
x q[2];
rz(0.092904776) q[3];
sx q[3];
rz(-1.5094525) q[3];
sx q[3];
rz(0.35862229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1702599) q[2];
sx q[2];
rz(-0.33997619) q[2];
sx q[2];
rz(0.60410947) q[2];
rz(2.4938834) q[3];
sx q[3];
rz(-0.52102399) q[3];
sx q[3];
rz(-2.6807396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.1243509) q[0];
sx q[0];
rz(-1.0931953) q[0];
sx q[0];
rz(-2.7830615) q[0];
rz(-0.56309807) q[1];
sx q[1];
rz(-0.9022572) q[1];
sx q[1];
rz(-2.8360352) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76680032) q[0];
sx q[0];
rz(-1.2106441) q[0];
sx q[0];
rz(2.4712579) q[0];
x q[1];
rz(2.9592909) q[2];
sx q[2];
rz(-1.7392225) q[2];
sx q[2];
rz(-1.4252848) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8628775) q[1];
sx q[1];
rz(-0.76932615) q[1];
sx q[1];
rz(-2.3946986) q[1];
rz(-pi) q[2];
rz(0.73718585) q[3];
sx q[3];
rz(-3.0372826) q[3];
sx q[3];
rz(2.7953469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.87865889) q[2];
sx q[2];
rz(-0.96076751) q[2];
sx q[2];
rz(-0.052138694) q[2];
rz(-0.42929286) q[3];
sx q[3];
rz(-1.41058) q[3];
sx q[3];
rz(-0.2521635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6843863) q[0];
sx q[0];
rz(-2.6173213) q[0];
sx q[0];
rz(-2.9065409) q[0];
rz(2.5382407) q[1];
sx q[1];
rz(-2.3124606) q[1];
sx q[1];
rz(-0.57650173) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2430778) q[0];
sx q[0];
rz(-2.2295526) q[0];
sx q[0];
rz(0.9131477) q[0];
x q[1];
rz(-1.2320249) q[2];
sx q[2];
rz(-0.82451754) q[2];
sx q[2];
rz(0.19245779) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9707905) q[1];
sx q[1];
rz(-1.9681381) q[1];
sx q[1];
rz(0.34797619) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0888388) q[3];
sx q[3];
rz(-1.4432943) q[3];
sx q[3];
rz(-1.8021999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.06360402) q[2];
sx q[2];
rz(-2.3006907) q[2];
sx q[2];
rz(1.7951175) q[2];
rz(-0.50955647) q[3];
sx q[3];
rz(-1.2600803) q[3];
sx q[3];
rz(0.30663651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6154196) q[0];
sx q[0];
rz(-1.0831447) q[0];
sx q[0];
rz(2.7857067) q[0];
rz(0.7016167) q[1];
sx q[1];
rz(-2.9569148) q[1];
sx q[1];
rz(0.96431771) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50011694) q[0];
sx q[0];
rz(-1.7044412) q[0];
sx q[0];
rz(2.0290613) q[0];
rz(-pi) q[1];
rz(1.6379328) q[2];
sx q[2];
rz(-2.8481641) q[2];
sx q[2];
rz(0.98240438) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0810534) q[1];
sx q[1];
rz(-2.6463228) q[1];
sx q[1];
rz(-3.1002863) q[1];
rz(-pi) q[2];
rz(-1.6526493) q[3];
sx q[3];
rz(-1.2602031) q[3];
sx q[3];
rz(-1.6037343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4926766) q[2];
sx q[2];
rz(-0.70768386) q[2];
sx q[2];
rz(-2.3259582) q[2];
rz(-2.2260769) q[3];
sx q[3];
rz(-0.86796498) q[3];
sx q[3];
rz(-0.22667949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8409519) q[0];
sx q[0];
rz(-0.17423593) q[0];
sx q[0];
rz(-1.9006282) q[0];
rz(-0.10494431) q[1];
sx q[1];
rz(-0.34067708) q[1];
sx q[1];
rz(-0.45936432) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37840415) q[0];
sx q[0];
rz(-1.8109461) q[0];
sx q[0];
rz(1.2021007) q[0];
rz(-pi) q[1];
rz(0.26920374) q[2];
sx q[2];
rz(-1.891231) q[2];
sx q[2];
rz(-1.073369) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1314754) q[1];
sx q[1];
rz(-1.3595716) q[1];
sx q[1];
rz(-2.696871) q[1];
x q[2];
rz(1.7107443) q[3];
sx q[3];
rz(-1.0010747) q[3];
sx q[3];
rz(2.8960814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0505872) q[2];
sx q[2];
rz(-0.688474) q[2];
sx q[2];
rz(-3.0509994) q[2];
rz(2.374384) q[3];
sx q[3];
rz(-1.6588666) q[3];
sx q[3];
rz(1.4513133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7200658) q[0];
sx q[0];
rz(-2.2179715) q[0];
sx q[0];
rz(1.5371171) q[0];
rz(0.9556669) q[1];
sx q[1];
rz(-1.6803398) q[1];
sx q[1];
rz(-2.59424) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.132453) q[0];
sx q[0];
rz(-1.4910772) q[0];
sx q[0];
rz(3.135677) q[0];
rz(-2.6112767) q[2];
sx q[2];
rz(-1.2387453) q[2];
sx q[2];
rz(2.1830993) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0136614) q[1];
sx q[1];
rz(-1.3359002) q[1];
sx q[1];
rz(0.79551819) q[1];
x q[2];
rz(-0.71228446) q[3];
sx q[3];
rz(-1.9654719) q[3];
sx q[3];
rz(0.29721949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0305816) q[2];
sx q[2];
rz(-2.9647398) q[2];
sx q[2];
rz(-0.64876968) q[2];
rz(0.23160058) q[3];
sx q[3];
rz(-2.2571371) q[3];
sx q[3];
rz(0.086656682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7975174) q[0];
sx q[0];
rz(-1.9301013) q[0];
sx q[0];
rz(1.7420266) q[0];
rz(-2.4823785) q[1];
sx q[1];
rz(-1.2792239) q[1];
sx q[1];
rz(-1.4469133) q[1];
rz(3.0921583) q[2];
sx q[2];
rz(-1.7988322) q[2];
sx q[2];
rz(-3.1222406) q[2];
rz(-3.0869879) q[3];
sx q[3];
rz(-1.7163897) q[3];
sx q[3];
rz(-0.39672273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
