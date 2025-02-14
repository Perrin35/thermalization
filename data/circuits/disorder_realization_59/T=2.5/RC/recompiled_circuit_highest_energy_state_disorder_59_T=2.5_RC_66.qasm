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
rz(2.4177999) q[0];
sx q[0];
rz(-1.1136709) q[0];
sx q[0];
rz(2.1768575) q[0];
rz(-0.073702987) q[1];
sx q[1];
rz(-0.6414203) q[1];
sx q[1];
rz(-1.4944271) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7429356) q[0];
sx q[0];
rz(-1.9558859) q[0];
sx q[0];
rz(-0.72011098) q[0];
rz(-1.2624717) q[2];
sx q[2];
rz(-2.8092896) q[2];
sx q[2];
rz(-2.6118956) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.19077483) q[1];
sx q[1];
rz(-2.4147646) q[1];
sx q[1];
rz(0.72587691) q[1];
rz(0.046207059) q[3];
sx q[3];
rz(-2.0997092) q[3];
sx q[3];
rz(-2.5387704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73704314) q[2];
sx q[2];
rz(-1.5207542) q[2];
sx q[2];
rz(-2.9250195) q[2];
rz(2.6259322) q[3];
sx q[3];
rz(-1.0786846) q[3];
sx q[3];
rz(-0.083958538) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11230042) q[0];
sx q[0];
rz(-3.0700505) q[0];
sx q[0];
rz(-1.6610425) q[0];
rz(0.59313613) q[1];
sx q[1];
rz(-2.1430404) q[1];
sx q[1];
rz(-0.28233972) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7337228) q[0];
sx q[0];
rz(-1.5998915) q[0];
sx q[0];
rz(-0.12761527) q[0];
rz(-pi) q[1];
rz(1.3725874) q[2];
sx q[2];
rz(-1.4396126) q[2];
sx q[2];
rz(-3.1114357) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2609387) q[1];
sx q[1];
rz(-1.0323413) q[1];
sx q[1];
rz(1.7314584) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7670338) q[3];
sx q[3];
rz(-1.1653596) q[3];
sx q[3];
rz(-2.4949511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4234408) q[2];
sx q[2];
rz(-2.580692) q[2];
sx q[2];
rz(2.1924428) q[2];
rz(-0.078350457) q[3];
sx q[3];
rz(-1.4698942) q[3];
sx q[3];
rz(1.2300434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80701989) q[0];
sx q[0];
rz(-3.0834575) q[0];
sx q[0];
rz(-0.19244254) q[0];
rz(-0.9737393) q[1];
sx q[1];
rz(-2.1111919) q[1];
sx q[1];
rz(1.8691501) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6405435) q[0];
sx q[0];
rz(-1.8132117) q[0];
sx q[0];
rz(2.4512176) q[0];
rz(-2.9203833) q[2];
sx q[2];
rz(-0.84748805) q[2];
sx q[2];
rz(-0.71585594) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.13779894) q[1];
sx q[1];
rz(-1.8013211) q[1];
sx q[1];
rz(0.25584883) q[1];
rz(-1.1149241) q[3];
sx q[3];
rz(-1.6041363) q[3];
sx q[3];
rz(-3.03862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5729313) q[2];
sx q[2];
rz(-0.91721407) q[2];
sx q[2];
rz(0.60079637) q[2];
rz(2.3066547) q[3];
sx q[3];
rz(-2.2227414) q[3];
sx q[3];
rz(0.44343534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.7610953) q[0];
sx q[0];
rz(-1.1325855) q[0];
sx q[0];
rz(1.4196716) q[0];
rz(-2.8555866) q[1];
sx q[1];
rz(-1.1347457) q[1];
sx q[1];
rz(0.51925117) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4409281) q[0];
sx q[0];
rz(-1.6608149) q[0];
sx q[0];
rz(0.16914455) q[0];
rz(-pi) q[1];
rz(-1.6748739) q[2];
sx q[2];
rz(-1.7517183) q[2];
sx q[2];
rz(-1.6044001) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1545461) q[1];
sx q[1];
rz(-2.0813165) q[1];
sx q[1];
rz(-1.6493919) q[1];
rz(-0.052922225) q[3];
sx q[3];
rz(-2.5949946) q[3];
sx q[3];
rz(1.4846168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.1661735) q[2];
sx q[2];
rz(-0.68048802) q[2];
sx q[2];
rz(2.5892881) q[2];
rz(2.35899) q[3];
sx q[3];
rz(-1.3026594) q[3];
sx q[3];
rz(-0.18463126) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2833589) q[0];
sx q[0];
rz(-2.5193546) q[0];
sx q[0];
rz(-1.2828113) q[0];
rz(-1.917631) q[1];
sx q[1];
rz(-1.6061648) q[1];
sx q[1];
rz(-0.70972365) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2252475) q[0];
sx q[0];
rz(-3.1053251) q[0];
sx q[0];
rz(2.4929918) q[0];
x q[1];
rz(0.46825306) q[2];
sx q[2];
rz(-1.5893468) q[2];
sx q[2];
rz(-1.5520688) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.97204527) q[1];
sx q[1];
rz(-1.2575713) q[1];
sx q[1];
rz(0.07446988) q[1];
rz(-pi) q[2];
rz(2.9622795) q[3];
sx q[3];
rz(-1.9783522) q[3];
sx q[3];
rz(-0.61791622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84511406) q[2];
sx q[2];
rz(-2.5402386) q[2];
sx q[2];
rz(2.5884886) q[2];
rz(-0.84135711) q[3];
sx q[3];
rz(-2.0935121) q[3];
sx q[3];
rz(1.1355737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76310054) q[0];
sx q[0];
rz(-1.0569514) q[0];
sx q[0];
rz(-0.77504778) q[0];
rz(1.4729602) q[1];
sx q[1];
rz(-0.62455636) q[1];
sx q[1];
rz(0.83736173) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96356364) q[0];
sx q[0];
rz(-0.63321165) q[0];
sx q[0];
rz(-1.4073611) q[0];
rz(2.6858575) q[2];
sx q[2];
rz(-1.7349458) q[2];
sx q[2];
rz(-1.8388621) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0114767) q[1];
sx q[1];
rz(-1.1728047) q[1];
sx q[1];
rz(2.6813862) q[1];
rz(1.9894137) q[3];
sx q[3];
rz(-1.3411203) q[3];
sx q[3];
rz(-0.21462378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.9019292) q[2];
sx q[2];
rz(-1.9749125) q[2];
sx q[2];
rz(0.76266328) q[2];
rz(2.934382) q[3];
sx q[3];
rz(-1.8800294) q[3];
sx q[3];
rz(-2.5337849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36829456) q[0];
sx q[0];
rz(-2.2070856) q[0];
sx q[0];
rz(-1.6873129) q[0];
rz(-1.8621209) q[1];
sx q[1];
rz(-2.2984633) q[1];
sx q[1];
rz(-1.4131193) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93122497) q[0];
sx q[0];
rz(-2.702091) q[0];
sx q[0];
rz(2.6511433) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57369487) q[2];
sx q[2];
rz(-0.92373935) q[2];
sx q[2];
rz(1.1653125) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3154732) q[1];
sx q[1];
rz(-0.54448381) q[1];
sx q[1];
rz(-0.61750169) q[1];
rz(-pi) q[2];
rz(0.90274324) q[3];
sx q[3];
rz(-0.53865005) q[3];
sx q[3];
rz(-0.99942943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9965245) q[2];
sx q[2];
rz(-2.2729496) q[2];
sx q[2];
rz(-2.9800912) q[2];
rz(2.4072371) q[3];
sx q[3];
rz(-2.2752094) q[3];
sx q[3];
rz(1.3041147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9228009) q[0];
sx q[0];
rz(-0.53023338) q[0];
sx q[0];
rz(-0.23736048) q[0];
rz(-2.3573719) q[1];
sx q[1];
rz(-1.3619245) q[1];
sx q[1];
rz(1.4201737) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6415757) q[0];
sx q[0];
rz(-2.5529914) q[0];
sx q[0];
rz(-0.55843784) q[0];
rz(-pi) q[1];
rz(1.8397485) q[2];
sx q[2];
rz(-2.4351845) q[2];
sx q[2];
rz(-3.067467) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7569468) q[1];
sx q[1];
rz(-1.8638041) q[1];
sx q[1];
rz(-2.5811367) q[1];
rz(-pi) q[2];
rz(-2.3444207) q[3];
sx q[3];
rz(-2.4690921) q[3];
sx q[3];
rz(3.1117126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4935442) q[2];
sx q[2];
rz(-1.250114) q[2];
sx q[2];
rz(-1.2250712) q[2];
rz(1.8152292) q[3];
sx q[3];
rz(-0.95329657) q[3];
sx q[3];
rz(-1.9476604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5305283) q[0];
sx q[0];
rz(-0.079212991) q[0];
sx q[0];
rz(-0.63842574) q[0];
rz(-0.05052677) q[1];
sx q[1];
rz(-1.2082929) q[1];
sx q[1];
rz(-2.5343177) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80562185) q[0];
sx q[0];
rz(-0.50862776) q[0];
sx q[0];
rz(0.87368272) q[0];
rz(1.6462244) q[2];
sx q[2];
rz(-0.68623073) q[2];
sx q[2];
rz(-0.33160147) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5386419) q[1];
sx q[1];
rz(-1.2455229) q[1];
sx q[1];
rz(0.074217721) q[1];
rz(-1.7631986) q[3];
sx q[3];
rz(-1.7524613) q[3];
sx q[3];
rz(-0.84027973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.66024572) q[2];
sx q[2];
rz(-1.1911743) q[2];
sx q[2];
rz(-2.353239) q[2];
rz(-0.83769074) q[3];
sx q[3];
rz(-2.3205784) q[3];
sx q[3];
rz(1.2825509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7343219) q[0];
sx q[0];
rz(-2.4096074) q[0];
sx q[0];
rz(-0.5087854) q[0];
rz(1.5599627) q[1];
sx q[1];
rz(-1.0204693) q[1];
sx q[1];
rz(0.4090974) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86544308) q[0];
sx q[0];
rz(-1.7245657) q[0];
sx q[0];
rz(-1.3332644) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81014388) q[2];
sx q[2];
rz(-2.499047) q[2];
sx q[2];
rz(-1.0660397) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0171788) q[1];
sx q[1];
rz(-0.22089566) q[1];
sx q[1];
rz(-0.048103766) q[1];
rz(1.9473829) q[3];
sx q[3];
rz(-1.9344182) q[3];
sx q[3];
rz(-1.0417494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6009377) q[2];
sx q[2];
rz(-1.6161852) q[2];
sx q[2];
rz(-1.8589004) q[2];
rz(-0.013785275) q[3];
sx q[3];
rz(-2.0163592) q[3];
sx q[3];
rz(-1.8085326) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8448821) q[0];
sx q[0];
rz(-1.7322576) q[0];
sx q[0];
rz(-2.5175293) q[0];
rz(2.2804672) q[1];
sx q[1];
rz(-2.594941) q[1];
sx q[1];
rz(-3.0088967) q[1];
rz(0.3202318) q[2];
sx q[2];
rz(-1.474517) q[2];
sx q[2];
rz(1.3040645) q[2];
rz(2.224773) q[3];
sx q[3];
rz(-1.0448169) q[3];
sx q[3];
rz(-1.4257594) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
