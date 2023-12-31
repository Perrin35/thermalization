OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2258423) q[0];
sx q[0];
rz(-0.031210829) q[0];
sx q[0];
rz(2.6565235) q[0];
rz(-2.354061) q[1];
sx q[1];
rz(-2.1252316) q[1];
sx q[1];
rz(-2.7273942) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038923351) q[0];
sx q[0];
rz(-1.9063213) q[0];
sx q[0];
rz(1.9468007) q[0];
rz(1.0010927) q[2];
sx q[2];
rz(-1.4909407) q[2];
sx q[2];
rz(0.18730883) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84250433) q[1];
sx q[1];
rz(-0.43049225) q[1];
sx q[1];
rz(1.7249785) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26478404) q[3];
sx q[3];
rz(-1.4720819) q[3];
sx q[3];
rz(2.7717154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2177314) q[2];
sx q[2];
rz(-1.2552746) q[2];
sx q[2];
rz(3.1100173) q[2];
rz(-1.8850373) q[3];
sx q[3];
rz(-0.50285181) q[3];
sx q[3];
rz(-0.52662915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8388222) q[0];
sx q[0];
rz(-1.4571723) q[0];
sx q[0];
rz(0.1698499) q[0];
rz(-2.4376712) q[1];
sx q[1];
rz(-1.0715276) q[1];
sx q[1];
rz(0.53952113) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8115494) q[0];
sx q[0];
rz(-2.0199611) q[0];
sx q[0];
rz(-3.0898068) q[0];
rz(-pi) q[1];
rz(-0.45692921) q[2];
sx q[2];
rz(-0.77605844) q[2];
sx q[2];
rz(1.4571783) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5121465) q[1];
sx q[1];
rz(-2.0680032) q[1];
sx q[1];
rz(-2.9260103) q[1];
rz(0.86124729) q[3];
sx q[3];
rz(-2.0112787) q[3];
sx q[3];
rz(1.4351821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4743621) q[2];
sx q[2];
rz(-1.2381866) q[2];
sx q[2];
rz(0.24307069) q[2];
rz(2.4754751) q[3];
sx q[3];
rz(-0.56454286) q[3];
sx q[3];
rz(1.2438141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3617525) q[0];
sx q[0];
rz(-0.11479522) q[0];
sx q[0];
rz(-0.4483805) q[0];
rz(-1.7547296) q[1];
sx q[1];
rz(-1.9877537) q[1];
sx q[1];
rz(0.2562491) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.019779531) q[0];
sx q[0];
rz(-0.67987961) q[0];
sx q[0];
rz(-2.7133184) q[0];
x q[1];
rz(-2.0492378) q[2];
sx q[2];
rz(-1.1698327) q[2];
sx q[2];
rz(-1.5329597) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.2910034) q[1];
sx q[1];
rz(-2.2271419) q[1];
sx q[1];
rz(2.7235051) q[1];
rz(-pi) q[2];
x q[2];
rz(2.933421) q[3];
sx q[3];
rz(-2.9609207) q[3];
sx q[3];
rz(-2.4908623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3391352) q[2];
sx q[2];
rz(-2.4352303) q[2];
sx q[2];
rz(0.57470542) q[2];
rz(1.3556708) q[3];
sx q[3];
rz(-1.1698497) q[3];
sx q[3];
rz(-1.908196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.213585) q[0];
sx q[0];
rz(-1.7127697) q[0];
sx q[0];
rz(-2.8821049) q[0];
rz(-1.9909987) q[1];
sx q[1];
rz(-1.3535627) q[1];
sx q[1];
rz(-2.4096699) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4703092) q[0];
sx q[0];
rz(-1.0571612) q[0];
sx q[0];
rz(-1.0259823) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.133693) q[2];
sx q[2];
rz(-2.4198654) q[2];
sx q[2];
rz(-1.114811) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.96958292) q[1];
sx q[1];
rz(-1.3130377) q[1];
sx q[1];
rz(1.3700563) q[1];
x q[2];
rz(2.6287574) q[3];
sx q[3];
rz(-0.25179112) q[3];
sx q[3];
rz(-0.73392111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4449473) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(-0.15110061) q[2];
rz(-2.5949196) q[3];
sx q[3];
rz(-2.0918545) q[3];
sx q[3];
rz(2.7643519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54995173) q[0];
sx q[0];
rz(-2.0786091) q[0];
sx q[0];
rz(1.2623825) q[0];
rz(-1.4683912) q[1];
sx q[1];
rz(-0.60931283) q[1];
sx q[1];
rz(-0.79777065) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7724458) q[0];
sx q[0];
rz(-1.6909084) q[0];
sx q[0];
rz(-0.0025047501) q[0];
rz(-0.61879976) q[2];
sx q[2];
rz(-0.36703645) q[2];
sx q[2];
rz(2.7235081) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3623912) q[1];
sx q[1];
rz(-0.63450846) q[1];
sx q[1];
rz(-1.7290551) q[1];
x q[2];
rz(-0.4425211) q[3];
sx q[3];
rz(-1.8707152) q[3];
sx q[3];
rz(2.4181441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.795934) q[2];
sx q[2];
rz(-0.63085932) q[2];
sx q[2];
rz(-2.8395555) q[2];
rz(-1.9942412) q[3];
sx q[3];
rz(-1.4619504) q[3];
sx q[3];
rz(-2.5938477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40925947) q[0];
sx q[0];
rz(-0.091826037) q[0];
sx q[0];
rz(1.9858032) q[0];
rz(-2.0571158) q[1];
sx q[1];
rz(-0.98025727) q[1];
sx q[1];
rz(3.0715122) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2698343) q[0];
sx q[0];
rz(-1.7239128) q[0];
sx q[0];
rz(-1.397875) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1224535) q[2];
sx q[2];
rz(-2.2746804) q[2];
sx q[2];
rz(2.4920419) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8225704) q[1];
sx q[1];
rz(-0.59148568) q[1];
sx q[1];
rz(2.849008) q[1];
x q[2];
rz(-1.4962247) q[3];
sx q[3];
rz(-1.6061022) q[3];
sx q[3];
rz(-1.1535742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.30248102) q[2];
sx q[2];
rz(-1.182686) q[2];
sx q[2];
rz(-0.62136674) q[2];
rz(-1.4403884) q[3];
sx q[3];
rz(-2.6337603) q[3];
sx q[3];
rz(2.5682209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086534111) q[0];
sx q[0];
rz(-1.4165514) q[0];
sx q[0];
rz(0.85987464) q[0];
rz(-1.9372008) q[1];
sx q[1];
rz(-2.269373) q[1];
sx q[1];
rz(3.133657) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36695776) q[0];
sx q[0];
rz(-1.8263706) q[0];
sx q[0];
rz(-0.5704244) q[0];
rz(-pi) q[1];
rz(1.0646348) q[2];
sx q[2];
rz(-2.3828265) q[2];
sx q[2];
rz(0.90492349) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4850033) q[1];
sx q[1];
rz(-2.2409229) q[1];
sx q[1];
rz(-2.2336002) q[1];
x q[2];
rz(0.90585917) q[3];
sx q[3];
rz(-2.6568036) q[3];
sx q[3];
rz(1.0672027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.69592151) q[2];
sx q[2];
rz(-1.2233223) q[2];
sx q[2];
rz(2.725214) q[2];
rz(-1.3683866) q[3];
sx q[3];
rz(-1.2973283) q[3];
sx q[3];
rz(-2.2369475) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1798582) q[0];
sx q[0];
rz(-2.8680153) q[0];
sx q[0];
rz(2.7767048) q[0];
rz(0.94003135) q[1];
sx q[1];
rz(-0.53661984) q[1];
sx q[1];
rz(1.5023124) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68740326) q[0];
sx q[0];
rz(-0.25035509) q[0];
sx q[0];
rz(-0.040391163) q[0];
x q[1];
rz(-2.6697568) q[2];
sx q[2];
rz(-0.89315692) q[2];
sx q[2];
rz(-3.133528) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9598436) q[1];
sx q[1];
rz(-1.0105003) q[1];
sx q[1];
rz(0.77397857) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1542529) q[3];
sx q[3];
rz(-0.22437469) q[3];
sx q[3];
rz(2.8191872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6442948) q[2];
sx q[2];
rz(-2.6361894) q[2];
sx q[2];
rz(1.0650744) q[2];
rz(0.30125695) q[3];
sx q[3];
rz(-1.4091636) q[3];
sx q[3];
rz(1.6206954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.168468) q[0];
sx q[0];
rz(-0.081806101) q[0];
sx q[0];
rz(-2.705943) q[0];
rz(1.7565953) q[1];
sx q[1];
rz(-0.47043097) q[1];
sx q[1];
rz(-0.41697821) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57396736) q[0];
sx q[0];
rz(-1.5331393) q[0];
sx q[0];
rz(-0.91659878) q[0];
x q[1];
rz(2.142799) q[2];
sx q[2];
rz(-0.98888328) q[2];
sx q[2];
rz(-1.4243766) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0028487) q[1];
sx q[1];
rz(-1.4321623) q[1];
sx q[1];
rz(2.0744051) q[1];
rz(-pi) q[2];
rz(3.0005089) q[3];
sx q[3];
rz(-1.699563) q[3];
sx q[3];
rz(-1.5878549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0372662) q[2];
sx q[2];
rz(-1.4929079) q[2];
sx q[2];
rz(-2.9157675) q[2];
rz(2.9337692) q[3];
sx q[3];
rz(-2.4184629) q[3];
sx q[3];
rz(0.58661714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41480961) q[0];
sx q[0];
rz(-2.8746958) q[0];
sx q[0];
rz(1.5243994) q[0];
rz(-0.95364755) q[1];
sx q[1];
rz(-1.8667659) q[1];
sx q[1];
rz(-1.3226002) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029862558) q[0];
sx q[0];
rz(-2.7435281) q[0];
sx q[0];
rz(1.8882621) q[0];
rz(-pi) q[1];
rz(1.0820504) q[2];
sx q[2];
rz(-2.6419123) q[2];
sx q[2];
rz(-0.79364712) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.10510124) q[1];
sx q[1];
rz(-2.0836012) q[1];
sx q[1];
rz(2.3829616) q[1];
rz(-pi) q[2];
rz(1.3721458) q[3];
sx q[3];
rz(-2.47654) q[3];
sx q[3];
rz(2.8952451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7913197) q[2];
sx q[2];
rz(-2.0952756) q[2];
sx q[2];
rz(-1.0160758) q[2];
rz(1.2223876) q[3];
sx q[3];
rz(-0.17761579) q[3];
sx q[3];
rz(2.5861752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9983457) q[0];
sx q[0];
rz(-0.9220985) q[0];
sx q[0];
rz(-2.0621598) q[0];
rz(-1.3636419) q[1];
sx q[1];
rz(-1.2294055) q[1];
sx q[1];
rz(1.3235863) q[1];
rz(0.017756391) q[2];
sx q[2];
rz(-2.6323071) q[2];
sx q[2];
rz(1.4321362) q[2];
rz(1.3397459) q[3];
sx q[3];
rz(-1.4828724) q[3];
sx q[3];
rz(2.8433269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
