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
rz(2.6707668) q[0];
sx q[0];
rz(-0.45007053) q[0];
sx q[0];
rz(-0.39029628) q[0];
rz(0.14416873) q[1];
sx q[1];
rz(-1.498797) q[1];
sx q[1];
rz(2.1013451) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9594355) q[0];
sx q[0];
rz(-2.6819026) q[0];
sx q[0];
rz(-0.88031405) q[0];
x q[1];
rz(2.7778012) q[2];
sx q[2];
rz(-1.4254839) q[2];
sx q[2];
rz(2.7940968) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7465248) q[1];
sx q[1];
rz(-1.3641502) q[1];
sx q[1];
rz(0.8916436) q[1];
rz(-pi) q[2];
x q[2];
rz(0.066180996) q[3];
sx q[3];
rz(-2.3621763) q[3];
sx q[3];
rz(-2.4090637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.9425923) q[2];
sx q[2];
rz(-2.5358443) q[2];
sx q[2];
rz(-1.5462297) q[2];
rz(0.47510251) q[3];
sx q[3];
rz(-0.64517704) q[3];
sx q[3];
rz(-1.9861541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5652745) q[0];
sx q[0];
rz(-2.1529038) q[0];
sx q[0];
rz(-0.43462547) q[0];
rz(-2.077153) q[1];
sx q[1];
rz(-0.39355215) q[1];
sx q[1];
rz(0.89148608) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28434174) q[0];
sx q[0];
rz(-1.9280757) q[0];
sx q[0];
rz(-2.7140359) q[0];
rz(-2.5847748) q[2];
sx q[2];
rz(-0.80792431) q[2];
sx q[2];
rz(-2.9064532) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.621755) q[1];
sx q[1];
rz(-0.72417604) q[1];
sx q[1];
rz(-2.4228079) q[1];
rz(3.1117571) q[3];
sx q[3];
rz(-2.2298457) q[3];
sx q[3];
rz(1.713879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.50856227) q[2];
sx q[2];
rz(-0.56266251) q[2];
sx q[2];
rz(2.0821345) q[2];
rz(-1.2262454) q[3];
sx q[3];
rz(-1.5303333) q[3];
sx q[3];
rz(-3.0219769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0934963) q[0];
sx q[0];
rz(-0.39830783) q[0];
sx q[0];
rz(-0.68921971) q[0];
rz(-0.41220176) q[1];
sx q[1];
rz(-2.0239794) q[1];
sx q[1];
rz(1.9788007) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9927558) q[0];
sx q[0];
rz(-2.798838) q[0];
sx q[0];
rz(-1.1008016) q[0];
rz(-pi) q[1];
rz(2.1751129) q[2];
sx q[2];
rz(-2.2508224) q[2];
sx q[2];
rz(0.71875188) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.38123576) q[1];
sx q[1];
rz(-0.79264005) q[1];
sx q[1];
rz(2.4355678) q[1];
rz(-pi) q[2];
rz(-0.80358867) q[3];
sx q[3];
rz(-2.7492964) q[3];
sx q[3];
rz(0.70480777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4031561) q[2];
sx q[2];
rz(-1.7313892) q[2];
sx q[2];
rz(2.9031244) q[2];
rz(-0.23369914) q[3];
sx q[3];
rz(-0.77478474) q[3];
sx q[3];
rz(0.15271798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23976633) q[0];
sx q[0];
rz(-0.16401839) q[0];
sx q[0];
rz(-2.8357586) q[0];
rz(-0.26890525) q[1];
sx q[1];
rz(-2.3482359) q[1];
sx q[1];
rz(2.7893524) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9304795) q[0];
sx q[0];
rz(-1.5310627) q[0];
sx q[0];
rz(1.7046403) q[0];
x q[1];
rz(-0.9446804) q[2];
sx q[2];
rz(-1.7292103) q[2];
sx q[2];
rz(2.2974654) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.70275408) q[1];
sx q[1];
rz(-1.9620336) q[1];
sx q[1];
rz(-2.4719098) q[1];
rz(1.3237185) q[3];
sx q[3];
rz(-2.3436147) q[3];
sx q[3];
rz(-1.1751323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8504146) q[2];
sx q[2];
rz(-1.520227) q[2];
sx q[2];
rz(0.16695437) q[2];
rz(-0.11635612) q[3];
sx q[3];
rz(-2.6226624) q[3];
sx q[3];
rz(1.8714582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6312234) q[0];
sx q[0];
rz(-0.33741697) q[0];
sx q[0];
rz(1.1153197) q[0];
rz(1.2315617) q[1];
sx q[1];
rz(-1.7111338) q[1];
sx q[1];
rz(-0.11428741) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63681215) q[0];
sx q[0];
rz(-2.0858602) q[0];
sx q[0];
rz(-2.5367141) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2670934) q[2];
sx q[2];
rz(-1.8371008) q[2];
sx q[2];
rz(1.5271306) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.28262269) q[1];
sx q[1];
rz(-0.77853528) q[1];
sx q[1];
rz(-0.60311879) q[1];
x q[2];
rz(0.70602368) q[3];
sx q[3];
rz(-1.5750275) q[3];
sx q[3];
rz(-2.8805681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0564702) q[2];
sx q[2];
rz(-1.2135999) q[2];
sx q[2];
rz(1.6950133) q[2];
rz(-0.96212402) q[3];
sx q[3];
rz(-2.646793) q[3];
sx q[3];
rz(1.3303293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1591448) q[0];
sx q[0];
rz(-1.5626937) q[0];
sx q[0];
rz(2.6084117) q[0];
rz(1.606733) q[1];
sx q[1];
rz(-2.3695562) q[1];
sx q[1];
rz(0.74347043) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4449544) q[0];
sx q[0];
rz(-1.4916723) q[0];
sx q[0];
rz(-0.29283578) q[0];
rz(-pi) q[1];
x q[1];
rz(1.942988) q[2];
sx q[2];
rz(-0.57665885) q[2];
sx q[2];
rz(1.8221591) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4308692) q[1];
sx q[1];
rz(-2.1942881) q[1];
sx q[1];
rz(-1.2602978) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33437403) q[3];
sx q[3];
rz(-0.97914234) q[3];
sx q[3];
rz(0.95229545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6712436) q[2];
sx q[2];
rz(-0.84731421) q[2];
sx q[2];
rz(2.642855) q[2];
rz(-2.4399452) q[3];
sx q[3];
rz(-0.68877733) q[3];
sx q[3];
rz(0.65830314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43211234) q[0];
sx q[0];
rz(-1.9146336) q[0];
sx q[0];
rz(-0.099932583) q[0];
rz(1.2808895) q[1];
sx q[1];
rz(-0.51191267) q[1];
sx q[1];
rz(0.73922149) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59935887) q[0];
sx q[0];
rz(-1.9606435) q[0];
sx q[0];
rz(-0.43655661) q[0];
rz(2.0002236) q[2];
sx q[2];
rz(-1.5985367) q[2];
sx q[2];
rz(0.77125473) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2249965) q[1];
sx q[1];
rz(-2.1526622) q[1];
sx q[1];
rz(-0.47486979) q[1];
x q[2];
rz(-1.5815063) q[3];
sx q[3];
rz(-2.3751343) q[3];
sx q[3];
rz(-0.86957726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7382536) q[2];
sx q[2];
rz(-0.6311987) q[2];
sx q[2];
rz(1.5367907) q[2];
rz(1.712045) q[3];
sx q[3];
rz(-1.4934544) q[3];
sx q[3];
rz(0.79609377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-1.0159863) q[0];
sx q[0];
rz(-2.7637389) q[0];
sx q[0];
rz(1.1146389) q[0];
rz(-0.58427018) q[1];
sx q[1];
rz(-2.22157) q[1];
sx q[1];
rz(-3.027473) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9586351) q[0];
sx q[0];
rz(-2.171423) q[0];
sx q[0];
rz(-1.22249) q[0];
rz(-pi) q[1];
rz(2.1340964) q[2];
sx q[2];
rz(-0.85342583) q[2];
sx q[2];
rz(-0.12441758) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0863959) q[1];
sx q[1];
rz(-0.82406161) q[1];
sx q[1];
rz(2.1419129) q[1];
x q[2];
rz(1.3033673) q[3];
sx q[3];
rz(-1.6839613) q[3];
sx q[3];
rz(-0.51089493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4383661) q[2];
sx q[2];
rz(-2.4195636) q[2];
sx q[2];
rz(-2.2528193) q[2];
rz(-1.5416386) q[3];
sx q[3];
rz(-1.2069353) q[3];
sx q[3];
rz(-2.3204939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.090076598) q[0];
sx q[0];
rz(-0.30802825) q[0];
sx q[0];
rz(-2.8544881) q[0];
rz(-1.2039394) q[1];
sx q[1];
rz(-2.2724889) q[1];
sx q[1];
rz(-1.513011) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8912761) q[0];
sx q[0];
rz(-0.34252942) q[0];
sx q[0];
rz(2.9129759) q[0];
rz(-pi) q[1];
rz(2.6937655) q[2];
sx q[2];
rz(-0.97140233) q[2];
sx q[2];
rz(0.87132711) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2681458) q[1];
sx q[1];
rz(-1.6153187) q[1];
sx q[1];
rz(2.4067307) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8230667) q[3];
sx q[3];
rz(-0.0031832198) q[3];
sx q[3];
rz(-1.1109655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5055351) q[2];
sx q[2];
rz(-1.1670185) q[2];
sx q[2];
rz(0.30910811) q[2];
rz(-1.2067893) q[3];
sx q[3];
rz(-2.759582) q[3];
sx q[3];
rz(0.25930723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(-1.9647656) q[0];
sx q[0];
rz(-0.72787705) q[0];
sx q[0];
rz(0.65560174) q[0];
rz(2.5128095) q[1];
sx q[1];
rz(-1.7318334) q[1];
sx q[1];
rz(1.7579196) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0198251) q[0];
sx q[0];
rz(-3.077113) q[0];
sx q[0];
rz(-2.4013258) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0655772) q[2];
sx q[2];
rz(-2.3411334) q[2];
sx q[2];
rz(-2.5474974) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0023345) q[1];
sx q[1];
rz(-2.9911925) q[1];
sx q[1];
rz(0.12087442) q[1];
x q[2];
rz(0.25383935) q[3];
sx q[3];
rz(-2.62818) q[3];
sx q[3];
rz(1.9771345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76929602) q[2];
sx q[2];
rz(-1.5089704) q[2];
sx q[2];
rz(-2.8995635) q[2];
rz(2.5567143) q[3];
sx q[3];
rz(-2.4681028) q[3];
sx q[3];
rz(0.52403319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.723421) q[0];
sx q[0];
rz(-2.5192498) q[0];
sx q[0];
rz(2.584516) q[0];
rz(-2.2611025) q[1];
sx q[1];
rz(-1.7387895) q[1];
sx q[1];
rz(-2.1008076) q[1];
rz(-3.1351719) q[2];
sx q[2];
rz(-1.3646135) q[2];
sx q[2];
rz(-2.7529181) q[2];
rz(2.1370173) q[3];
sx q[3];
rz(-1.4907111) q[3];
sx q[3];
rz(0.78540594) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
