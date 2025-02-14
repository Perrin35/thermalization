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
rz(-0.47082585) q[0];
sx q[0];
rz(3.5916632) q[0];
sx q[0];
rz(9.8150742) q[0];
rz(0.14416873) q[1];
sx q[1];
rz(-1.498797) q[1];
sx q[1];
rz(2.1013451) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9269104) q[0];
sx q[0];
rz(-1.9198863) q[0];
sx q[0];
rz(0.30544282) q[0];
x q[1];
rz(-1.7261271) q[2];
sx q[2];
rz(-1.2110146) q[2];
sx q[2];
rz(1.8632165) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2147602) q[1];
sx q[1];
rz(-2.4364987) q[1];
sx q[1];
rz(-1.8929204) q[1];
rz(-pi) q[2];
rz(1.6360498) q[3];
sx q[3];
rz(-2.3480519) q[3];
sx q[3];
rz(2.3161567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1990004) q[2];
sx q[2];
rz(-2.5358443) q[2];
sx q[2];
rz(-1.595363) q[2];
rz(2.6664901) q[3];
sx q[3];
rz(-0.64517704) q[3];
sx q[3];
rz(-1.1554385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5763181) q[0];
sx q[0];
rz(-2.1529038) q[0];
sx q[0];
rz(-0.43462547) q[0];
rz(-1.0644396) q[1];
sx q[1];
rz(-2.7480405) q[1];
sx q[1];
rz(-2.2501066) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1284244) q[0];
sx q[0];
rz(-1.1718318) q[0];
sx q[0];
rz(-1.1815039) q[0];
x q[1];
rz(2.5847748) q[2];
sx q[2];
rz(-2.3336683) q[2];
sx q[2];
rz(-2.9064532) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52917186) q[1];
sx q[1];
rz(-2.0222221) q[1];
sx q[1];
rz(0.58731095) q[1];
rz(-pi) q[2];
rz(2.2300612) q[3];
sx q[3];
rz(-1.5472104) q[3];
sx q[3];
rz(3.0167836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6330304) q[2];
sx q[2];
rz(-0.56266251) q[2];
sx q[2];
rz(-2.0821345) q[2];
rz(1.2262454) q[3];
sx q[3];
rz(-1.5303333) q[3];
sx q[3];
rz(-0.1196158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0934963) q[0];
sx q[0];
rz(-0.39830783) q[0];
sx q[0];
rz(-0.68921971) q[0];
rz(2.7293909) q[1];
sx q[1];
rz(-2.0239794) q[1];
sx q[1];
rz(-1.162792) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6434518) q[0];
sx q[0];
rz(-1.8751133) q[0];
sx q[0];
rz(2.9813719) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77670375) q[2];
sx q[2];
rz(-1.113184) q[2];
sx q[2];
rz(0.44241487) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.65022578) q[1];
sx q[1];
rz(-1.0904445) q[1];
sx q[1];
rz(-0.65746376) q[1];
rz(1.8602636) q[3];
sx q[3];
rz(-1.839387) q[3];
sx q[3];
rz(-0.1381607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.73843655) q[2];
sx q[2];
rz(-1.4102035) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23976633) q[0];
sx q[0];
rz(-0.16401839) q[0];
sx q[0];
rz(-0.30583403) q[0];
rz(2.8726874) q[1];
sx q[1];
rz(-0.79335672) q[1];
sx q[1];
rz(-2.7893524) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0687843) q[0];
sx q[0];
rz(-3.0020092) q[0];
sx q[0];
rz(-1.8603345) q[0];
rz(-pi) q[1];
rz(0.19465372) q[2];
sx q[2];
rz(-0.95370624) q[2];
sx q[2];
rz(0.84026166) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4388386) q[1];
sx q[1];
rz(-1.9620336) q[1];
sx q[1];
rz(0.66968285) q[1];
x q[2];
rz(1.3237185) q[3];
sx q[3];
rz(-2.3436147) q[3];
sx q[3];
rz(1.9664604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.29117808) q[2];
sx q[2];
rz(-1.6213657) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51036924) q[0];
sx q[0];
rz(-0.33741697) q[0];
sx q[0];
rz(1.1153197) q[0];
rz(1.2315617) q[1];
sx q[1];
rz(-1.7111338) q[1];
sx q[1];
rz(3.0273052) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5530722) q[0];
sx q[0];
rz(-2.3686396) q[0];
sx q[0];
rz(-2.3585178) q[0];
x q[1];
rz(1.2670934) q[2];
sx q[2];
rz(-1.3044918) q[2];
sx q[2];
rz(1.5271306) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.28262269) q[1];
sx q[1];
rz(-0.77853528) q[1];
sx q[1];
rz(2.5384739) q[1];
x q[2];
rz(1.565236) q[3];
sx q[3];
rz(-2.2768124) q[3];
sx q[3];
rz(1.8354285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.085122434) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1591448) q[0];
sx q[0];
rz(-1.5626937) q[0];
sx q[0];
rz(0.53318095) q[0];
rz(-1.606733) q[1];
sx q[1];
rz(-0.77203647) q[1];
sx q[1];
rz(0.74347043) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6966382) q[0];
sx q[0];
rz(-1.4916723) q[0];
sx q[0];
rz(0.29283578) q[0];
rz(-pi) q[1];
rz(2.1155223) q[2];
sx q[2];
rz(-1.3711978) q[2];
sx q[2];
rz(-0.064909086) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0452776) q[1];
sx q[1];
rz(-1.320134) q[1];
sx q[1];
rz(2.4947007) q[1];
x q[2];
rz(2.8072186) q[3];
sx q[3];
rz(-0.97914234) q[3];
sx q[3];
rz(-0.95229545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6712436) q[2];
sx q[2];
rz(-2.2942784) q[2];
sx q[2];
rz(-0.49873763) q[2];
rz(-0.70164743) q[3];
sx q[3];
rz(-2.4528153) q[3];
sx q[3];
rz(0.65830314) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7094803) q[0];
sx q[0];
rz(-1.9146336) q[0];
sx q[0];
rz(3.0416601) q[0];
rz(1.8607032) q[1];
sx q[1];
rz(-0.51191267) q[1];
sx q[1];
rz(2.4023712) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5422338) q[0];
sx q[0];
rz(-1.1809491) q[0];
sx q[0];
rz(2.705036) q[0];
rz(-pi) q[1];
rz(-2.0002236) q[2];
sx q[2];
rz(-1.5985367) q[2];
sx q[2];
rz(2.3703379) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.92956738) q[1];
sx q[1];
rz(-1.1788569) q[1];
sx q[1];
rz(0.93393109) q[1];
rz(-3.1312814) q[3];
sx q[3];
rz(-0.80439321) q[3];
sx q[3];
rz(-2.2571486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40333906) q[2];
sx q[2];
rz(-0.6311987) q[2];
sx q[2];
rz(1.6048019) q[2];
rz(1.4295476) q[3];
sx q[3];
rz(-1.4934544) q[3];
sx q[3];
rz(-0.79609377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0159863) q[0];
sx q[0];
rz(-0.37785372) q[0];
sx q[0];
rz(-2.0269537) q[0];
rz(-0.58427018) q[1];
sx q[1];
rz(-2.22157) q[1];
sx q[1];
rz(0.11411962) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18543359) q[0];
sx q[0];
rz(-1.8562278) q[0];
sx q[0];
rz(2.5117842) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.592351) q[2];
sx q[2];
rz(-2.261544) q[2];
sx q[2];
rz(-0.88976414) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.32875672) q[1];
sx q[1];
rz(-2.2362773) q[1];
sx q[1];
rz(2.6129641) q[1];
x q[2];
rz(1.1646284) q[3];
sx q[3];
rz(-0.28985786) q[3];
sx q[3];
rz(-1.4508307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.70322651) q[2];
sx q[2];
rz(-0.72202903) q[2];
sx q[2];
rz(0.88877338) q[2];
rz(-1.599954) q[3];
sx q[3];
rz(-1.9346574) q[3];
sx q[3];
rz(-2.3204939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0515161) q[0];
sx q[0];
rz(-2.8335644) q[0];
sx q[0];
rz(0.28710452) q[0];
rz(-1.2039394) q[1];
sx q[1];
rz(-0.86910373) q[1];
sx q[1];
rz(1.513011) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6490895) q[0];
sx q[0];
rz(-1.2375298) q[0];
sx q[0];
rz(-1.6514342) q[0];
x q[1];
rz(2.6937655) q[2];
sx q[2];
rz(-0.97140233) q[2];
sx q[2];
rz(-2.2702655) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.87344681) q[1];
sx q[1];
rz(-1.6153187) q[1];
sx q[1];
rz(2.4067307) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1407981) q[3];
sx q[3];
rz(-1.5677139) q[3];
sx q[3];
rz(2.2828988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.63605753) q[2];
sx q[2];
rz(-1.9745741) q[2];
sx q[2];
rz(2.8324845) q[2];
rz(-1.2067893) q[3];
sx q[3];
rz(-2.759582) q[3];
sx q[3];
rz(0.25930723) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1768271) q[0];
sx q[0];
rz(-2.4137156) q[0];
sx q[0];
rz(-2.4859909) q[0];
rz(2.5128095) q[1];
sx q[1];
rz(-1.4097593) q[1];
sx q[1];
rz(1.383673) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1882598) q[0];
sx q[0];
rz(-1.6142705) q[0];
sx q[0];
rz(3.093958) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6864917) q[2];
sx q[2];
rz(-2.2544207) q[2];
sx q[2];
rz(-0.065082642) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.261499) q[1];
sx q[1];
rz(-1.7200909) q[1];
sx q[1];
rz(-1.5525251) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8877533) q[3];
sx q[3];
rz(-0.51341265) q[3];
sx q[3];
rz(-1.1644582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.76929602) q[2];
sx q[2];
rz(-1.5089704) q[2];
sx q[2];
rz(0.24202913) q[2];
rz(2.5567143) q[3];
sx q[3];
rz(-2.4681028) q[3];
sx q[3];
rz(0.52403319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.723421) q[0];
sx q[0];
rz(-2.5192498) q[0];
sx q[0];
rz(2.584516) q[0];
rz(0.88049018) q[1];
sx q[1];
rz(-1.7387895) q[1];
sx q[1];
rz(-2.1008076) q[1];
rz(0.006420709) q[2];
sx q[2];
rz(-1.3646135) q[2];
sx q[2];
rz(-2.7529181) q[2];
rz(1.4222894) q[3];
sx q[3];
rz(-2.5703493) q[3];
sx q[3];
rz(-0.66019365) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
