OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.00081113022) q[0];
sx q[0];
rz(-0.82511628) q[0];
sx q[0];
rz(0.42186475) q[0];
rz(2.3946664) q[1];
sx q[1];
rz(-0.59023017) q[1];
sx q[1];
rz(-2.3205369) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8418772) q[0];
sx q[0];
rz(-1.4869191) q[0];
sx q[0];
rz(1.5257349) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9157588) q[2];
sx q[2];
rz(-1.7712812) q[2];
sx q[2];
rz(-3.0124913) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.62748903) q[1];
sx q[1];
rz(-1.8541876) q[1];
sx q[1];
rz(-1.9691739) q[1];
x q[2];
rz(1.9299292) q[3];
sx q[3];
rz(-0.36595038) q[3];
sx q[3];
rz(2.8173878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.39516285) q[2];
sx q[2];
rz(-0.91150993) q[2];
sx q[2];
rz(-2.5968623) q[2];
rz(-1.644545) q[3];
sx q[3];
rz(-1.112273) q[3];
sx q[3];
rz(-1.8018319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6853952) q[0];
sx q[0];
rz(-2.39769) q[0];
sx q[0];
rz(2.5681382) q[0];
rz(0.042512976) q[1];
sx q[1];
rz(-1.7492234) q[1];
sx q[1];
rz(-1.858985) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1853582) q[0];
sx q[0];
rz(-1.6969674) q[0];
sx q[0];
rz(-1.5500359) q[0];
x q[1];
rz(2.8473275) q[2];
sx q[2];
rz(-1.3252826) q[2];
sx q[2];
rz(-2.9835999) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.01067082) q[1];
sx q[1];
rz(-1.0464484) q[1];
sx q[1];
rz(-1.4078543) q[1];
rz(0.62915985) q[3];
sx q[3];
rz(-1.0608309) q[3];
sx q[3];
rz(-0.30889084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.87806988) q[2];
sx q[2];
rz(-1.9922549) q[2];
sx q[2];
rz(-0.44352201) q[2];
rz(-1.589132) q[3];
sx q[3];
rz(-0.3905206) q[3];
sx q[3];
rz(0.21903567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012101128) q[0];
sx q[0];
rz(-2.2624367) q[0];
sx q[0];
rz(-0.40170676) q[0];
rz(1.4178287) q[1];
sx q[1];
rz(-0.44770733) q[1];
sx q[1];
rz(2.0254501) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4075466) q[0];
sx q[0];
rz(-0.73483682) q[0];
sx q[0];
rz(-0.46291344) q[0];
rz(-pi) q[1];
rz(2.9878408) q[2];
sx q[2];
rz(-0.88482403) q[2];
sx q[2];
rz(-1.5123715) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1765343) q[1];
sx q[1];
rz(-1.1557587) q[1];
sx q[1];
rz(-1.0227636) q[1];
rz(-pi) q[2];
rz(-0.12401993) q[3];
sx q[3];
rz(-2.5895666) q[3];
sx q[3];
rz(-0.13904143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.22905722) q[2];
sx q[2];
rz(-0.29650649) q[2];
sx q[2];
rz(-2.1336446) q[2];
rz(-1.6352765) q[3];
sx q[3];
rz(-1.9624036) q[3];
sx q[3];
rz(-2.262871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85663831) q[0];
sx q[0];
rz(-3.1143739) q[0];
sx q[0];
rz(0.83754367) q[0];
rz(1.8800927) q[1];
sx q[1];
rz(-1.3976169) q[1];
sx q[1];
rz(1.8998442) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.094114) q[0];
sx q[0];
rz(-2.6015208) q[0];
sx q[0];
rz(-2.5940226) q[0];
rz(0.9622039) q[2];
sx q[2];
rz(-1.7458916) q[2];
sx q[2];
rz(2.1661557) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.11644046) q[1];
sx q[1];
rz(-2.5106627) q[1];
sx q[1];
rz(-2.4063944) q[1];
rz(-pi) q[2];
rz(0.37184985) q[3];
sx q[3];
rz(-1.0222553) q[3];
sx q[3];
rz(-1.8601102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1602829) q[2];
sx q[2];
rz(-2.1872988) q[2];
sx q[2];
rz(2.1295638) q[2];
rz(2.1435598) q[3];
sx q[3];
rz(-2.4055552) q[3];
sx q[3];
rz(-1.4724822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07974682) q[0];
sx q[0];
rz(-2.340402) q[0];
sx q[0];
rz(-2.8629942) q[0];
rz(2.8246236) q[1];
sx q[1];
rz(-1.3895037) q[1];
sx q[1];
rz(0.46151361) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6742548) q[0];
sx q[0];
rz(-1.6188356) q[0];
sx q[0];
rz(1.7397299) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1198844) q[2];
sx q[2];
rz(-2.4442992) q[2];
sx q[2];
rz(-1.0233699) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6205004) q[1];
sx q[1];
rz(-0.25687718) q[1];
sx q[1];
rz(1.1939373) q[1];
rz(-0.096624537) q[3];
sx q[3];
rz(-2.1139675) q[3];
sx q[3];
rz(-0.66978776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3658112) q[2];
sx q[2];
rz(-1.3373988) q[2];
sx q[2];
rz(-2.865045) q[2];
rz(-0.60025275) q[3];
sx q[3];
rz(-2.4252031) q[3];
sx q[3];
rz(-0.53441179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9944331) q[0];
sx q[0];
rz(-0.58013791) q[0];
sx q[0];
rz(-2.4603727) q[0];
rz(-2.6874806) q[1];
sx q[1];
rz(-1.9845767) q[1];
sx q[1];
rz(0.62201321) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9227064) q[0];
sx q[0];
rz(-1.6657296) q[0];
sx q[0];
rz(-0.19849284) q[0];
rz(0.76592457) q[2];
sx q[2];
rz(-0.799338) q[2];
sx q[2];
rz(2.9098985) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6890735) q[1];
sx q[1];
rz(-1.9770685) q[1];
sx q[1];
rz(2.3157332) q[1];
x q[2];
rz(-0.92414738) q[3];
sx q[3];
rz(-2.7139276) q[3];
sx q[3];
rz(0.88811737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5155718) q[2];
sx q[2];
rz(-2.5280648) q[2];
sx q[2];
rz(2.386509) q[2];
rz(-0.10609047) q[3];
sx q[3];
rz(-1.2664436) q[3];
sx q[3];
rz(0.46060002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047091529) q[0];
sx q[0];
rz(-0.61151183) q[0];
sx q[0];
rz(-3.0754572) q[0];
rz(3.0905837) q[1];
sx q[1];
rz(-2.3936733) q[1];
sx q[1];
rz(-1.8775108) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0361657) q[0];
sx q[0];
rz(-1.9308865) q[0];
sx q[0];
rz(-2.197916) q[0];
x q[1];
rz(1.8813921) q[2];
sx q[2];
rz(-0.87970886) q[2];
sx q[2];
rz(-2.4533009) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.79489975) q[1];
sx q[1];
rz(-2.3155619) q[1];
sx q[1];
rz(-0.54462437) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1751647) q[3];
sx q[3];
rz(-1.1671178) q[3];
sx q[3];
rz(0.56433557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2833726) q[2];
sx q[2];
rz(-1.1274575) q[2];
sx q[2];
rz(2.9952725) q[2];
rz(-2.537651) q[3];
sx q[3];
rz(-0.54655176) q[3];
sx q[3];
rz(-1.4124136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.716575) q[0];
sx q[0];
rz(-1.9442433) q[0];
sx q[0];
rz(3.0498411) q[0];
rz(3.0692406) q[1];
sx q[1];
rz(-1.3183343) q[1];
sx q[1];
rz(-2.4718557) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2429072) q[0];
sx q[0];
rz(-1.4353485) q[0];
sx q[0];
rz(-0.29512914) q[0];
rz(-1.7013427) q[2];
sx q[2];
rz(-1.1645082) q[2];
sx q[2];
rz(0.31995904) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.32423985) q[1];
sx q[1];
rz(-1.8756556) q[1];
sx q[1];
rz(1.7225811) q[1];
x q[2];
rz(-2.9961139) q[3];
sx q[3];
rz(-2.0091741) q[3];
sx q[3];
rz(2.9117025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.79121315) q[2];
sx q[2];
rz(-2.4243441) q[2];
sx q[2];
rz(-2.4199602) q[2];
rz(0.77389884) q[3];
sx q[3];
rz(-2.1589203) q[3];
sx q[3];
rz(-0.53988808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4197107) q[0];
sx q[0];
rz(-1.1996491) q[0];
sx q[0];
rz(2.6019959) q[0];
rz(0.78060141) q[1];
sx q[1];
rz(-1.2465979) q[1];
sx q[1];
rz(2.7194729) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1264421) q[0];
sx q[0];
rz(-1.0202304) q[0];
sx q[0];
rz(-2.2534739) q[0];
rz(-pi) q[1];
rz(-1.6306179) q[2];
sx q[2];
rz(-2.0622396) q[2];
sx q[2];
rz(-2.7697542) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6680702) q[1];
sx q[1];
rz(-2.0501185) q[1];
sx q[1];
rz(-3.098686) q[1];
rz(-pi) q[2];
x q[2];
rz(0.067764564) q[3];
sx q[3];
rz(-0.58971206) q[3];
sx q[3];
rz(0.36524352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0427986) q[2];
sx q[2];
rz(-2.1197987) q[2];
sx q[2];
rz(-0.54879028) q[2];
rz(0.15268606) q[3];
sx q[3];
rz(-2.1137674) q[3];
sx q[3];
rz(0.71167439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76950908) q[0];
sx q[0];
rz(-2.7476269) q[0];
sx q[0];
rz(-0.60419303) q[0];
rz(-1.4671885) q[1];
sx q[1];
rz(-1.7908275) q[1];
sx q[1];
rz(1.1473354) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1579188) q[0];
sx q[0];
rz(-0.78494736) q[0];
sx q[0];
rz(0.60529937) q[0];
rz(-pi) q[1];
rz(-0.7355392) q[2];
sx q[2];
rz(-1.8076287) q[2];
sx q[2];
rz(-0.21422591) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5512451) q[1];
sx q[1];
rz(-2.6585124) q[1];
sx q[1];
rz(-3.1084968) q[1];
rz(-pi) q[2];
rz(-0.34965054) q[3];
sx q[3];
rz(-1.1929301) q[3];
sx q[3];
rz(0.99724619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.6818105) q[2];
sx q[2];
rz(-1.1638389) q[2];
sx q[2];
rz(-0.60010827) q[2];
rz(2.4188304) q[3];
sx q[3];
rz(-0.99812752) q[3];
sx q[3];
rz(2.7616937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4031274) q[0];
sx q[0];
rz(-1.6275788) q[0];
sx q[0];
rz(1.7049261) q[0];
rz(-1.5858831) q[1];
sx q[1];
rz(-1.7441505) q[1];
sx q[1];
rz(2.2081262) q[1];
rz(-0.49028291) q[2];
sx q[2];
rz(-0.96302196) q[2];
sx q[2];
rz(-1.8204126) q[2];
rz(2.4248471) q[3];
sx q[3];
rz(-0.32201736) q[3];
sx q[3];
rz(1.6797408) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
