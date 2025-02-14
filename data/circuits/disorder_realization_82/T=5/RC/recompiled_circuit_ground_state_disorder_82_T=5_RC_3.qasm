OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6717186) q[0];
sx q[0];
rz(-1.1550386) q[0];
sx q[0];
rz(1.4299097) q[0];
rz(-2.6861796) q[1];
sx q[1];
rz(-2.5502584) q[1];
sx q[1];
rz(1.1270181) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4658964) q[0];
sx q[0];
rz(-1.5592335) q[0];
sx q[0];
rz(1.3402864) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8363709) q[2];
sx q[2];
rz(-1.2392534) q[2];
sx q[2];
rz(-1.8198554) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.110636) q[1];
sx q[1];
rz(-2.1345761) q[1];
sx q[1];
rz(-2.8085676) q[1];
rz(1.5529049) q[3];
sx q[3];
rz(-1.761529) q[3];
sx q[3];
rz(-2.408556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.11898018) q[2];
sx q[2];
rz(-2.4369414) q[2];
sx q[2];
rz(1.4878147) q[2];
rz(2.7704499) q[3];
sx q[3];
rz(-2.0149588) q[3];
sx q[3];
rz(0.77905542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.30598518) q[0];
sx q[0];
rz(-1.3779409) q[0];
sx q[0];
rz(0.66407472) q[0];
rz(2.5750419) q[1];
sx q[1];
rz(-1.5357176) q[1];
sx q[1];
rz(0.18633349) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.60255) q[0];
sx q[0];
rz(-1.7930174) q[0];
sx q[0];
rz(-2.1150388) q[0];
rz(-2.8736224) q[2];
sx q[2];
rz(-1.9391962) q[2];
sx q[2];
rz(0.024193833) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1809606) q[1];
sx q[1];
rz(-0.69968984) q[1];
sx q[1];
rz(-1.8262499) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0746344) q[3];
sx q[3];
rz(-1.0506588) q[3];
sx q[3];
rz(2.0926856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.32249195) q[2];
sx q[2];
rz(-1.2995316) q[2];
sx q[2];
rz(-1.5796278) q[2];
rz(1.9942079) q[3];
sx q[3];
rz(-2.8473144) q[3];
sx q[3];
rz(-1.7779721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086562432) q[0];
sx q[0];
rz(-1.4921621) q[0];
sx q[0];
rz(0.30558875) q[0];
rz(-1.2809523) q[1];
sx q[1];
rz(-2.8260904) q[1];
sx q[1];
rz(1.618128) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20210914) q[0];
sx q[0];
rz(-1.3393219) q[0];
sx q[0];
rz(0.37136308) q[0];
x q[1];
rz(1.5554713) q[2];
sx q[2];
rz(-1.5256679) q[2];
sx q[2];
rz(-3.0588104) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9429607) q[1];
sx q[1];
rz(-0.80201282) q[1];
sx q[1];
rz(0.50870163) q[1];
x q[2];
rz(-1.2382069) q[3];
sx q[3];
rz(-2.5032515) q[3];
sx q[3];
rz(-0.35143839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.65172255) q[2];
sx q[2];
rz(-1.8946596) q[2];
sx q[2];
rz(0.22001246) q[2];
rz(-0.079719933) q[3];
sx q[3];
rz(-2.1646175) q[3];
sx q[3];
rz(-0.93130934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7109969) q[0];
sx q[0];
rz(-1.3871223) q[0];
sx q[0];
rz(0.0084477607) q[0];
rz(1.9795817) q[1];
sx q[1];
rz(-1.9879257) q[1];
sx q[1];
rz(-2.0268424) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24766391) q[0];
sx q[0];
rz(-1.8801483) q[0];
sx q[0];
rz(0.19750144) q[0];
x q[1];
rz(0.77728072) q[2];
sx q[2];
rz(-2.4156164) q[2];
sx q[2];
rz(0.41997806) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6709958) q[1];
sx q[1];
rz(-0.45491114) q[1];
sx q[1];
rz(-0.90683337) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9825806) q[3];
sx q[3];
rz(-1.1819541) q[3];
sx q[3];
rz(1.6568615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4343425) q[2];
sx q[2];
rz(-2.12314) q[2];
sx q[2];
rz(-0.84558359) q[2];
rz(2.8905458) q[3];
sx q[3];
rz(-1.8699402) q[3];
sx q[3];
rz(2.3959851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.71395981) q[0];
sx q[0];
rz(-0.41249713) q[0];
sx q[0];
rz(-1.0300256) q[0];
rz(-2.5469942) q[1];
sx q[1];
rz(-1.4232891) q[1];
sx q[1];
rz(-1.3535708) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44424442) q[0];
sx q[0];
rz(-1.5500907) q[0];
sx q[0];
rz(-0.022788825) q[0];
rz(-pi) q[1];
rz(-1.3859904) q[2];
sx q[2];
rz(-1.8880185) q[2];
sx q[2];
rz(-0.76390195) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36921484) q[1];
sx q[1];
rz(-2.2352043) q[1];
sx q[1];
rz(-2.496705) q[1];
x q[2];
rz(-0.69820893) q[3];
sx q[3];
rz(-0.72091736) q[3];
sx q[3];
rz(1.9424734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2305962) q[2];
sx q[2];
rz(-1.5065008) q[2];
sx q[2];
rz(-0.050749151) q[2];
rz(1.1132318) q[3];
sx q[3];
rz(-1.9751996) q[3];
sx q[3];
rz(-0.39065233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0942866) q[0];
sx q[0];
rz(-1.0599437) q[0];
sx q[0];
rz(-2.768709) q[0];
rz(1.3039543) q[1];
sx q[1];
rz(-1.2645489) q[1];
sx q[1];
rz(1.0252999) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4639125) q[0];
sx q[0];
rz(-0.45651528) q[0];
sx q[0];
rz(2.9677261) q[0];
x q[1];
rz(-2.9766259) q[2];
sx q[2];
rz(-1.8968762) q[2];
sx q[2];
rz(-3.0086089) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2596332) q[1];
sx q[1];
rz(-0.72364961) q[1];
sx q[1];
rz(-1.8745443) q[1];
rz(-pi) q[2];
rz(2.0399569) q[3];
sx q[3];
rz(-0.46814748) q[3];
sx q[3];
rz(-0.461138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.23124) q[2];
sx q[2];
rz(-2.9515036) q[2];
sx q[2];
rz(-1.4300038) q[2];
rz(1.6010239) q[3];
sx q[3];
rz(-2.4547596) q[3];
sx q[3];
rz(-1.6704667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1533399) q[0];
sx q[0];
rz(-1.8415469) q[0];
sx q[0];
rz(2.7100995) q[0];
rz(-1.8776114) q[1];
sx q[1];
rz(-1.5411721) q[1];
sx q[1];
rz(-2.4914609) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1558485) q[0];
sx q[0];
rz(-1.8580282) q[0];
sx q[0];
rz(1.4850856) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6986786) q[2];
sx q[2];
rz(-2.4543216) q[2];
sx q[2];
rz(2.5812347) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.11870594) q[1];
sx q[1];
rz(-2.431015) q[1];
sx q[1];
rz(1.7479244) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0286719) q[3];
sx q[3];
rz(-1.5776537) q[3];
sx q[3];
rz(2.6049861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.55107895) q[2];
sx q[2];
rz(-1.7503909) q[2];
sx q[2];
rz(3.1371269) q[2];
rz(3.1320069) q[3];
sx q[3];
rz(-1.3349345) q[3];
sx q[3];
rz(-0.34346223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9850605) q[0];
sx q[0];
rz(-0.28230202) q[0];
sx q[0];
rz(0.16047934) q[0];
rz(1.3849244) q[1];
sx q[1];
rz(-0.79912186) q[1];
sx q[1];
rz(2.8327732) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8978116) q[0];
sx q[0];
rz(-1.7574991) q[0];
sx q[0];
rz(3.0152882) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7032531) q[2];
sx q[2];
rz(-0.76073863) q[2];
sx q[2];
rz(2.1957514) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.908573) q[1];
sx q[1];
rz(-1.9601788) q[1];
sx q[1];
rz(1.1638457) q[1];
x q[2];
rz(-2.4229523) q[3];
sx q[3];
rz(-2.6860533) q[3];
sx q[3];
rz(2.9589341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.34362346) q[2];
sx q[2];
rz(-0.60911959) q[2];
sx q[2];
rz(1.5999751) q[2];
rz(-0.68495098) q[3];
sx q[3];
rz(-1.5545605) q[3];
sx q[3];
rz(1.5233013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-2.5110382) q[0];
sx q[0];
rz(-1.4383974) q[0];
sx q[0];
rz(2.5901929) q[0];
rz(0.4772056) q[1];
sx q[1];
rz(-0.91608945) q[1];
sx q[1];
rz(-2.3068857) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4906368) q[0];
sx q[0];
rz(-1.0792562) q[0];
sx q[0];
rz(-1.0440582) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4407225) q[2];
sx q[2];
rz(-1.5972023) q[2];
sx q[2];
rz(0.7757265) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1795802) q[1];
sx q[1];
rz(-2.5287345) q[1];
sx q[1];
rz(1.5108766) q[1];
x q[2];
rz(-2.1797997) q[3];
sx q[3];
rz(-2.0941705) q[3];
sx q[3];
rz(-2.9030622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7029552) q[2];
sx q[2];
rz(-2.1795887) q[2];
sx q[2];
rz(0.23923242) q[2];
rz(2.8171825) q[3];
sx q[3];
rz(-1.4374461) q[3];
sx q[3];
rz(1.5391866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5156373) q[0];
sx q[0];
rz(-0.13787585) q[0];
sx q[0];
rz(2.4249518) q[0];
rz(0.58249885) q[1];
sx q[1];
rz(-1.3526252) q[1];
sx q[1];
rz(-2.8315721) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4449883) q[0];
sx q[0];
rz(-1.2221081) q[0];
sx q[0];
rz(-0.83072386) q[0];
x q[1];
rz(2.3815709) q[2];
sx q[2];
rz(-0.54229247) q[2];
sx q[2];
rz(-0.063864683) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1056052) q[1];
sx q[1];
rz(-1.9275371) q[1];
sx q[1];
rz(2.3272661) q[1];
rz(-pi) q[2];
rz(2.5955947) q[3];
sx q[3];
rz(-1.1359478) q[3];
sx q[3];
rz(-2.7194589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6325355) q[2];
sx q[2];
rz(-2.5639503) q[2];
sx q[2];
rz(1.4198111) q[2];
rz(-1.696473) q[3];
sx q[3];
rz(-2.4308379) q[3];
sx q[3];
rz(-2.3783309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9025018) q[0];
sx q[0];
rz(-0.9466753) q[0];
sx q[0];
rz(-1.7539903) q[0];
rz(-1.6613962) q[1];
sx q[1];
rz(-1.4085242) q[1];
sx q[1];
rz(-2.9396802) q[1];
rz(1.4954064) q[2];
sx q[2];
rz(-1.8373377) q[2];
sx q[2];
rz(-1.1129825) q[2];
rz(0.46764163) q[3];
sx q[3];
rz(-2.1976039) q[3];
sx q[3];
rz(-0.80898367) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
