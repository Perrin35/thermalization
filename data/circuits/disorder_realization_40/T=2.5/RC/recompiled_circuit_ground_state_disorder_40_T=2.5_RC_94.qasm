OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3482398) q[0];
sx q[0];
rz(4.6849522) q[0];
sx q[0];
rz(7.9230928) q[0];
rz(1.7465254) q[1];
sx q[1];
rz(-2.7434064) q[1];
sx q[1];
rz(-0.15287486) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9119551) q[0];
sx q[0];
rz(-0.16235936) q[0];
sx q[0];
rz(2.1549757) q[0];
rz(-pi) q[1];
rz(-0.35043535) q[2];
sx q[2];
rz(-1.4154684) q[2];
sx q[2];
rz(0.84463476) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6907883) q[1];
sx q[1];
rz(-1.6310454) q[1];
sx q[1];
rz(-1.7598399) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46061116) q[3];
sx q[3];
rz(-0.95129993) q[3];
sx q[3];
rz(-1.7509489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.07831002) q[2];
sx q[2];
rz(-1.4274884) q[2];
sx q[2];
rz(-2.5417292) q[2];
rz(2.6317224) q[3];
sx q[3];
rz(-2.8190835) q[3];
sx q[3];
rz(-2.9738284) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.584499) q[0];
sx q[0];
rz(-2.2826513) q[0];
sx q[0];
rz(2.9993045) q[0];
rz(1.8498869) q[1];
sx q[1];
rz(-2.6085491) q[1];
sx q[1];
rz(-0.60089111) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73616926) q[0];
sx q[0];
rz(-1.7975419) q[0];
sx q[0];
rz(-1.1419683) q[0];
rz(-pi) q[1];
rz(-0.63670701) q[2];
sx q[2];
rz(-1.2629303) q[2];
sx q[2];
rz(-2.740048) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8411202) q[1];
sx q[1];
rz(-0.71492793) q[1];
sx q[1];
rz(-0.27473533) q[1];
x q[2];
rz(2.8138312) q[3];
sx q[3];
rz(-2.6034688) q[3];
sx q[3];
rz(-1.6445352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.94704023) q[2];
sx q[2];
rz(-1.1636473) q[2];
sx q[2];
rz(0.8645424) q[2];
rz(-2.7583127) q[3];
sx q[3];
rz(-2.7195103) q[3];
sx q[3];
rz(0.88467902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8466723) q[0];
sx q[0];
rz(-2.8546794) q[0];
sx q[0];
rz(-2.3160146) q[0];
rz(-0.83746743) q[1];
sx q[1];
rz(-1.0191963) q[1];
sx q[1];
rz(0.79140633) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.237903) q[0];
sx q[0];
rz(-2.3775953) q[0];
sx q[0];
rz(-1.6873515) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7010274) q[2];
sx q[2];
rz(-1.669906) q[2];
sx q[2];
rz(2.2203022) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1371932) q[1];
sx q[1];
rz(-1.2257842) q[1];
sx q[1];
rz(0.54917224) q[1];
x q[2];
rz(0.84816459) q[3];
sx q[3];
rz(-2.5301683) q[3];
sx q[3];
rz(1.320672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5241663) q[2];
sx q[2];
rz(-2.7238621) q[2];
sx q[2];
rz(1.2512655) q[2];
rz(-1.1795801) q[3];
sx q[3];
rz(-2.0056632) q[3];
sx q[3];
rz(-2.484926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0254822) q[0];
sx q[0];
rz(-1.9132834) q[0];
sx q[0];
rz(-2.3226698) q[0];
rz(-3.0602835) q[1];
sx q[1];
rz(-2.5550877) q[1];
sx q[1];
rz(-0.78048817) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9355192) q[0];
sx q[0];
rz(-1.9447826) q[0];
sx q[0];
rz(-1.5078578) q[0];
x q[1];
rz(0.92018668) q[2];
sx q[2];
rz(-1.4942385) q[2];
sx q[2];
rz(2.6849942) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3230715) q[1];
sx q[1];
rz(-1.8667606) q[1];
sx q[1];
rz(0.37533111) q[1];
rz(-0.31663043) q[3];
sx q[3];
rz(-1.1854611) q[3];
sx q[3];
rz(-1.6485293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9235733) q[2];
sx q[2];
rz(-1.4998481) q[2];
sx q[2];
rz(1.2884595) q[2];
rz(-1.6740547) q[3];
sx q[3];
rz(-0.54687423) q[3];
sx q[3];
rz(2.7047777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.9913637) q[0];
sx q[0];
rz(-0.65058351) q[0];
sx q[0];
rz(-0.22317602) q[0];
rz(-0.079809345) q[1];
sx q[1];
rz(-0.97750074) q[1];
sx q[1];
rz(-0.83597437) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8337934) q[0];
sx q[0];
rz(-0.5029486) q[0];
sx q[0];
rz(2.0849243) q[0];
rz(-pi) q[1];
rz(-2.6537544) q[2];
sx q[2];
rz(-1.9291995) q[2];
sx q[2];
rz(-0.20939669) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.32415307) q[1];
sx q[1];
rz(-2.0803703) q[1];
sx q[1];
rz(-0.47199366) q[1];
rz(0.23578819) q[3];
sx q[3];
rz(-2.1683768) q[3];
sx q[3];
rz(1.2553297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.23433805) q[2];
sx q[2];
rz(-0.85504389) q[2];
sx q[2];
rz(-1.3912531) q[2];
rz(1.2587345) q[3];
sx q[3];
rz(-1.3842868) q[3];
sx q[3];
rz(2.7425227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2788972) q[0];
sx q[0];
rz(-2.6628222) q[0];
sx q[0];
rz(-0.71267772) q[0];
rz(1.124294) q[1];
sx q[1];
rz(-1.5702039) q[1];
sx q[1];
rz(2.1052776) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22651214) q[0];
sx q[0];
rz(-2.4944502) q[0];
sx q[0];
rz(2.9467366) q[0];
rz(-pi) q[1];
rz(0.89547248) q[2];
sx q[2];
rz(-2.7717441) q[2];
sx q[2];
rz(-2.0167882) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8935298) q[1];
sx q[1];
rz(-1.7772733) q[1];
sx q[1];
rz(0.96915396) q[1];
rz(-2.2676663) q[3];
sx q[3];
rz(-1.7146661) q[3];
sx q[3];
rz(1.2269536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0599646) q[2];
sx q[2];
rz(-2.5689503) q[2];
sx q[2];
rz(2.6247978) q[2];
rz(-1.9948657) q[3];
sx q[3];
rz(-2.1490993) q[3];
sx q[3];
rz(0.049588047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0674755) q[0];
sx q[0];
rz(-0.68454409) q[0];
sx q[0];
rz(1.462498) q[0];
rz(2.0430203) q[1];
sx q[1];
rz(-1.699828) q[1];
sx q[1];
rz(2.3625653) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50139105) q[0];
sx q[0];
rz(-3.0872652) q[0];
sx q[0];
rz(-3.1110248) q[0];
rz(-2.2279694) q[2];
sx q[2];
rz(-0.69013287) q[2];
sx q[2];
rz(-2.4724378) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.15493079) q[1];
sx q[1];
rz(-0.26954654) q[1];
sx q[1];
rz(-0.69828548) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1356462) q[3];
sx q[3];
rz(-2.3823934) q[3];
sx q[3];
rz(0.58461207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.74759185) q[2];
sx q[2];
rz(-2.438075) q[2];
sx q[2];
rz(0.35959378) q[2];
rz(-2.8408585) q[3];
sx q[3];
rz(-2.3263003) q[3];
sx q[3];
rz(-2.1451779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71036285) q[0];
sx q[0];
rz(-1.2082986) q[0];
sx q[0];
rz(-2.8408458) q[0];
rz(-2.6077479) q[1];
sx q[1];
rz(-2.0678935) q[1];
sx q[1];
rz(-3.0278382) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14989242) q[0];
sx q[0];
rz(-1.0731369) q[0];
sx q[0];
rz(-2.0927621) q[0];
rz(-2.5328296) q[2];
sx q[2];
rz(-1.5396313) q[2];
sx q[2];
rz(-1.6682012) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2054322) q[1];
sx q[1];
rz(-2.0179085) q[1];
sx q[1];
rz(2.4287249) q[1];
rz(-2.9748125) q[3];
sx q[3];
rz(-1.9403701) q[3];
sx q[3];
rz(-2.5931234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.45812312) q[2];
sx q[2];
rz(-0.14294954) q[2];
sx q[2];
rz(2.6911823) q[2];
rz(0.9149552) q[3];
sx q[3];
rz(-2.2136642) q[3];
sx q[3];
rz(-1.3373059) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27110711) q[0];
sx q[0];
rz(-0.73942375) q[0];
sx q[0];
rz(-0.41573218) q[0];
rz(2.7102026) q[1];
sx q[1];
rz(-0.58512551) q[1];
sx q[1];
rz(2.364667) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8857461) q[0];
sx q[0];
rz(-1.5836459) q[0];
sx q[0];
rz(1.6240516) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3007352) q[2];
sx q[2];
rz(-1.0094202) q[2];
sx q[2];
rz(-2.4955668) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.61762702) q[1];
sx q[1];
rz(-1.0203779) q[1];
sx q[1];
rz(1.074076) q[1];
rz(-2.611972) q[3];
sx q[3];
rz(-2.6130015) q[3];
sx q[3];
rz(1.1666672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1315769) q[2];
sx q[2];
rz(-0.78993979) q[2];
sx q[2];
rz(1.5124403) q[2];
rz(-0.43859628) q[3];
sx q[3];
rz(-0.83493835) q[3];
sx q[3];
rz(2.6252873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55385357) q[0];
sx q[0];
rz(-0.53936154) q[0];
sx q[0];
rz(-1.825765) q[0];
rz(-0.081789628) q[1];
sx q[1];
rz(-1.1809228) q[1];
sx q[1];
rz(-1.2710424) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1493001) q[0];
sx q[0];
rz(-1.1470511) q[0];
sx q[0];
rz(-0.032024327) q[0];
x q[1];
rz(-0.84065915) q[2];
sx q[2];
rz(-1.5224384) q[2];
sx q[2];
rz(0.73634597) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0942022) q[1];
sx q[1];
rz(-1.8592111) q[1];
sx q[1];
rz(2.3408195) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5457235) q[3];
sx q[3];
rz(-1.4891324) q[3];
sx q[3];
rz(-1.1829528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.79909331) q[2];
sx q[2];
rz(-1.7010331) q[2];
sx q[2];
rz(-1.052617) q[2];
rz(1.8140225) q[3];
sx q[3];
rz(-1.1955639) q[3];
sx q[3];
rz(0.50935203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4159745) q[0];
sx q[0];
rz(-1.7617891) q[0];
sx q[0];
rz(-1.2431385) q[0];
rz(1.7983408) q[1];
sx q[1];
rz(-0.46011283) q[1];
sx q[1];
rz(-2.7759001) q[1];
rz(2.9467498) q[2];
sx q[2];
rz(-1.101782) q[2];
sx q[2];
rz(-0.89141104) q[2];
rz(2.7300077) q[3];
sx q[3];
rz(-1.5829908) q[3];
sx q[3];
rz(-0.79394059) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
