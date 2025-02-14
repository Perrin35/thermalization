OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7656443) q[0];
sx q[0];
rz(-0.41416895) q[0];
sx q[0];
rz(2.3400657) q[0];
rz(-0.0030567788) q[1];
sx q[1];
rz(-0.86441511) q[1];
sx q[1];
rz(3.0473696) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1202336) q[0];
sx q[0];
rz(-0.84478837) q[0];
sx q[0];
rz(-1.057522) q[0];
rz(-pi) q[1];
rz(-1.8434486) q[2];
sx q[2];
rz(-1.5412207) q[2];
sx q[2];
rz(2.3675413) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7677418) q[1];
sx q[1];
rz(-2.5818425) q[1];
sx q[1];
rz(2.1327726) q[1];
rz(1.5685969) q[3];
sx q[3];
rz(-1.4181976) q[3];
sx q[3];
rz(-0.27657498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.5376544) q[2];
sx q[2];
rz(-0.64017576) q[2];
sx q[2];
rz(-1.5113277) q[2];
rz(-0.18621914) q[3];
sx q[3];
rz(-2.7112466) q[3];
sx q[3];
rz(2.490624) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52361012) q[0];
sx q[0];
rz(-1.1815434) q[0];
sx q[0];
rz(0.13667662) q[0];
rz(-3.0467721) q[1];
sx q[1];
rz(-0.46351981) q[1];
sx q[1];
rz(-2.3725841) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.626132) q[0];
sx q[0];
rz(-2.7410586) q[0];
sx q[0];
rz(0.35450165) q[0];
x q[1];
rz(1.9273754) q[2];
sx q[2];
rz(-2.2327023) q[2];
sx q[2];
rz(1.8476768) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0598073) q[1];
sx q[1];
rz(-1.4084653) q[1];
sx q[1];
rz(2.3947706) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8212998) q[3];
sx q[3];
rz(-0.13938306) q[3];
sx q[3];
rz(1.2047307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39917055) q[2];
sx q[2];
rz(-0.48226446) q[2];
sx q[2];
rz(1.4313618) q[2];
rz(2.6509905) q[3];
sx q[3];
rz(-0.49121818) q[3];
sx q[3];
rz(0.77559364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2630149) q[0];
sx q[0];
rz(-2.7624625) q[0];
sx q[0];
rz(-1.0947134) q[0];
rz(-1.3021944) q[1];
sx q[1];
rz(-0.98919386) q[1];
sx q[1];
rz(3.1375695) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99056315) q[0];
sx q[0];
rz(-0.73712611) q[0];
sx q[0];
rz(-1.3746136) q[0];
rz(-pi) q[1];
rz(-2.3962069) q[2];
sx q[2];
rz(-0.89085397) q[2];
sx q[2];
rz(2.5542575) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0672863) q[1];
sx q[1];
rz(-1.8278711) q[1];
sx q[1];
rz(-1.6021614) q[1];
x q[2];
rz(-2.4278042) q[3];
sx q[3];
rz(-1.3828438) q[3];
sx q[3];
rz(-0.97602188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5502988) q[2];
sx q[2];
rz(-0.58612263) q[2];
sx q[2];
rz(0.44814056) q[2];
rz(-2.6552933) q[3];
sx q[3];
rz(-2.2794162) q[3];
sx q[3];
rz(-2.015131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9124741) q[0];
sx q[0];
rz(-1.7121226) q[0];
sx q[0];
rz(0.66977704) q[0];
rz(-1.229333) q[1];
sx q[1];
rz(-2.6826617) q[1];
sx q[1];
rz(-0.545048) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9733423) q[0];
sx q[0];
rz(-2.0297883) q[0];
sx q[0];
rz(-2.5516872) q[0];
rz(2.2322234) q[2];
sx q[2];
rz(-2.0880359) q[2];
sx q[2];
rz(-1.6632207) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.8432359) q[1];
sx q[1];
rz(-1.3562172) q[1];
sx q[1];
rz(3.0797348) q[1];
rz(-pi) q[2];
rz(1.0944975) q[3];
sx q[3];
rz(-0.83134388) q[3];
sx q[3];
rz(3.1325213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5547319) q[2];
sx q[2];
rz(-1.6394337) q[2];
sx q[2];
rz(-0.16793212) q[2];
rz(1.2083017) q[3];
sx q[3];
rz(-0.30253634) q[3];
sx q[3];
rz(-1.4471853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5912882) q[0];
sx q[0];
rz(-1.2606324) q[0];
sx q[0];
rz(3.1377129) q[0];
rz(0.30715352) q[1];
sx q[1];
rz(-2.3852564) q[1];
sx q[1];
rz(0.53877962) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2821747) q[0];
sx q[0];
rz(-0.27056405) q[0];
sx q[0];
rz(0.55563386) q[0];
rz(-0.42947763) q[2];
sx q[2];
rz(-1.6971849) q[2];
sx q[2];
rz(2.1548558) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9591892) q[1];
sx q[1];
rz(-0.99970523) q[1];
sx q[1];
rz(1.1491999) q[1];
rz(-0.49541028) q[3];
sx q[3];
rz(-0.79548478) q[3];
sx q[3];
rz(-2.0919679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6946081) q[2];
sx q[2];
rz(-1.0993404) q[2];
sx q[2];
rz(-1.2896607) q[2];
rz(-1.5223711) q[3];
sx q[3];
rz(-0.74519849) q[3];
sx q[3];
rz(-1.6600018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8412142) q[0];
sx q[0];
rz(-2.2024246) q[0];
sx q[0];
rz(-3.0389431) q[0];
rz(-2.4094021) q[1];
sx q[1];
rz(-2.3468572) q[1];
sx q[1];
rz(2.5644459) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6711376) q[0];
sx q[0];
rz(-1.5744493) q[0];
sx q[0];
rz(1.5994344) q[0];
rz(-pi) q[1];
rz(0.58697613) q[2];
sx q[2];
rz(-1.4329301) q[2];
sx q[2];
rz(1.1816813) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0476372) q[1];
sx q[1];
rz(-3.006662) q[1];
sx q[1];
rz(0.62420242) q[1];
rz(-pi) q[2];
rz(-0.34586819) q[3];
sx q[3];
rz(-2.3551999) q[3];
sx q[3];
rz(0.47084558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.131549) q[2];
sx q[2];
rz(-0.71655822) q[2];
sx q[2];
rz(1.5232167) q[2];
rz(0.067954436) q[3];
sx q[3];
rz(-2.8165635) q[3];
sx q[3];
rz(-2.3006191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.0129358) q[0];
sx q[0];
rz(-2.1067297) q[0];
sx q[0];
rz(-2.3701684) q[0];
rz(0.5386638) q[1];
sx q[1];
rz(-1.795105) q[1];
sx q[1];
rz(-1.0661941) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5911884) q[0];
sx q[0];
rz(-2.8771391) q[0];
sx q[0];
rz(-1.868684) q[0];
rz(0.23287878) q[2];
sx q[2];
rz(-1.0693197) q[2];
sx q[2];
rz(-1.6362783) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.33246189) q[1];
sx q[1];
rz(-1.5576502) q[1];
sx q[1];
rz(2.9244208) q[1];
x q[2];
rz(0.66460426) q[3];
sx q[3];
rz(-0.41575894) q[3];
sx q[3];
rz(0.98008728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6923339) q[2];
sx q[2];
rz(-0.80654311) q[2];
sx q[2];
rz(-0.44245693) q[2];
rz(0.19065204) q[3];
sx q[3];
rz(-0.72398829) q[3];
sx q[3];
rz(2.382522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(0.2451179) q[0];
sx q[0];
rz(-2.9560095) q[0];
sx q[0];
rz(-1.8316487) q[0];
rz(-0.39778057) q[1];
sx q[1];
rz(-0.82292992) q[1];
sx q[1];
rz(1.7479053) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01263857) q[0];
sx q[0];
rz(-2.0067747) q[0];
sx q[0];
rz(-2.1672744) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2692503) q[2];
sx q[2];
rz(-1.796495) q[2];
sx q[2];
rz(1.7762453) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.950395) q[1];
sx q[1];
rz(-1.1069595) q[1];
sx q[1];
rz(1.2359124) q[1];
x q[2];
rz(-3.0354064) q[3];
sx q[3];
rz(-1.7683889) q[3];
sx q[3];
rz(-1.9492409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1132249) q[2];
sx q[2];
rz(-2.4312879) q[2];
sx q[2];
rz(-0.073471546) q[2];
rz(2.4469589) q[3];
sx q[3];
rz(-1.2725384) q[3];
sx q[3];
rz(1.0276444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6756814) q[0];
sx q[0];
rz(-0.2008734) q[0];
sx q[0];
rz(-0.95941108) q[0];
rz(-0.77964669) q[1];
sx q[1];
rz(-1.0028853) q[1];
sx q[1];
rz(-2.8693105) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36764964) q[0];
sx q[0];
rz(-3.127357) q[0];
sx q[0];
rz(1.1971426) q[0];
x q[1];
rz(-0.86560005) q[2];
sx q[2];
rz(-1.8097449) q[2];
sx q[2];
rz(-1.0098863) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.694937) q[1];
sx q[1];
rz(-2.4866101) q[1];
sx q[1];
rz(-0.86465624) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97314463) q[3];
sx q[3];
rz(-2.2242507) q[3];
sx q[3];
rz(-2.0367095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.73725629) q[2];
sx q[2];
rz(-1.5055089) q[2];
sx q[2];
rz(-0.20757248) q[2];
rz(-0.10457822) q[3];
sx q[3];
rz(-2.6360377) q[3];
sx q[3];
rz(-1.6149717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.680147) q[0];
sx q[0];
rz(-1.9141645) q[0];
sx q[0];
rz(1.0567868) q[0];
rz(-2.5959004) q[1];
sx q[1];
rz(-0.97815198) q[1];
sx q[1];
rz(0.32870764) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0595039) q[0];
sx q[0];
rz(-0.57841136) q[0];
sx q[0];
rz(2.2057981) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6255837) q[2];
sx q[2];
rz(-1.4984926) q[2];
sx q[2];
rz(-2.876578) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.199177) q[1];
sx q[1];
rz(-0.34245447) q[1];
sx q[1];
rz(1.6204932) q[1];
rz(-pi) q[2];
rz(0.77158934) q[3];
sx q[3];
rz(-1.7696915) q[3];
sx q[3];
rz(-1.0071013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9823965) q[2];
sx q[2];
rz(-0.13267645) q[2];
sx q[2];
rz(-1.1245493) q[2];
rz(0.076796181) q[3];
sx q[3];
rz(-0.43282893) q[3];
sx q[3];
rz(-2.2176149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9818253) q[0];
sx q[0];
rz(-1.2564909) q[0];
sx q[0];
rz(1.5782574) q[0];
rz(0.81623296) q[1];
sx q[1];
rz(-1.4030133) q[1];
sx q[1];
rz(-1.0907008) q[1];
rz(-2.6098567) q[2];
sx q[2];
rz(-2.3219863) q[2];
sx q[2];
rz(1.815997) q[2];
rz(0.40533752) q[3];
sx q[3];
rz(-1.5673076) q[3];
sx q[3];
rz(0.81893541) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
