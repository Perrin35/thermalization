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
rz(1.243408) q[0];
sx q[0];
rz(4.7060407) q[0];
sx q[0];
rz(10.336182) q[0];
rz(-0.66943327) q[1];
sx q[1];
rz(-0.27639204) q[1];
sx q[1];
rz(-2.3120094) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2612635) q[0];
sx q[0];
rz(-1.1733049) q[0];
sx q[0];
rz(-0.3573768) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8943823) q[2];
sx q[2];
rz(-1.4805549) q[2];
sx q[2];
rz(-1.5238786) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1948283) q[1];
sx q[1];
rz(-2.2948537) q[1];
sx q[1];
rz(-1.9763234) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31949701) q[3];
sx q[3];
rz(-0.45630672) q[3];
sx q[3];
rz(-0.97684492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7726045) q[2];
sx q[2];
rz(-0.96869865) q[2];
sx q[2];
rz(1.9770835) q[2];
rz(0.72689593) q[3];
sx q[3];
rz(-1.8098857) q[3];
sx q[3];
rz(3.0708142) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1233391) q[0];
sx q[0];
rz(-0.66272074) q[0];
sx q[0];
rz(-2.8501999) q[0];
rz(-2.5413051) q[1];
sx q[1];
rz(-0.95708668) q[1];
sx q[1];
rz(2.9765863) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8931029) q[0];
sx q[0];
rz(-1.6418253) q[0];
sx q[0];
rz(-1.8336589) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1218614) q[2];
sx q[2];
rz(-1.4693207) q[2];
sx q[2];
rz(-1.5047274) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.2272039) q[1];
sx q[1];
rz(-1.640135) q[1];
sx q[1];
rz(0.68080618) q[1];
rz(-pi) q[2];
rz(-0.21221186) q[3];
sx q[3];
rz(-2.8877875) q[3];
sx q[3];
rz(-2.5830944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0773641) q[2];
sx q[2];
rz(-0.95244971) q[2];
sx q[2];
rz(2.9020818) q[2];
rz(-0.32660487) q[3];
sx q[3];
rz(-0.17290792) q[3];
sx q[3];
rz(-3.0467564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8851682) q[0];
sx q[0];
rz(-0.435985) q[0];
sx q[0];
rz(-1.5727795) q[0];
rz(-2.7478711) q[1];
sx q[1];
rz(-2.5865793) q[1];
sx q[1];
rz(-2.6655925) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8825553) q[0];
sx q[0];
rz(-0.60094423) q[0];
sx q[0];
rz(1.9447692) q[0];
x q[1];
rz(-3.1114678) q[2];
sx q[2];
rz(-1.7576412) q[2];
sx q[2];
rz(-0.060299071) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.16533599) q[1];
sx q[1];
rz(-1.1603198) q[1];
sx q[1];
rz(0.93595502) q[1];
rz(-pi) q[2];
rz(-0.24781052) q[3];
sx q[3];
rz(-1.2427287) q[3];
sx q[3];
rz(2.1472665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8126882) q[2];
sx q[2];
rz(-0.71949553) q[2];
sx q[2];
rz(2.1997814) q[2];
rz(-1.4224667) q[3];
sx q[3];
rz(-2.7673281) q[3];
sx q[3];
rz(2.4678738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068168966) q[0];
sx q[0];
rz(-2.9015559) q[0];
sx q[0];
rz(-1.6085251) q[0];
rz(0.34670058) q[1];
sx q[1];
rz(-1.0418714) q[1];
sx q[1];
rz(0.18992058) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1153206) q[0];
sx q[0];
rz(-0.85912017) q[0];
sx q[0];
rz(1.9188966) q[0];
x q[1];
rz(-2.047891) q[2];
sx q[2];
rz(-1.8162842) q[2];
sx q[2];
rz(-2.7141822) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4081215) q[1];
sx q[1];
rz(-1.6279393) q[1];
sx q[1];
rz(0.96459324) q[1];
rz(-2.3307822) q[3];
sx q[3];
rz(-0.56104198) q[3];
sx q[3];
rz(-1.1277792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3612264) q[2];
sx q[2];
rz(-2.2915338) q[2];
sx q[2];
rz(0.38193199) q[2];
rz(-0.10968883) q[3];
sx q[3];
rz(-1.0332801) q[3];
sx q[3];
rz(-3.0221353) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0383976) q[0];
sx q[0];
rz(-1.9586451) q[0];
sx q[0];
rz(2.2947445) q[0];
rz(0.13941828) q[1];
sx q[1];
rz(-0.81652313) q[1];
sx q[1];
rz(0.74904186) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8183621) q[0];
sx q[0];
rz(-1.7182171) q[0];
sx q[0];
rz(2.4458103) q[0];
rz(-pi) q[1];
rz(-0.48495082) q[2];
sx q[2];
rz(-1.6963939) q[2];
sx q[2];
rz(1.4664949) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3917424) q[1];
sx q[1];
rz(-0.23255177) q[1];
sx q[1];
rz(2.560023) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0723125) q[3];
sx q[3];
rz(-1.4677047) q[3];
sx q[3];
rz(-0.13971288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2060812) q[2];
sx q[2];
rz(-1.8847909) q[2];
sx q[2];
rz(2.0743267) q[2];
rz(2.4308128) q[3];
sx q[3];
rz(-1.658541) q[3];
sx q[3];
rz(-1.425364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34273219) q[0];
sx q[0];
rz(-1.6628168) q[0];
sx q[0];
rz(0.96858281) q[0];
rz(-2.4505278) q[1];
sx q[1];
rz(-1.3422809) q[1];
sx q[1];
rz(1.8341281) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2742776) q[0];
sx q[0];
rz(-2.4558892) q[0];
sx q[0];
rz(-0.4088485) q[0];
x q[1];
rz(0.5796104) q[2];
sx q[2];
rz(-0.34380772) q[2];
sx q[2];
rz(2.0909617) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.56388748) q[1];
sx q[1];
rz(-1.1839377) q[1];
sx q[1];
rz(0.056864077) q[1];
x q[2];
rz(-0.97882459) q[3];
sx q[3];
rz(-2.8205964) q[3];
sx q[3];
rz(-2.5227566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.4751733) q[2];
sx q[2];
rz(-0.23491493) q[2];
sx q[2];
rz(1.2972181) q[2];
rz(1.3951067) q[3];
sx q[3];
rz(-1.8078943) q[3];
sx q[3];
rz(-1.5313799) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.34201) q[0];
sx q[0];
rz(-2.3260249) q[0];
sx q[0];
rz(-0.94014257) q[0];
rz(-0.6483342) q[1];
sx q[1];
rz(-1.9377361) q[1];
sx q[1];
rz(0.26888332) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20737442) q[0];
sx q[0];
rz(-1.8158127) q[0];
sx q[0];
rz(-0.55565683) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2598632) q[2];
sx q[2];
rz(-2.249392) q[2];
sx q[2];
rz(0.95695615) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2317048) q[1];
sx q[1];
rz(-0.99867601) q[1];
sx q[1];
rz(2.8673299) q[1];
rz(-pi) q[2];
rz(2.4869789) q[3];
sx q[3];
rz(-2.6179492) q[3];
sx q[3];
rz(-2.3414827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.0073283422) q[2];
sx q[2];
rz(-1.2790044) q[2];
sx q[2];
rz(-0.68592611) q[2];
rz(-0.736233) q[3];
sx q[3];
rz(-2.1015621) q[3];
sx q[3];
rz(0.22549103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69970423) q[0];
sx q[0];
rz(-1.2208953) q[0];
sx q[0];
rz(-2.3928483) q[0];
rz(3.0486095) q[1];
sx q[1];
rz(-0.75983202) q[1];
sx q[1];
rz(2.0704796) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7629246) q[0];
sx q[0];
rz(-2.0708186) q[0];
sx q[0];
rz(1.0944233) q[0];
x q[1];
rz(0.57141177) q[2];
sx q[2];
rz(-2.4021752) q[2];
sx q[2];
rz(2.6020056) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.4680571) q[1];
sx q[1];
rz(-0.92294611) q[1];
sx q[1];
rz(0.38783973) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0893965) q[3];
sx q[3];
rz(-1.1602656) q[3];
sx q[3];
rz(-3.0018788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45212713) q[2];
sx q[2];
rz(-2.5048544) q[2];
sx q[2];
rz(0.20991906) q[2];
rz(0.12957761) q[3];
sx q[3];
rz(-2.1942873) q[3];
sx q[3];
rz(-1.5773704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0196446) q[0];
sx q[0];
rz(-0.26109281) q[0];
sx q[0];
rz(0.013539465) q[0];
rz(0.53567046) q[1];
sx q[1];
rz(-0.56709138) q[1];
sx q[1];
rz(-3.1279235) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1197311) q[0];
sx q[0];
rz(-0.80799229) q[0];
sx q[0];
rz(0.31532212) q[0];
rz(-2.057836) q[2];
sx q[2];
rz(-1.9197316) q[2];
sx q[2];
rz(1.9395804) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0108384) q[1];
sx q[1];
rz(-1.2474244) q[1];
sx q[1];
rz(0.73079988) q[1];
x q[2];
rz(-1.1982273) q[3];
sx q[3];
rz(-2.2856224) q[3];
sx q[3];
rz(1.4095588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.51836625) q[2];
sx q[2];
rz(-1.7609111) q[2];
sx q[2];
rz(-1.7784485) q[2];
rz(2.3220883) q[3];
sx q[3];
rz(-0.47911152) q[3];
sx q[3];
rz(0.68377408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0226456) q[0];
sx q[0];
rz(-0.85856694) q[0];
sx q[0];
rz(1.6534506) q[0];
rz(0.34291521) q[1];
sx q[1];
rz(-1.5838793) q[1];
sx q[1];
rz(0.75103474) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2011178) q[0];
sx q[0];
rz(-1.507326) q[0];
sx q[0];
rz(-1.8497784) q[0];
rz(-1.9575167) q[2];
sx q[2];
rz(-2.6474075) q[2];
sx q[2];
rz(1.7750211) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4655059) q[1];
sx q[1];
rz(-1.4526396) q[1];
sx q[1];
rz(-2.5398269) q[1];
x q[2];
rz(3.1084177) q[3];
sx q[3];
rz(-1.4806595) q[3];
sx q[3];
rz(1.7189004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.670383) q[2];
sx q[2];
rz(-2.1977916) q[2];
sx q[2];
rz(0.38718265) q[2];
rz(-2.6662628) q[3];
sx q[3];
rz(-1.9742222) q[3];
sx q[3];
rz(-2.3502684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8925856) q[0];
sx q[0];
rz(-0.96207608) q[0];
sx q[0];
rz(0.67208653) q[0];
rz(-2.7723906) q[1];
sx q[1];
rz(-0.75405706) q[1];
sx q[1];
rz(-1.3934607) q[1];
rz(0.71699981) q[2];
sx q[2];
rz(-1.7115657) q[2];
sx q[2];
rz(-1.905237) q[2];
rz(-3.0283916) q[3];
sx q[3];
rz(-1.8424215) q[3];
sx q[3];
rz(2.7383534) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
