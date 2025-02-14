OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31792274) q[0];
sx q[0];
rz(1.6190785) q[0];
sx q[0];
rz(9.4879307) q[0];
rz(-1.3178648) q[1];
sx q[1];
rz(-2.7422819) q[1];
sx q[1];
rz(0.77594405) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5634034) q[0];
sx q[0];
rz(-3.0886123) q[0];
sx q[0];
rz(0.53081675) q[0];
x q[1];
rz(2.6444998) q[2];
sx q[2];
rz(-2.1462206) q[2];
sx q[2];
rz(-2.6445553) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9777955) q[1];
sx q[1];
rz(-1.7479291) q[1];
sx q[1];
rz(0.85834412) q[1];
rz(0.41764655) q[3];
sx q[3];
rz(-0.54134901) q[3];
sx q[3];
rz(-0.82994313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.45106384) q[2];
sx q[2];
rz(-1.7454001) q[2];
sx q[2];
rz(-2.6032791) q[2];
rz(-0.29317835) q[3];
sx q[3];
rz(-1.5436951) q[3];
sx q[3];
rz(-1.8831016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20741589) q[0];
sx q[0];
rz(-0.84686142) q[0];
sx q[0];
rz(0.22879452) q[0];
rz(-3.0711807) q[1];
sx q[1];
rz(-1.2733302) q[1];
sx q[1];
rz(1.061903) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5492229) q[0];
sx q[0];
rz(-1.8350936) q[0];
sx q[0];
rz(-2.6068674) q[0];
x q[1];
rz(-2.0309762) q[2];
sx q[2];
rz(-2.0041582) q[2];
sx q[2];
rz(0.037511911) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6390266) q[1];
sx q[1];
rz(-0.4430534) q[1];
sx q[1];
rz(0.73213099) q[1];
x q[2];
rz(2.6999989) q[3];
sx q[3];
rz(-2.674235) q[3];
sx q[3];
rz(0.87477126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.215302) q[2];
sx q[2];
rz(-2.6965202) q[2];
sx q[2];
rz(3.0253547) q[2];
rz(0.91533533) q[3];
sx q[3];
rz(-1.4696308) q[3];
sx q[3];
rz(2.5168929) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6621458) q[0];
sx q[0];
rz(-1.5597458) q[0];
sx q[0];
rz(0.33552718) q[0];
rz(-1.7299995) q[1];
sx q[1];
rz(-0.84452191) q[1];
sx q[1];
rz(1.1988877) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9379282) q[0];
sx q[0];
rz(-2.0988582) q[0];
sx q[0];
rz(0.33007921) q[0];
x q[1];
rz(-2.6012318) q[2];
sx q[2];
rz(-2.4746404) q[2];
sx q[2];
rz(1.6467384) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4863805) q[1];
sx q[1];
rz(-0.81334121) q[1];
sx q[1];
rz(1.7310623) q[1];
x q[2];
rz(-2.8923404) q[3];
sx q[3];
rz(-1.0947168) q[3];
sx q[3];
rz(0.92633807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.18232839) q[2];
sx q[2];
rz(-0.79572833) q[2];
sx q[2];
rz(-2.0286593) q[2];
rz(0.5111323) q[3];
sx q[3];
rz(-1.3730349) q[3];
sx q[3];
rz(-0.50890499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48753259) q[0];
sx q[0];
rz(-0.075981058) q[0];
sx q[0];
rz(-2.7677166) q[0];
rz(-1.182425) q[1];
sx q[1];
rz(-1.1610718) q[1];
sx q[1];
rz(-2.9808796) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99237011) q[0];
sx q[0];
rz(-2.9806031) q[0];
sx q[0];
rz(-0.54246728) q[0];
x q[1];
rz(2.2437975) q[2];
sx q[2];
rz(-0.89859573) q[2];
sx q[2];
rz(0.77122818) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.070878167) q[1];
sx q[1];
rz(-1.1859815) q[1];
sx q[1];
rz(1.0139731) q[1];
rz(-0.31779685) q[3];
sx q[3];
rz(-1.8565912) q[3];
sx q[3];
rz(2.5366549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.60761991) q[2];
sx q[2];
rz(-2.0325568) q[2];
sx q[2];
rz(2.4268761) q[2];
rz(-2.886582) q[3];
sx q[3];
rz(-1.3304354) q[3];
sx q[3];
rz(-2.3431006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0890927) q[0];
sx q[0];
rz(-2.8369501) q[0];
sx q[0];
rz(0.76328817) q[0];
rz(-2.8870562) q[1];
sx q[1];
rz(-1.3724047) q[1];
sx q[1];
rz(1.4483784) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3445374) q[0];
sx q[0];
rz(-1.2303924) q[0];
sx q[0];
rz(-1.1974687) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15621846) q[2];
sx q[2];
rz(-2.3153981) q[2];
sx q[2];
rz(-0.35075364) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9336945) q[1];
sx q[1];
rz(-0.95732821) q[1];
sx q[1];
rz(-1.3529569) q[1];
rz(1.17286) q[3];
sx q[3];
rz(-0.89943652) q[3];
sx q[3];
rz(-1.9111247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.280507) q[2];
sx q[2];
rz(-2.0528767) q[2];
sx q[2];
rz(-1.0293993) q[2];
rz(1.5329125) q[3];
sx q[3];
rz(-2.2423988) q[3];
sx q[3];
rz(2.5077584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7480046) q[0];
sx q[0];
rz(-0.91359502) q[0];
sx q[0];
rz(2.9314801) q[0];
rz(-0.70136079) q[1];
sx q[1];
rz(-1.7751834) q[1];
sx q[1];
rz(-1.1022386) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1195591) q[0];
sx q[0];
rz(-0.73444542) q[0];
sx q[0];
rz(2.2950315) q[0];
rz(-pi) q[1];
rz(3.1386078) q[2];
sx q[2];
rz(-2.4793803) q[2];
sx q[2];
rz(2.4110297) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.452207) q[1];
sx q[1];
rz(-1.8836792) q[1];
sx q[1];
rz(2.41928) q[1];
rz(-pi) q[2];
rz(-2.5999448) q[3];
sx q[3];
rz(-1.4337501) q[3];
sx q[3];
rz(-0.26462072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8474951) q[2];
sx q[2];
rz(-1.5832486) q[2];
sx q[2];
rz(-1.4964649) q[2];
rz(1.7485113) q[3];
sx q[3];
rz(-2.2025509) q[3];
sx q[3];
rz(0.67162544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.3944655) q[0];
sx q[0];
rz(-1.324993) q[0];
sx q[0];
rz(2.6649244) q[0];
rz(1.7851625) q[1];
sx q[1];
rz(-1.8442267) q[1];
sx q[1];
rz(2.2043601) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.40537) q[0];
sx q[0];
rz(-0.69044603) q[0];
sx q[0];
rz(1.6274981) q[0];
rz(-0.071120485) q[2];
sx q[2];
rz(-2.4223339) q[2];
sx q[2];
rz(-0.90380441) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2768749) q[1];
sx q[1];
rz(-1.9559033) q[1];
sx q[1];
rz(2.6996946) q[1];
rz(-0.44853278) q[3];
sx q[3];
rz(-1.6584907) q[3];
sx q[3];
rz(2.0929071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53200191) q[2];
sx q[2];
rz(-1.3401745) q[2];
sx q[2];
rz(-1.0756005) q[2];
rz(-0.39014751) q[3];
sx q[3];
rz(-1.8307999) q[3];
sx q[3];
rz(2.3384317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53668642) q[0];
sx q[0];
rz(-1.0746047) q[0];
sx q[0];
rz(-2.5965776) q[0];
rz(2.1489428) q[1];
sx q[1];
rz(-2.1558709) q[1];
sx q[1];
rz(-1.4535905) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6777991) q[0];
sx q[0];
rz(-0.92831573) q[0];
sx q[0];
rz(2.5911261) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9729068) q[2];
sx q[2];
rz(-0.14384362) q[2];
sx q[2];
rz(0.84022249) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7632227) q[1];
sx q[1];
rz(-1.6838131) q[1];
sx q[1];
rz(-2.3341353) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5239632) q[3];
sx q[3];
rz(-1.5107656) q[3];
sx q[3];
rz(-1.5200652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40889007) q[2];
sx q[2];
rz(-1.7444892) q[2];
sx q[2];
rz(-2.8453541) q[2];
rz(-0.57684165) q[3];
sx q[3];
rz(-2.7448765) q[3];
sx q[3];
rz(1.2805987) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7883564) q[0];
sx q[0];
rz(-0.78859538) q[0];
sx q[0];
rz(-0.22407918) q[0];
rz(-0.60399404) q[1];
sx q[1];
rz(-0.9915587) q[1];
sx q[1];
rz(1.6557065) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1093134) q[0];
sx q[0];
rz(-1.3242189) q[0];
sx q[0];
rz(0.37301491) q[0];
x q[1];
rz(2.9913556) q[2];
sx q[2];
rz(-1.3753819) q[2];
sx q[2];
rz(-3.0216211) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5331577) q[1];
sx q[1];
rz(-0.67109334) q[1];
sx q[1];
rz(-0.41678269) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6194686) q[3];
sx q[3];
rz(-2.587954) q[3];
sx q[3];
rz(1.6185135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.47933444) q[2];
sx q[2];
rz(-1.0666288) q[2];
sx q[2];
rz(0.38702854) q[2];
rz(0.28338638) q[3];
sx q[3];
rz(-2.3895388) q[3];
sx q[3];
rz(-2.7391105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0249483) q[0];
sx q[0];
rz(-2.8706757) q[0];
sx q[0];
rz(1.1234294) q[0];
rz(-1.4943538) q[1];
sx q[1];
rz(-1.6554183) q[1];
sx q[1];
rz(2.2656238) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4433243) q[0];
sx q[0];
rz(-2.351555) q[0];
sx q[0];
rz(-0.36418551) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89089762) q[2];
sx q[2];
rz(-0.65717893) q[2];
sx q[2];
rz(0.79492043) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3265423) q[1];
sx q[1];
rz(-2.2382394) q[1];
sx q[1];
rz(-2.8615281) q[1];
x q[2];
rz(-0.42966162) q[3];
sx q[3];
rz(-1.7335827) q[3];
sx q[3];
rz(1.2186039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4103849) q[2];
sx q[2];
rz(-1.7816252) q[2];
sx q[2];
rz(-3.0044921) q[2];
rz(-1.7588663) q[3];
sx q[3];
rz(-0.9226678) q[3];
sx q[3];
rz(-1.2285129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.726534) q[0];
sx q[0];
rz(-2.6714323) q[0];
sx q[0];
rz(2.5191125) q[0];
rz(2.7525735) q[1];
sx q[1];
rz(-0.35094378) q[1];
sx q[1];
rz(2.6019179) q[1];
rz(-0.48609235) q[2];
sx q[2];
rz(-1.9584502) q[2];
sx q[2];
rz(2.9050585) q[2];
rz(-0.20122726) q[3];
sx q[3];
rz(-1.22898) q[3];
sx q[3];
rz(-2.1635319) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
