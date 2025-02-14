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
rz(-1.5771447) q[0];
sx q[0];
rz(0.91140437) q[0];
rz(2.4721594) q[1];
sx q[1];
rz(-2.8652006) q[1];
sx q[1];
rz(-0.82958329) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54691168) q[0];
sx q[0];
rz(-1.24238) q[0];
sx q[0];
rz(-1.9921147) q[0];
rz(-0.24721036) q[2];
sx q[2];
rz(-1.6610378) q[2];
sx q[2];
rz(-1.5238786) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9467643) q[1];
sx q[1];
rz(-2.2948537) q[1];
sx q[1];
rz(1.1652693) q[1];
rz(-pi) q[2];
rz(-0.31949701) q[3];
sx q[3];
rz(-0.45630672) q[3];
sx q[3];
rz(-0.97684492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.36898819) q[2];
sx q[2];
rz(-0.96869865) q[2];
sx q[2];
rz(1.1645092) q[2];
rz(-2.4146967) q[3];
sx q[3];
rz(-1.3317069) q[3];
sx q[3];
rz(-3.0708142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0182536) q[0];
sx q[0];
rz(-0.66272074) q[0];
sx q[0];
rz(0.29139274) q[0];
rz(0.60028752) q[1];
sx q[1];
rz(-0.95708668) q[1];
sx q[1];
rz(-0.16500638) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5801738) q[0];
sx q[0];
rz(-2.8695171) q[0];
sx q[0];
rz(1.3035357) q[0];
rz(-pi) q[1];
rz(2.0197312) q[2];
sx q[2];
rz(-1.4693207) q[2];
sx q[2];
rz(1.6368653) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.2272039) q[1];
sx q[1];
rz(-1.5014577) q[1];
sx q[1];
rz(0.68080618) q[1];
rz(-pi) q[2];
rz(-0.24834541) q[3];
sx q[3];
rz(-1.5178866) q[3];
sx q[3];
rz(-1.9236881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0773641) q[2];
sx q[2];
rz(-0.95244971) q[2];
sx q[2];
rz(0.23951086) q[2];
rz(0.32660487) q[3];
sx q[3];
rz(-0.17290792) q[3];
sx q[3];
rz(-0.094836205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8851682) q[0];
sx q[0];
rz(-0.435985) q[0];
sx q[0];
rz(-1.5727795) q[0];
rz(2.7478711) q[1];
sx q[1];
rz(-0.55501333) q[1];
sx q[1];
rz(0.47600019) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25903738) q[0];
sx q[0];
rz(-2.5406484) q[0];
sx q[0];
rz(-1.9447692) q[0];
rz(-pi) q[1];
rz(-0.030124859) q[2];
sx q[2];
rz(-1.3839514) q[2];
sx q[2];
rz(3.0812936) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2326374) q[1];
sx q[1];
rz(-0.74027762) q[1];
sx q[1];
rz(-0.93772447) q[1];
rz(-pi) q[2];
rz(-2.1952403) q[3];
sx q[3];
rz(-2.7331684) q[3];
sx q[3];
rz(-0.32865903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8126882) q[2];
sx q[2];
rz(-2.4220971) q[2];
sx q[2];
rz(-2.1997814) q[2];
rz(1.719126) q[3];
sx q[3];
rz(-0.37426451) q[3];
sx q[3];
rz(0.67371887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0734237) q[0];
sx q[0];
rz(-0.24003679) q[0];
sx q[0];
rz(-1.5330676) q[0];
rz(0.34670058) q[1];
sx q[1];
rz(-1.0418714) q[1];
sx q[1];
rz(-2.9516721) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1153206) q[0];
sx q[0];
rz(-0.85912017) q[0];
sx q[0];
rz(1.222696) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0713351) q[2];
sx q[2];
rz(-0.53218666) q[2];
sx q[2];
rz(2.4376873) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7334711) q[1];
sx q[1];
rz(-1.6279393) q[1];
sx q[1];
rz(0.96459324) q[1];
x q[2];
rz(0.40855405) q[3];
sx q[3];
rz(-1.9667278) q[3];
sx q[3];
rz(0.2847288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7803663) q[2];
sx q[2];
rz(-2.2915338) q[2];
sx q[2];
rz(-2.7596607) q[2];
rz(0.10968883) q[3];
sx q[3];
rz(-2.1083125) q[3];
sx q[3];
rz(-3.0221353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1031951) q[0];
sx q[0];
rz(-1.1829475) q[0];
sx q[0];
rz(0.84684816) q[0];
rz(-3.0021744) q[1];
sx q[1];
rz(-0.81652313) q[1];
sx q[1];
rz(-2.3925508) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7199166) q[0];
sx q[0];
rz(-2.4329209) q[0];
sx q[0];
rz(-2.9139374) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48495082) q[2];
sx q[2];
rz(-1.4451988) q[2];
sx q[2];
rz(1.6750977) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74816675) q[1];
sx q[1];
rz(-1.4438549) q[1];
sx q[1];
rz(-0.19537651) q[1];
rz(-pi) q[2];
rz(-1.7827634) q[3];
sx q[3];
rz(-2.6304768) q[3];
sx q[3];
rz(-1.2455452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2060812) q[2];
sx q[2];
rz(-1.8847909) q[2];
sx q[2];
rz(-2.0743267) q[2];
rz(0.71077985) q[3];
sx q[3];
rz(-1.658541) q[3];
sx q[3];
rz(-1.7162286) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34273219) q[0];
sx q[0];
rz(-1.4787759) q[0];
sx q[0];
rz(-0.96858281) q[0];
rz(2.4505278) q[1];
sx q[1];
rz(-1.3422809) q[1];
sx q[1];
rz(-1.8341281) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1145612) q[0];
sx q[0];
rz(-1.8252715) q[0];
sx q[0];
rz(-0.64395321) q[0];
rz(1.3771626) q[2];
sx q[2];
rz(-1.2848952) q[2];
sx q[2];
rz(-1.658197) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5777052) q[1];
sx q[1];
rz(-1.1839377) q[1];
sx q[1];
rz(-3.0847286) q[1];
rz(-pi) q[2];
rz(-2.9581466) q[3];
sx q[3];
rz(-1.3058834) q[3];
sx q[3];
rz(1.2353171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.4751733) q[2];
sx q[2];
rz(-2.9066777) q[2];
sx q[2];
rz(-1.8443745) q[2];
rz(1.3951067) q[3];
sx q[3];
rz(-1.3336983) q[3];
sx q[3];
rz(1.5313799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.7995826) q[0];
sx q[0];
rz(-0.81556773) q[0];
sx q[0];
rz(0.94014257) q[0];
rz(-0.6483342) q[1];
sx q[1];
rz(-1.2038566) q[1];
sx q[1];
rz(2.8727093) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4057345) q[0];
sx q[0];
rz(-0.60204217) q[0];
sx q[0];
rz(-2.6989537) q[0];
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
rz(-1.9098879) q[1];
sx q[1];
rz(-2.1429166) q[1];
sx q[1];
rz(0.27426274) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.712065) q[3];
sx q[3];
rz(-1.2614354) q[3];
sx q[3];
rz(-1.3573028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.0073283422) q[2];
sx q[2];
rz(-1.8625883) q[2];
sx q[2];
rz(-0.68592611) q[2];
rz(2.4053597) q[3];
sx q[3];
rz(-1.0400306) q[3];
sx q[3];
rz(2.9161016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69970423) q[0];
sx q[0];
rz(-1.9206973) q[0];
sx q[0];
rz(2.3928483) q[0];
rz(3.0486095) q[1];
sx q[1];
rz(-2.3817606) q[1];
sx q[1];
rz(1.071113) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37866805) q[0];
sx q[0];
rz(-1.0707741) q[0];
sx q[0];
rz(-1.0944233) q[0];
rz(-pi) q[1];
rz(2.0290211) q[2];
sx q[2];
rz(-2.1734218) q[2];
sx q[2];
rz(-0.17652179) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.4680571) q[1];
sx q[1];
rz(-2.2186465) q[1];
sx q[1];
rz(0.38783973) q[1];
x q[2];
rz(-1.6900916) q[3];
sx q[3];
rz(-0.41364851) q[3];
sx q[3];
rz(-3.1320437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.45212713) q[2];
sx q[2];
rz(-0.63673821) q[2];
sx q[2];
rz(0.20991906) q[2];
rz(-3.012015) q[3];
sx q[3];
rz(-2.1942873) q[3];
sx q[3];
rz(-1.5773704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.121948) q[0];
sx q[0];
rz(-0.26109281) q[0];
sx q[0];
rz(0.013539465) q[0];
rz(2.6059222) q[1];
sx q[1];
rz(-2.5745013) q[1];
sx q[1];
rz(0.013669107) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5608198) q[0];
sx q[0];
rz(-2.3285064) q[0];
sx q[0];
rz(-1.2570501) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91000189) q[2];
sx q[2];
rz(-2.5507413) q[2];
sx q[2];
rz(0.94205085) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0108384) q[1];
sx q[1];
rz(-1.2474244) q[1];
sx q[1];
rz(-0.73079988) q[1];
rz(-pi) q[2];
rz(-1.1982273) q[3];
sx q[3];
rz(-0.85597023) q[3];
sx q[3];
rz(-1.4095588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.51836625) q[2];
sx q[2];
rz(-1.3806815) q[2];
sx q[2];
rz(-1.3631442) q[2];
rz(-2.3220883) q[3];
sx q[3];
rz(-0.47911152) q[3];
sx q[3];
rz(2.4578186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.118947) q[0];
sx q[0];
rz(-0.85856694) q[0];
sx q[0];
rz(-1.6534506) q[0];
rz(-2.7986774) q[1];
sx q[1];
rz(-1.5838793) q[1];
sx q[1];
rz(-2.3905579) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5540051) q[0];
sx q[0];
rz(-2.8556654) q[0];
sx q[0];
rz(1.343973) q[0];
rz(-0.20047156) q[2];
sx q[2];
rz(-1.1159917) q[2];
sx q[2];
rz(-1.7998296) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.955525) q[1];
sx q[1];
rz(-2.1677818) q[1];
sx q[1];
rz(1.7138193) q[1];
rz(-pi) q[2];
rz(-1.6609825) q[3];
sx q[3];
rz(-1.5377561) q[3];
sx q[3];
rz(-0.15109135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.670383) q[2];
sx q[2];
rz(-2.1977916) q[2];
sx q[2];
rz(0.38718265) q[2];
rz(-0.47532982) q[3];
sx q[3];
rz(-1.9742222) q[3];
sx q[3];
rz(2.3502684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8925856) q[0];
sx q[0];
rz(-2.1795166) q[0];
sx q[0];
rz(-2.4695061) q[0];
rz(-0.36920209) q[1];
sx q[1];
rz(-2.3875356) q[1];
sx q[1];
rz(1.748132) q[1];
rz(1.7566219) q[2];
sx q[2];
rz(-0.86238774) q[2];
sx q[2];
rz(2.9288616) q[2];
rz(-1.297507) q[3];
sx q[3];
rz(-1.4617625) q[3];
sx q[3];
rz(1.1980496) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
