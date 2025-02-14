OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.1104133) q[0];
sx q[0];
rz(-2.2194982) q[0];
sx q[0];
rz(-2.3715012) q[0];
rz(0.71647477) q[1];
sx q[1];
rz(-2.1991576) q[1];
sx q[1];
rz(0.64185774) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3549558) q[0];
sx q[0];
rz(-1.6842951) q[0];
sx q[0];
rz(2.9391975) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0014512295) q[2];
sx q[2];
rz(-1.5719169) q[2];
sx q[2];
rz(-1.4933153) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9145467) q[1];
sx q[1];
rz(-1.7265336) q[1];
sx q[1];
rz(2.7776633) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3083616) q[3];
sx q[3];
rz(-2.0249484) q[3];
sx q[3];
rz(-1.7851225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.98051071) q[2];
sx q[2];
rz(-1.0970205) q[2];
sx q[2];
rz(0.80992997) q[2];
rz(3.1071281) q[3];
sx q[3];
rz(-2.4813215) q[3];
sx q[3];
rz(3.0380429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.506839) q[0];
sx q[0];
rz(-2.6320808) q[0];
sx q[0];
rz(2.4822045) q[0];
rz(1.4393282) q[1];
sx q[1];
rz(-1.5123475) q[1];
sx q[1];
rz(2.6606681) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.07671) q[0];
sx q[0];
rz(-1.1336898) q[0];
sx q[0];
rz(-2.9124898) q[0];
rz(-pi) q[1];
rz(-0.73591097) q[2];
sx q[2];
rz(-0.99756587) q[2];
sx q[2];
rz(0.57511759) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9240504) q[1];
sx q[1];
rz(-0.21044193) q[1];
sx q[1];
rz(-1.7173052) q[1];
rz(-pi) q[2];
rz(2.0388076) q[3];
sx q[3];
rz(-2.0986522) q[3];
sx q[3];
rz(0.46609391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3769569) q[2];
sx q[2];
rz(-0.83676338) q[2];
sx q[2];
rz(-2.5118206) q[2];
rz(-1.9624286) q[3];
sx q[3];
rz(-2.4383014) q[3];
sx q[3];
rz(2.0298957) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69029194) q[0];
sx q[0];
rz(-0.010882219) q[0];
sx q[0];
rz(-0.96845281) q[0];
rz(-0.14006607) q[1];
sx q[1];
rz(-1.3532956) q[1];
sx q[1];
rz(2.5586939) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.934865) q[0];
sx q[0];
rz(-1.7375653) q[0];
sx q[0];
rz(1.5251446) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3084564) q[2];
sx q[2];
rz(-2.3916187) q[2];
sx q[2];
rz(1.394608) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1212226) q[1];
sx q[1];
rz(-2.2285151) q[1];
sx q[1];
rz(-2.4037564) q[1];
rz(1.8898592) q[3];
sx q[3];
rz(-1.3033861) q[3];
sx q[3];
rz(1.9134932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6110903) q[2];
sx q[2];
rz(-3.0641596) q[2];
sx q[2];
rz(2.9275242) q[2];
rz(-2.8392082) q[3];
sx q[3];
rz(-2.3979135) q[3];
sx q[3];
rz(0.42588699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9532303) q[0];
sx q[0];
rz(-1.0091877) q[0];
sx q[0];
rz(2.0971712) q[0];
rz(2.2638679) q[1];
sx q[1];
rz(-1.5526155) q[1];
sx q[1];
rz(-2.9728319) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2450352) q[0];
sx q[0];
rz(-1.3447133) q[0];
sx q[0];
rz(-1.0624136) q[0];
x q[1];
rz(1.020476) q[2];
sx q[2];
rz(-0.3786217) q[2];
sx q[2];
rz(1.0892717) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8873621) q[1];
sx q[1];
rz(-1.7866644) q[1];
sx q[1];
rz(2.5309023) q[1];
x q[2];
rz(2.2764858) q[3];
sx q[3];
rz(-1.5388515) q[3];
sx q[3];
rz(2.3695996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3622482) q[2];
sx q[2];
rz(-1.8579973) q[2];
sx q[2];
rz(-1.0603504) q[2];
rz(2.849071) q[3];
sx q[3];
rz(-2.4157603) q[3];
sx q[3];
rz(3.0574851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5523858) q[0];
sx q[0];
rz(-0.64738208) q[0];
sx q[0];
rz(-0.86828434) q[0];
rz(2.5153416) q[1];
sx q[1];
rz(-1.3887082) q[1];
sx q[1];
rz(-1.7832322) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6407961) q[0];
sx q[0];
rz(-1.8587041) q[0];
sx q[0];
rz(-1.4782216) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3542074) q[2];
sx q[2];
rz(-1.1973477) q[2];
sx q[2];
rz(0.83080705) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.89857453) q[1];
sx q[1];
rz(-2.822091) q[1];
sx q[1];
rz(-2.363149) q[1];
rz(-pi) q[2];
rz(0.56760889) q[3];
sx q[3];
rz(-2.4205812) q[3];
sx q[3];
rz(0.22338671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2061578) q[2];
sx q[2];
rz(-0.81659603) q[2];
sx q[2];
rz(-0.52331501) q[2];
rz(0.37007904) q[3];
sx q[3];
rz(-0.78032929) q[3];
sx q[3];
rz(-1.3453329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0321781) q[0];
sx q[0];
rz(-2.9258756) q[0];
sx q[0];
rz(2.7874462) q[0];
rz(-0.94193637) q[1];
sx q[1];
rz(-1.4553921) q[1];
sx q[1];
rz(-1.8249493) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2433853) q[0];
sx q[0];
rz(-1.6750209) q[0];
sx q[0];
rz(-2.6263155) q[0];
rz(-pi) q[1];
rz(-2.5923827) q[2];
sx q[2];
rz(-1.5592279) q[2];
sx q[2];
rz(2.9988097) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5806953) q[1];
sx q[1];
rz(-1.744701) q[1];
sx q[1];
rz(-1.6328801) q[1];
rz(-pi) q[2];
rz(-1.8527669) q[3];
sx q[3];
rz(-0.94579711) q[3];
sx q[3];
rz(-0.19240141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3256623) q[2];
sx q[2];
rz(-2.7553835) q[2];
sx q[2];
rz(0.33829921) q[2];
rz(0.47652388) q[3];
sx q[3];
rz(-2.3881113) q[3];
sx q[3];
rz(-2.8009955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4395831) q[0];
sx q[0];
rz(-1.5583353) q[0];
sx q[0];
rz(2.6038792) q[0];
rz(-1.733755) q[1];
sx q[1];
rz(-0.45999637) q[1];
sx q[1];
rz(2.5163311) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94948927) q[0];
sx q[0];
rz(-1.2117627) q[0];
sx q[0];
rz(2.5395995) q[0];
rz(-pi) q[1];
x q[1];
rz(2.602732) q[2];
sx q[2];
rz(-1.5075353) q[2];
sx q[2];
rz(0.45743313) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.75692528) q[1];
sx q[1];
rz(-2.1672492) q[1];
sx q[1];
rz(1.8310962) q[1];
x q[2];
rz(-1.2042562) q[3];
sx q[3];
rz(-0.82895422) q[3];
sx q[3];
rz(-0.77223611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9739146) q[2];
sx q[2];
rz(-2.6748071) q[2];
sx q[2];
rz(-1.7612877) q[2];
rz(-0.4365094) q[3];
sx q[3];
rz(-1.0218388) q[3];
sx q[3];
rz(0.71353394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9647144) q[0];
sx q[0];
rz(-2.5202993) q[0];
sx q[0];
rz(-0.51625133) q[0];
rz(-0.77975726) q[1];
sx q[1];
rz(-0.98194352) q[1];
sx q[1];
rz(1.0468743) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6927495) q[0];
sx q[0];
rz(-1.7147281) q[0];
sx q[0];
rz(-0.30152614) q[0];
x q[1];
rz(0.5811695) q[2];
sx q[2];
rz(-1.1343252) q[2];
sx q[2];
rz(-1.798686) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6987933) q[1];
sx q[1];
rz(-1.6644786) q[1];
sx q[1];
rz(-3.0832861) q[1];
rz(0.67892142) q[3];
sx q[3];
rz(-2.553294) q[3];
sx q[3];
rz(-0.029840851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.20389916) q[2];
sx q[2];
rz(-0.1940618) q[2];
sx q[2];
rz(0.95721179) q[2];
rz(2.8474478) q[3];
sx q[3];
rz(-1.1295986) q[3];
sx q[3];
rz(2.8637776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99943632) q[0];
sx q[0];
rz(-0.23040982) q[0];
sx q[0];
rz(2.9545422) q[0];
rz(-2.710178) q[1];
sx q[1];
rz(-0.43771935) q[1];
sx q[1];
rz(-1.4923219) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0566643) q[0];
sx q[0];
rz(-1.5001552) q[0];
sx q[0];
rz(-1.4991374) q[0];
rz(-pi) q[1];
rz(-1.0643105) q[2];
sx q[2];
rz(-0.37831719) q[2];
sx q[2];
rz(1.0418721) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.3304676) q[1];
sx q[1];
rz(-0.44422517) q[1];
sx q[1];
rz(-1.1796239) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1077376) q[3];
sx q[3];
rz(-2.9599422) q[3];
sx q[3];
rz(2.8665115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6734068) q[2];
sx q[2];
rz(-0.91172051) q[2];
sx q[2];
rz(-2.1569596) q[2];
rz(2.6268688) q[3];
sx q[3];
rz(-0.52557164) q[3];
sx q[3];
rz(1.908186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.89727) q[0];
sx q[0];
rz(-1.4800625) q[0];
sx q[0];
rz(2.3829714) q[0];
rz(1.1933391) q[1];
sx q[1];
rz(-1.9494282) q[1];
sx q[1];
rz(1.4512482) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7465736) q[0];
sx q[0];
rz(-1.7749471) q[0];
sx q[0];
rz(1.8042807) q[0];
rz(1.7733712) q[2];
sx q[2];
rz(-1.612101) q[2];
sx q[2];
rz(-1.6172723) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2047233) q[1];
sx q[1];
rz(-1.544049) q[1];
sx q[1];
rz(-2.1863947) q[1];
rz(-pi) q[2];
rz(-2.8192807) q[3];
sx q[3];
rz(-1.1813643) q[3];
sx q[3];
rz(-0.001231391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54185581) q[2];
sx q[2];
rz(-2.2208322) q[2];
sx q[2];
rz(2.3403781) q[2];
rz(-2.2090705) q[3];
sx q[3];
rz(-1.9263809) q[3];
sx q[3];
rz(0.41509375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5634609) q[0];
sx q[0];
rz(-1.7628071) q[0];
sx q[0];
rz(2.5264869) q[0];
rz(-2.8201132) q[1];
sx q[1];
rz(-2.165806) q[1];
sx q[1];
rz(1.4478366) q[1];
rz(0.64340016) q[2];
sx q[2];
rz(-0.20558029) q[2];
sx q[2];
rz(-2.554648) q[2];
rz(-0.23033167) q[3];
sx q[3];
rz(-1.8445704) q[3];
sx q[3];
rz(-0.45382378) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
