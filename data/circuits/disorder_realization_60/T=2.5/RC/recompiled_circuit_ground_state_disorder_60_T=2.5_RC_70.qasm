OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.085091703) q[0];
sx q[0];
rz(-0.29274517) q[0];
sx q[0];
rz(2.2226287) q[0];
rz(2.7594944) q[1];
sx q[1];
rz(-0.16799071) q[1];
sx q[1];
rz(2.0007432) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98812449) q[0];
sx q[0];
rz(-1.690295) q[0];
sx q[0];
rz(1.7480127) q[0];
rz(-2.8260255) q[2];
sx q[2];
rz(-1.3006388) q[2];
sx q[2];
rz(2.9569654) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4881058) q[1];
sx q[1];
rz(-1.4380903) q[1];
sx q[1];
rz(0.38487558) q[1];
rz(2.010072) q[3];
sx q[3];
rz(-1.695493) q[3];
sx q[3];
rz(1.6623868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4484278) q[2];
sx q[2];
rz(-1.0591256) q[2];
sx q[2];
rz(1.9654174) q[2];
rz(0.042393427) q[3];
sx q[3];
rz(-1.8040801) q[3];
sx q[3];
rz(-0.5347518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6723044) q[0];
sx q[0];
rz(-0.35457087) q[0];
sx q[0];
rz(0.18445036) q[0];
rz(3.0448659) q[1];
sx q[1];
rz(-1.5238785) q[1];
sx q[1];
rz(-1.2299889) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7589129) q[0];
sx q[0];
rz(-1.935674) q[0];
sx q[0];
rz(1.7880102) q[0];
x q[1];
rz(-3.0646992) q[2];
sx q[2];
rz(-0.56083365) q[2];
sx q[2];
rz(-0.70181134) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.36791554) q[1];
sx q[1];
rz(-1.333548) q[1];
sx q[1];
rz(1.3648454) q[1];
rz(-pi) q[2];
x q[2];
rz(1.473068) q[3];
sx q[3];
rz(-0.37993452) q[3];
sx q[3];
rz(-1.987628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.35473287) q[2];
sx q[2];
rz(-2.730098) q[2];
sx q[2];
rz(0.011431781) q[2];
rz(-1.3139906) q[3];
sx q[3];
rz(-1.7766137) q[3];
sx q[3];
rz(2.438681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3979724) q[0];
sx q[0];
rz(-3.008606) q[0];
sx q[0];
rz(2.5308894) q[0];
rz(-0.01344219) q[1];
sx q[1];
rz(-1.2862658) q[1];
sx q[1];
rz(0.59649831) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4134408) q[0];
sx q[0];
rz(-0.96609173) q[0];
sx q[0];
rz(2.938943) q[0];
rz(-pi) q[1];
rz(0.30412401) q[2];
sx q[2];
rz(-2.1275188) q[2];
sx q[2];
rz(0.010526882) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.96287476) q[1];
sx q[1];
rz(-1.8824008) q[1];
sx q[1];
rz(-2.3942663) q[1];
x q[2];
rz(1.6241118) q[3];
sx q[3];
rz(-2.5870596) q[3];
sx q[3];
rz(-2.7729101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.67768031) q[2];
sx q[2];
rz(-1.8296506) q[2];
sx q[2];
rz(-0.00027351969) q[2];
rz(-1.6655946) q[3];
sx q[3];
rz(-0.6690343) q[3];
sx q[3];
rz(3.0345548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45604712) q[0];
sx q[0];
rz(-0.16655971) q[0];
sx q[0];
rz(1.1628994) q[0];
rz(-0.97745013) q[1];
sx q[1];
rz(-1.6012499) q[1];
sx q[1];
rz(-1.5553442) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.453131) q[0];
sx q[0];
rz(-1.5785839) q[0];
sx q[0];
rz(-1.5611737) q[0];
x q[1];
rz(-1.6155717) q[2];
sx q[2];
rz(-0.78592671) q[2];
sx q[2];
rz(-1.8375374) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6628389) q[1];
sx q[1];
rz(-0.49703076) q[1];
sx q[1];
rz(3.1113487) q[1];
x q[2];
rz(-1.076872) q[3];
sx q[3];
rz(-1.3528429) q[3];
sx q[3];
rz(-0.92011425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1996475) q[2];
sx q[2];
rz(-0.5117828) q[2];
sx q[2];
rz(-0.23435782) q[2];
rz(2.6394081) q[3];
sx q[3];
rz(-1.0415223) q[3];
sx q[3];
rz(-2.485399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31768826) q[0];
sx q[0];
rz(-2.3264139) q[0];
sx q[0];
rz(-1.748388) q[0];
rz(0.69270095) q[1];
sx q[1];
rz(-1.0857948) q[1];
sx q[1];
rz(1.0903953) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4596953) q[0];
sx q[0];
rz(-2.5715264) q[0];
sx q[0];
rz(1.9697813) q[0];
rz(-0.88700358) q[2];
sx q[2];
rz(-0.74216026) q[2];
sx q[2];
rz(0.61486828) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.322256) q[1];
sx q[1];
rz(-2.6988479) q[1];
sx q[1];
rz(-1.2135452) q[1];
rz(-pi) q[2];
rz(1.3133009) q[3];
sx q[3];
rz(-2.5867043) q[3];
sx q[3];
rz(-0.99770791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.580487) q[2];
sx q[2];
rz(-1.1539536) q[2];
sx q[2];
rz(0.53606501) q[2];
rz(0.29340336) q[3];
sx q[3];
rz(-1.2111726) q[3];
sx q[3];
rz(2.6611879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(0.90477657) q[0];
sx q[0];
rz(-0.95229709) q[0];
sx q[0];
rz(0.099040898) q[0];
rz(-1.0320484) q[1];
sx q[1];
rz(-2.2910304) q[1];
sx q[1];
rz(1.7031857) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.563614) q[0];
sx q[0];
rz(-0.7398842) q[0];
sx q[0];
rz(-1.9533964) q[0];
rz(-pi) q[1];
rz(2.8929404) q[2];
sx q[2];
rz(-0.2919582) q[2];
sx q[2];
rz(1.1077293) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.21675303) q[1];
sx q[1];
rz(-2.3509563) q[1];
sx q[1];
rz(-0.032921493) q[1];
x q[2];
rz(-2.8436207) q[3];
sx q[3];
rz(-1.4557585) q[3];
sx q[3];
rz(-1.486766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9859887) q[2];
sx q[2];
rz(-0.5039379) q[2];
sx q[2];
rz(-2.3206553) q[2];
rz(1.4486897) q[3];
sx q[3];
rz(-2.4235348) q[3];
sx q[3];
rz(-1.6528355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2445225) q[0];
sx q[0];
rz(-0.91461602) q[0];
sx q[0];
rz(0.50265092) q[0];
rz(-1.2333168) q[1];
sx q[1];
rz(-2.2063467) q[1];
sx q[1];
rz(-0.9800235) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9217958) q[0];
sx q[0];
rz(-1.3338519) q[0];
sx q[0];
rz(-1.8264997) q[0];
rz(-0.77620164) q[2];
sx q[2];
rz(-1.3598249) q[2];
sx q[2];
rz(2.7609776) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0539541) q[1];
sx q[1];
rz(-0.77869895) q[1];
sx q[1];
rz(-0.33233541) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84021826) q[3];
sx q[3];
rz(-1.7257236) q[3];
sx q[3];
rz(-1.8904101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22214733) q[2];
sx q[2];
rz(-1.28747) q[2];
sx q[2];
rz(1.7670828) q[2];
rz(1.8253271) q[3];
sx q[3];
rz(-0.85597435) q[3];
sx q[3];
rz(-2.159806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027243622) q[0];
sx q[0];
rz(-2.7482432) q[0];
sx q[0];
rz(-0.91841665) q[0];
rz(2.9367327) q[1];
sx q[1];
rz(-1.058895) q[1];
sx q[1];
rz(-1.6580261) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0201538) q[0];
sx q[0];
rz(-0.32628548) q[0];
sx q[0];
rz(-2.043528) q[0];
rz(3.0181418) q[2];
sx q[2];
rz(-2.3268893) q[2];
sx q[2];
rz(-1.3208117) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5402447) q[1];
sx q[1];
rz(-1.6588604) q[1];
sx q[1];
rz(-0.87012) q[1];
rz(1.5213883) q[3];
sx q[3];
rz(-2.2893583) q[3];
sx q[3];
rz(2.3895532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.001699) q[2];
sx q[2];
rz(-2.6504982) q[2];
sx q[2];
rz(1.3341382) q[2];
rz(0.88207465) q[3];
sx q[3];
rz(-0.69552723) q[3];
sx q[3];
rz(1.849256) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0598711) q[0];
sx q[0];
rz(-0.43161714) q[0];
sx q[0];
rz(3.1234142) q[0];
rz(-0.40953088) q[1];
sx q[1];
rz(-1.4298226) q[1];
sx q[1];
rz(-0.10083625) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5715953) q[0];
sx q[0];
rz(-0.93678189) q[0];
sx q[0];
rz(1.5279395) q[0];
rz(-0.11282044) q[2];
sx q[2];
rz(-0.42306468) q[2];
sx q[2];
rz(-3.0213838) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.455115) q[1];
sx q[1];
rz(-1.2077189) q[1];
sx q[1];
rz(2.8836714) q[1];
rz(-pi) q[2];
rz(-2.0193897) q[3];
sx q[3];
rz(-0.27315419) q[3];
sx q[3];
rz(0.34452439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8044901) q[2];
sx q[2];
rz(-0.10919315) q[2];
sx q[2];
rz(1.8936554) q[2];
rz(-1.2248056) q[3];
sx q[3];
rz(-2.2962511) q[3];
sx q[3];
rz(0.065751806) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9058022) q[0];
sx q[0];
rz(-1.4757272) q[0];
sx q[0];
rz(-0.32050785) q[0];
rz(0.39101741) q[1];
sx q[1];
rz(-1.1415488) q[1];
sx q[1];
rz(-1.5919707) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029767903) q[0];
sx q[0];
rz(-0.27720189) q[0];
sx q[0];
rz(-0.76407822) q[0];
rz(-0.51410772) q[2];
sx q[2];
rz(-1.7127345) q[2];
sx q[2];
rz(0.98279233) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84234778) q[1];
sx q[1];
rz(-1.4834373) q[1];
sx q[1];
rz(-1.2016988) q[1];
x q[2];
rz(-1.5783089) q[3];
sx q[3];
rz(-1.148461) q[3];
sx q[3];
rz(-2.0762805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0512507) q[2];
sx q[2];
rz(-2.0917459) q[2];
sx q[2];
rz(2.7412097) q[2];
rz(-2.9054437) q[3];
sx q[3];
rz(-0.103424) q[3];
sx q[3];
rz(2.1307438) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.878933) q[0];
sx q[0];
rz(-2.6317609) q[0];
sx q[0];
rz(1.2608933) q[0];
rz(0.48514584) q[1];
sx q[1];
rz(-2.3427675) q[1];
sx q[1];
rz(2.6642703) q[1];
rz(0.12564364) q[2];
sx q[2];
rz(-1.3640654) q[2];
sx q[2];
rz(-2.9403093) q[2];
rz(1.8386091) q[3];
sx q[3];
rz(-0.23614863) q[3];
sx q[3];
rz(1.9105259) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
