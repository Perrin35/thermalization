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
rz(-2.4251179) q[1];
sx q[1];
rz(-0.94243503) q[1];
sx q[1];
rz(2.4997349) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72015136) q[0];
sx q[0];
rz(-0.23166616) q[0];
sx q[0];
rz(2.6257444) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5719169) q[2];
sx q[2];
rz(-1.5722476) q[2];
sx q[2];
rz(3.0641132) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0986276) q[1];
sx q[1];
rz(-2.7471099) q[1];
sx q[1];
rz(-0.41542713) q[1];
x q[2];
rz(-2.3083616) q[3];
sx q[3];
rz(-1.1166443) q[3];
sx q[3];
rz(-1.7851225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1610819) q[2];
sx q[2];
rz(-2.0445721) q[2];
sx q[2];
rz(-2.3316627) q[2];
rz(-0.03446456) q[3];
sx q[3];
rz(-0.66027111) q[3];
sx q[3];
rz(-3.0380429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.506839) q[0];
sx q[0];
rz(-0.50951183) q[0];
sx q[0];
rz(0.65938812) q[0];
rz(1.4393282) q[1];
sx q[1];
rz(-1.5123475) q[1];
sx q[1];
rz(2.6606681) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.07671) q[0];
sx q[0];
rz(-1.1336898) q[0];
sx q[0];
rz(2.9124898) q[0];
rz(0.73591097) q[2];
sx q[2];
rz(-0.99756587) q[2];
sx q[2];
rz(2.5664751) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9316599) q[1];
sx q[1];
rz(-1.5402964) q[1];
sx q[1];
rz(1.7790487) q[1];
x q[2];
rz(-1.102785) q[3];
sx q[3];
rz(-1.0429405) q[3];
sx q[3];
rz(-0.46609391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3769569) q[2];
sx q[2];
rz(-0.83676338) q[2];
sx q[2];
rz(2.5118206) q[2];
rz(1.9624286) q[3];
sx q[3];
rz(-2.4383014) q[3];
sx q[3];
rz(-2.0298957) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69029194) q[0];
sx q[0];
rz(-0.010882219) q[0];
sx q[0];
rz(2.1731398) q[0];
rz(0.14006607) q[1];
sx q[1];
rz(-1.3532956) q[1];
sx q[1];
rz(0.5828988) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.934865) q[0];
sx q[0];
rz(-1.7375653) q[0];
sx q[0];
rz(1.616448) q[0];
x q[1];
rz(1.3084564) q[2];
sx q[2];
rz(-2.3916187) q[2];
sx q[2];
rz(1.7469847) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1212226) q[1];
sx q[1];
rz(-2.2285151) q[1];
sx q[1];
rz(-2.4037564) q[1];
x q[2];
rz(0.85286136) q[3];
sx q[3];
rz(-0.4133458) q[3];
sx q[3];
rz(-1.0173544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6110903) q[2];
sx q[2];
rz(-0.077433057) q[2];
sx q[2];
rz(2.9275242) q[2];
rz(2.8392082) q[3];
sx q[3];
rz(-0.74367911) q[3];
sx q[3];
rz(0.42588699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18836235) q[0];
sx q[0];
rz(-1.0091877) q[0];
sx q[0];
rz(2.0971712) q[0];
rz(0.87772477) q[1];
sx q[1];
rz(-1.5889771) q[1];
sx q[1];
rz(-2.9728319) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2450352) q[0];
sx q[0];
rz(-1.3447133) q[0];
sx q[0];
rz(-1.0624136) q[0];
rz(2.9364768) q[2];
sx q[2];
rz(-1.8913219) q[2];
sx q[2];
rz(-1.4687302) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8873621) q[1];
sx q[1];
rz(-1.7866644) q[1];
sx q[1];
rz(0.61069031) q[1];
rz(-2.2764858) q[3];
sx q[3];
rz(-1.5388515) q[3];
sx q[3];
rz(0.77199304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3622482) q[2];
sx q[2];
rz(-1.8579973) q[2];
sx q[2];
rz(2.0812422) q[2];
rz(-2.849071) q[3];
sx q[3];
rz(-0.7258324) q[3];
sx q[3];
rz(-0.084107548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58920687) q[0];
sx q[0];
rz(-2.4942106) q[0];
sx q[0];
rz(2.2733083) q[0];
rz(2.5153416) q[1];
sx q[1];
rz(-1.3887082) q[1];
sx q[1];
rz(-1.7832322) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0979484) q[0];
sx q[0];
rz(-1.4820423) q[0];
sx q[0];
rz(0.28907816) q[0];
rz(-pi) q[1];
rz(-0.38154885) q[2];
sx q[2];
rz(-1.7722436) q[2];
sx q[2];
rz(0.82009456) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7029801) q[1];
sx q[1];
rz(-1.3452521) q[1];
sx q[1];
rz(1.799052) q[1];
x q[2];
rz(-2.503848) q[3];
sx q[3];
rz(-1.9336091) q[3];
sx q[3];
rz(-0.90076288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9354349) q[2];
sx q[2];
rz(-0.81659603) q[2];
sx q[2];
rz(-0.52331501) q[2];
rz(-2.7715136) q[3];
sx q[3];
rz(-2.3612634) q[3];
sx q[3];
rz(1.3453329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(0.10941457) q[0];
sx q[0];
rz(-0.21571708) q[0];
sx q[0];
rz(-0.35414645) q[0];
rz(-0.94193637) q[1];
sx q[1];
rz(-1.6862005) q[1];
sx q[1];
rz(1.8249493) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2433853) q[0];
sx q[0];
rz(-1.4665717) q[0];
sx q[0];
rz(0.51527713) q[0];
rz(-pi) q[1];
rz(-1.5572335) q[2];
sx q[2];
rz(-1.0216273) q[2];
sx q[2];
rz(-1.4209335) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1209379) q[1];
sx q[1];
rz(-1.6319425) q[1];
sx q[1];
rz(-2.967359) q[1];
x q[2];
rz(-2.497358) q[3];
sx q[3];
rz(-1.3431864) q[3];
sx q[3];
rz(-1.5952974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3256623) q[2];
sx q[2];
rz(-2.7553835) q[2];
sx q[2];
rz(-2.8032934) q[2];
rz(-2.6650688) q[3];
sx q[3];
rz(-2.3881113) q[3];
sx q[3];
rz(-2.8009955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7020096) q[0];
sx q[0];
rz(-1.5832573) q[0];
sx q[0];
rz(-0.53771341) q[0];
rz(-1.733755) q[1];
sx q[1];
rz(-2.6815963) q[1];
sx q[1];
rz(0.62526155) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94948927) q[0];
sx q[0];
rz(-1.92983) q[0];
sx q[0];
rz(2.5395995) q[0];
x q[1];
rz(1.6444667) q[2];
sx q[2];
rz(-1.0331312) q[2];
sx q[2];
rz(-1.0755838) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1*pi/10) q[1];
sx q[1];
rz(-2.4971967) q[1];
sx q[1];
rz(0.36233904) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77615555) q[3];
sx q[3];
rz(-1.8381834) q[3];
sx q[3];
rz(-0.54477967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9739146) q[2];
sx q[2];
rz(-0.46678552) q[2];
sx q[2];
rz(1.3803049) q[2];
rz(-2.7050833) q[3];
sx q[3];
rz(-2.1197539) q[3];
sx q[3];
rz(0.71353394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1768782) q[0];
sx q[0];
rz(-0.62129337) q[0];
sx q[0];
rz(-0.51625133) q[0];
rz(2.3618354) q[1];
sx q[1];
rz(-2.1596491) q[1];
sx q[1];
rz(-1.0468743) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16653331) q[0];
sx q[0];
rz(-1.2724845) q[0];
sx q[0];
rz(1.7214283) q[0];
rz(-2.0798415) q[2];
sx q[2];
rz(-2.0915481) q[2];
sx q[2];
rz(-2.6428509) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0081353) q[1];
sx q[1];
rz(-1.6288469) q[1];
sx q[1];
rz(-1.4769555) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9674928) q[3];
sx q[3];
rz(-1.124212) q[3];
sx q[3];
rz(-2.4012938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.20389916) q[2];
sx q[2];
rz(-0.1940618) q[2];
sx q[2];
rz(0.95721179) q[2];
rz(-0.29414487) q[3];
sx q[3];
rz(-1.1295986) q[3];
sx q[3];
rz(2.8637776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99943632) q[0];
sx q[0];
rz(-2.9111828) q[0];
sx q[0];
rz(-0.18705046) q[0];
rz(-2.710178) q[1];
sx q[1];
rz(-2.7038733) q[1];
sx q[1];
rz(1.4923219) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6607912) q[0];
sx q[0];
rz(-1.4993164) q[0];
sx q[0];
rz(-3.0707703) q[0];
rz(-pi) q[1];
rz(1.0643105) q[2];
sx q[2];
rz(-0.37831719) q[2];
sx q[2];
rz(2.0997206) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.3304676) q[1];
sx q[1];
rz(-0.44422517) q[1];
sx q[1];
rz(1.9619688) q[1];
rz(-pi) q[2];
rz(1.727333) q[3];
sx q[3];
rz(-1.4782584) q[3];
sx q[3];
rz(-2.375556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.46818587) q[2];
sx q[2];
rz(-2.2298721) q[2];
sx q[2];
rz(2.1569596) q[2];
rz(-0.51472384) q[3];
sx q[3];
rz(-2.616021) q[3];
sx q[3];
rz(-1.908186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24432261) q[0];
sx q[0];
rz(-1.4800625) q[0];
sx q[0];
rz(0.75862128) q[0];
rz(1.9482535) q[1];
sx q[1];
rz(-1.1921644) q[1];
sx q[1];
rz(1.4512482) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8817026) q[0];
sx q[0];
rz(-2.8326747) q[0];
sx q[0];
rz(0.84084671) q[0];
rz(0.042165857) q[2];
sx q[2];
rz(-1.3683967) q[2];
sx q[2];
rz(3.1035977) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.38498653) q[1];
sx q[1];
rz(-2.1861417) q[1];
sx q[1];
rz(-0.032757515) q[1];
rz(-pi) q[2];
rz(1.162446) q[3];
sx q[3];
rz(-1.2733885) q[3];
sx q[3];
rz(-1.698146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.54185581) q[2];
sx q[2];
rz(-0.92076045) q[2];
sx q[2];
rz(-2.3403781) q[2];
rz(-2.2090705) q[3];
sx q[3];
rz(-1.2152117) q[3];
sx q[3];
rz(2.7264989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5781317) q[0];
sx q[0];
rz(-1.7628071) q[0];
sx q[0];
rz(2.5264869) q[0];
rz(-2.8201132) q[1];
sx q[1];
rz(-2.165806) q[1];
sx q[1];
rz(1.4478366) q[1];
rz(-1.6952487) q[2];
sx q[2];
rz(-1.4067408) q[2];
sx q[2];
rz(-0.066700145) q[2];
rz(1.8516171) q[3];
sx q[3];
rz(-1.7924037) q[3];
sx q[3];
rz(-2.0879346) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
