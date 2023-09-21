OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0302304) q[0];
sx q[0];
rz(-0.65519873) q[0];
sx q[0];
rz(0.97638786) q[0];
rz(2.2250277) q[1];
sx q[1];
rz(-2.38148) q[1];
sx q[1];
rz(-2.3488933) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8108609) q[0];
sx q[0];
rz(-0.90002492) q[0];
sx q[0];
rz(0.030681507) q[0];
rz(-1.3408469) q[2];
sx q[2];
rz(-1.0385498) q[2];
sx q[2];
rz(1.4872273) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.77289171) q[1];
sx q[1];
rz(-1.4979616) q[1];
sx q[1];
rz(-1.1198977) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84150746) q[3];
sx q[3];
rz(-0.91472799) q[3];
sx q[3];
rz(-0.63667242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.67316002) q[2];
sx q[2];
rz(-2.1940553) q[2];
sx q[2];
rz(2.6426278) q[2];
rz(0.4237825) q[3];
sx q[3];
rz(-1.2200004) q[3];
sx q[3];
rz(3.1288778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38133165) q[0];
sx q[0];
rz(-1.6853764) q[0];
sx q[0];
rz(0.99200845) q[0];
rz(-2.965062) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(-0.36546779) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7468837) q[0];
sx q[0];
rz(-2.611126) q[0];
sx q[0];
rz(-0.51325004) q[0];
x q[1];
rz(-0.89895504) q[2];
sx q[2];
rz(-1.0869173) q[2];
sx q[2];
rz(-1.3324141) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.47479113) q[1];
sx q[1];
rz(-0.7102237) q[1];
sx q[1];
rz(2.6067961) q[1];
rz(-pi) q[2];
rz(-1.69181) q[3];
sx q[3];
rz(-1.3711208) q[3];
sx q[3];
rz(-1.4724842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3084597) q[2];
sx q[2];
rz(-2.5434727) q[2];
sx q[2];
rz(2.0642521) q[2];
rz(2.9789553) q[3];
sx q[3];
rz(-1.2626303) q[3];
sx q[3];
rz(-1.8796273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66115528) q[0];
sx q[0];
rz(-0.85395956) q[0];
sx q[0];
rz(-0.61082947) q[0];
rz(-1.707533) q[1];
sx q[1];
rz(-1.779498) q[1];
sx q[1];
rz(-0.45062137) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6621646) q[0];
sx q[0];
rz(-3.056884) q[0];
sx q[0];
rz(1.1685632) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8350416) q[2];
sx q[2];
rz(-1.2129285) q[2];
sx q[2];
rz(1.4150261) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8038517) q[1];
sx q[1];
rz(-0.9461113) q[1];
sx q[1];
rz(1.1869663) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87617989) q[3];
sx q[3];
rz(-1.3324705) q[3];
sx q[3];
rz(2.5755469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0314363) q[2];
sx q[2];
rz(-1.0855731) q[2];
sx q[2];
rz(-1.191656) q[2];
rz(0.84447652) q[3];
sx q[3];
rz(-0.28544393) q[3];
sx q[3];
rz(-0.58117956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0464756) q[0];
sx q[0];
rz(-1.7549544) q[0];
sx q[0];
rz(2.0003831) q[0];
rz(1.8449239) q[1];
sx q[1];
rz(-2.7099687) q[1];
sx q[1];
rz(2.8362714) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6396219) q[0];
sx q[0];
rz(-1.0143711) q[0];
sx q[0];
rz(1.200838) q[0];
x q[1];
rz(-2.6778053) q[2];
sx q[2];
rz(-1.1897435) q[2];
sx q[2];
rz(0.87069893) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7118453) q[1];
sx q[1];
rz(-1.4652805) q[1];
sx q[1];
rz(2.3464416) q[1];
rz(0.78379811) q[3];
sx q[3];
rz(-0.67516203) q[3];
sx q[3];
rz(-0.53962196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2083464) q[2];
sx q[2];
rz(-1.5657921) q[2];
sx q[2];
rz(-2.7944881) q[2];
rz(0.38337213) q[3];
sx q[3];
rz(-0.88630668) q[3];
sx q[3];
rz(-1.0679519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40889302) q[0];
sx q[0];
rz(-0.86298958) q[0];
sx q[0];
rz(-0.53891671) q[0];
rz(-1.7319038) q[1];
sx q[1];
rz(-2.6795645) q[1];
sx q[1];
rz(-2.7358823) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6255394) q[0];
sx q[0];
rz(-2.2279983) q[0];
sx q[0];
rz(1.9198734) q[0];
rz(1.97498) q[2];
sx q[2];
rz(-0.6037854) q[2];
sx q[2];
rz(0.052978901) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5768482) q[1];
sx q[1];
rz(-1.9124568) q[1];
sx q[1];
rz(-2.8738352) q[1];
x q[2];
rz(1.3375072) q[3];
sx q[3];
rz(-2.3217215) q[3];
sx q[3];
rz(1.6398721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6695909) q[2];
sx q[2];
rz(-2.3386023) q[2];
sx q[2];
rz(-1.4542788) q[2];
rz(0.83868319) q[3];
sx q[3];
rz(-1.8554153) q[3];
sx q[3];
rz(-1.3735166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8307761) q[0];
sx q[0];
rz(-0.61284471) q[0];
sx q[0];
rz(2.7344761) q[0];
rz(0.85917568) q[1];
sx q[1];
rz(-1.5770864) q[1];
sx q[1];
rz(1.8242594) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9068127) q[0];
sx q[0];
rz(-0.54348511) q[0];
sx q[0];
rz(1.5275005) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98951927) q[2];
sx q[2];
rz(-2.1534854) q[2];
sx q[2];
rz(-1.6473824) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0504005) q[1];
sx q[1];
rz(-2.4064153) q[1];
sx q[1];
rz(1.5666195) q[1];
x q[2];
rz(-0.42767151) q[3];
sx q[3];
rz(-1.0276405) q[3];
sx q[3];
rz(-0.021051858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.09981) q[2];
sx q[2];
rz(-0.87974292) q[2];
sx q[2];
rz(2.2463634) q[2];
rz(-1.9334531) q[3];
sx q[3];
rz(-2.4619305) q[3];
sx q[3];
rz(-2.4334548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5773425) q[0];
sx q[0];
rz(-2.0944493) q[0];
sx q[0];
rz(-2.0136925) q[0];
rz(0.212542) q[1];
sx q[1];
rz(-1.3393341) q[1];
sx q[1];
rz(1.9099265) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82176498) q[0];
sx q[0];
rz(-1.2152831) q[0];
sx q[0];
rz(-2.6000644) q[0];
rz(-pi) q[1];
rz(-0.77881323) q[2];
sx q[2];
rz(-2.3328569) q[2];
sx q[2];
rz(0.29701172) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5994508) q[1];
sx q[1];
rz(-0.25240024) q[1];
sx q[1];
rz(-2.7155994) q[1];
rz(-pi) q[2];
rz(3.1197458) q[3];
sx q[3];
rz(-0.96641809) q[3];
sx q[3];
rz(1.5165839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.82905519) q[2];
sx q[2];
rz(-0.91784358) q[2];
sx q[2];
rz(0.60116872) q[2];
rz(0.51856315) q[3];
sx q[3];
rz(-0.31373599) q[3];
sx q[3];
rz(1.4275838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9940051) q[0];
sx q[0];
rz(-0.00031358263) q[0];
sx q[0];
rz(0.073154733) q[0];
rz(-2.1144497) q[1];
sx q[1];
rz(-1.606753) q[1];
sx q[1];
rz(-0.59246078) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97544599) q[0];
sx q[0];
rz(-1.6172292) q[0];
sx q[0];
rz(-0.59068824) q[0];
rz(-pi) q[1];
rz(2.5657797) q[2];
sx q[2];
rz(-2.0137219) q[2];
sx q[2];
rz(1.5720313) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5432424) q[1];
sx q[1];
rz(-2.228884) q[1];
sx q[1];
rz(2.7212218) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98832066) q[3];
sx q[3];
rz(-2.6675825) q[3];
sx q[3];
rz(1.7085027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5804194) q[2];
sx q[2];
rz(-1.6757123) q[2];
sx q[2];
rz(0.44211659) q[2];
rz(0.90028611) q[3];
sx q[3];
rz(-1.0628275) q[3];
sx q[3];
rz(-3.1282848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028037926) q[0];
sx q[0];
rz(-0.91911626) q[0];
sx q[0];
rz(-1.8369209) q[0];
rz(1.4470944) q[1];
sx q[1];
rz(-1.1728975) q[1];
sx q[1];
rz(0.56217271) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5931372) q[0];
sx q[0];
rz(-0.89567417) q[0];
sx q[0];
rz(-2.0329342) q[0];
rz(-pi) q[1];
x q[1];
rz(0.093900605) q[2];
sx q[2];
rz(-2.0996089) q[2];
sx q[2];
rz(1.7319861) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.11645424) q[1];
sx q[1];
rz(-1.0944195) q[1];
sx q[1];
rz(-0.78671793) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1702873) q[3];
sx q[3];
rz(-1.1937965) q[3];
sx q[3];
rz(-1.4124944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8582981) q[2];
sx q[2];
rz(-1.4026182) q[2];
sx q[2];
rz(2.1600058) q[2];
rz(-0.0043407241) q[3];
sx q[3];
rz(-1.0467632) q[3];
sx q[3];
rz(-0.59345746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-1.3146064) q[0];
sx q[0];
rz(-1.2207323) q[0];
sx q[0];
rz(-3.0944968) q[0];
rz(1.3866407) q[1];
sx q[1];
rz(-1.9853233) q[1];
sx q[1];
rz(3.093739) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78794107) q[0];
sx q[0];
rz(-1.5418058) q[0];
sx q[0];
rz(-2.2257462) q[0];
rz(1.4638897) q[2];
sx q[2];
rz(-2.6968743) q[2];
sx q[2];
rz(-2.5267548) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1298444) q[1];
sx q[1];
rz(-0.99600345) q[1];
sx q[1];
rz(2.6227436) q[1];
rz(-pi) q[2];
rz(-3.0129274) q[3];
sx q[3];
rz(-0.4056969) q[3];
sx q[3];
rz(2.2588244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.64426595) q[2];
sx q[2];
rz(-1.1844119) q[2];
sx q[2];
rz(2.2052374) q[2];
rz(2.3406773) q[3];
sx q[3];
rz(-0.51910669) q[3];
sx q[3];
rz(-0.67805964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3180278) q[0];
sx q[0];
rz(-0.93946539) q[0];
sx q[0];
rz(-2.9325624) q[0];
rz(-0.50280747) q[1];
sx q[1];
rz(-1.3192758) q[1];
sx q[1];
rz(-2.1154138) q[1];
rz(-0.035632523) q[2];
sx q[2];
rz(-0.55152668) q[2];
sx q[2];
rz(1.3338911) q[2];
rz(-0.83947635) q[3];
sx q[3];
rz(-2.6203733) q[3];
sx q[3];
rz(0.93422191) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];