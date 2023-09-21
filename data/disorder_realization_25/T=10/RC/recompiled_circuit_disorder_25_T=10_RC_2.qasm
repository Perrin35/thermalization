OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0796354) q[0];
sx q[0];
rz(-2.893653) q[0];
sx q[0];
rz(0.82013446) q[0];
rz(-3.1007383) q[1];
sx q[1];
rz(-0.78894579) q[1];
sx q[1];
rz(3.0537002) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.041931) q[0];
sx q[0];
rz(-0.66177216) q[0];
sx q[0];
rz(-2.3753138) q[0];
rz(-2.7152039) q[2];
sx q[2];
rz(-1.083056) q[2];
sx q[2];
rz(2.3418155) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3275798) q[1];
sx q[1];
rz(-2.6373133) q[1];
sx q[1];
rz(-0.0045279702) q[1];
x q[2];
rz(0.0058562972) q[3];
sx q[3];
rz(-1.8294888) q[3];
sx q[3];
rz(-2.4394025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1408046) q[2];
sx q[2];
rz(-1.452383) q[2];
sx q[2];
rz(2.5722356) q[2];
rz(-1.5287483) q[3];
sx q[3];
rz(-2.5053535) q[3];
sx q[3];
rz(1.3585842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59812087) q[0];
sx q[0];
rz(-2.9765029) q[0];
sx q[0];
rz(0.55364451) q[0];
rz(-1.9042632) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(-1.233261) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.935826) q[0];
sx q[0];
rz(-0.92139771) q[0];
sx q[0];
rz(2.4539024) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4428667) q[2];
sx q[2];
rz(-2.2413841) q[2];
sx q[2];
rz(0.48534457) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1439647) q[1];
sx q[1];
rz(-1.7768754) q[1];
sx q[1];
rz(-3.1254083) q[1];
rz(-0.29970701) q[3];
sx q[3];
rz(-2.4518659) q[3];
sx q[3];
rz(2.0852065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1374986) q[2];
sx q[2];
rz(-1.4592905) q[2];
sx q[2];
rz(-3.1393576) q[2];
rz(-0.8301174) q[3];
sx q[3];
rz(-2.3486962) q[3];
sx q[3];
rz(-1.2494276) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24901351) q[0];
sx q[0];
rz(-0.61616388) q[0];
sx q[0];
rz(-2.2900443) q[0];
rz(-2.3705204) q[1];
sx q[1];
rz(-1.4665736) q[1];
sx q[1];
rz(0.071203701) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93807968) q[0];
sx q[0];
rz(-2.8305614) q[0];
sx q[0];
rz(-0.87840338) q[0];
rz(-pi) q[1];
rz(0.26239563) q[2];
sx q[2];
rz(-1.1294282) q[2];
sx q[2];
rz(1.0052094) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60656602) q[1];
sx q[1];
rz(-2.2177794) q[1];
sx q[1];
rz(1.7391298) q[1];
x q[2];
rz(0.76020469) q[3];
sx q[3];
rz(-2.1225727) q[3];
sx q[3];
rz(2.9387568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3123793) q[2];
sx q[2];
rz(-0.92386121) q[2];
sx q[2];
rz(2.6181347) q[2];
rz(-0.33189014) q[3];
sx q[3];
rz(-3.0811946) q[3];
sx q[3];
rz(2.6383242) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7317384) q[0];
sx q[0];
rz(-2.8926352) q[0];
sx q[0];
rz(2.3160034) q[0];
rz(-1.9748953) q[1];
sx q[1];
rz(-0.86667934) q[1];
sx q[1];
rz(-1.694214) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3576413) q[0];
sx q[0];
rz(-0.28342993) q[0];
sx q[0];
rz(2.1901603) q[0];
rz(-pi) q[1];
rz(-3.1324603) q[2];
sx q[2];
rz(-1.7491241) q[2];
sx q[2];
rz(0.4724617) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.73562685) q[1];
sx q[1];
rz(-1.7298797) q[1];
sx q[1];
rz(-2.2707978) q[1];
rz(-pi) q[2];
rz(2.7778266) q[3];
sx q[3];
rz(-1.7978661) q[3];
sx q[3];
rz(2.098339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.46431413) q[2];
sx q[2];
rz(-1.7769287) q[2];
sx q[2];
rz(-0.14492598) q[2];
rz(-2.1302917) q[3];
sx q[3];
rz(-2.5144808) q[3];
sx q[3];
rz(3.1269126) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99575627) q[0];
sx q[0];
rz(-0.053481426) q[0];
sx q[0];
rz(-0.68403912) q[0];
rz(-1.9513291) q[1];
sx q[1];
rz(-0.66527706) q[1];
sx q[1];
rz(1.2129983) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1367462) q[0];
sx q[0];
rz(-2.6803225) q[0];
sx q[0];
rz(-2.524551) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62972516) q[2];
sx q[2];
rz(-1.0169528) q[2];
sx q[2];
rz(2.6385006) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7291521) q[1];
sx q[1];
rz(-1.3153331) q[1];
sx q[1];
rz(-1.4490119) q[1];
rz(-pi) q[2];
rz(1.3092242) q[3];
sx q[3];
rz(-1.7422758) q[3];
sx q[3];
rz(1.8985626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6325536) q[2];
sx q[2];
rz(-1.8703439) q[2];
sx q[2];
rz(0.63009134) q[2];
rz(0.57224327) q[3];
sx q[3];
rz(-2.4961491) q[3];
sx q[3];
rz(-2.2563289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1155788) q[0];
sx q[0];
rz(-0.84340874) q[0];
sx q[0];
rz(2.8552326) q[0];
rz(-2.5419366) q[1];
sx q[1];
rz(-1.1279794) q[1];
sx q[1];
rz(-0.58475959) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2858546) q[0];
sx q[0];
rz(-1.8830839) q[0];
sx q[0];
rz(0.62494846) q[0];
rz(-1.5969959) q[2];
sx q[2];
rz(-2.1850586) q[2];
sx q[2];
rz(-0.40528497) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4510182) q[1];
sx q[1];
rz(-1.283053) q[1];
sx q[1];
rz(2.2698127) q[1];
rz(-0.29295178) q[3];
sx q[3];
rz(-0.17360273) q[3];
sx q[3];
rz(-2.0041182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9770603) q[2];
sx q[2];
rz(-2.2571199) q[2];
sx q[2];
rz(-1.5131081) q[2];
rz(-2.1785054) q[3];
sx q[3];
rz(-1.77308) q[3];
sx q[3];
rz(-2.517038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088783711) q[0];
sx q[0];
rz(-1.7099021) q[0];
sx q[0];
rz(-0.54876304) q[0];
rz(-3.0896297) q[1];
sx q[1];
rz(-2.6497662) q[1];
sx q[1];
rz(2.5351977) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.279778) q[0];
sx q[0];
rz(-1.8068411) q[0];
sx q[0];
rz(-1.0311014) q[0];
rz(0.90784448) q[2];
sx q[2];
rz(-1.2693451) q[2];
sx q[2];
rz(0.19336685) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.062321669) q[1];
sx q[1];
rz(-0.81264001) q[1];
sx q[1];
rz(-0.15516849) q[1];
x q[2];
rz(1.1612732) q[3];
sx q[3];
rz(-1.2425353) q[3];
sx q[3];
rz(1.3271774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2018044) q[2];
sx q[2];
rz(-0.99210343) q[2];
sx q[2];
rz(2.6331804) q[2];
rz(-1.1172179) q[3];
sx q[3];
rz(-2.9682187) q[3];
sx q[3];
rz(1.6569051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.33928076) q[0];
sx q[0];
rz(-2.7957714) q[0];
sx q[0];
rz(-2.3876277) q[0];
rz(0.94999653) q[1];
sx q[1];
rz(-1.4818622) q[1];
sx q[1];
rz(-2.9139013) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84275093) q[0];
sx q[0];
rz(-1.1584873) q[0];
sx q[0];
rz(3.0639248) q[0];
rz(0.44616206) q[2];
sx q[2];
rz(-1.9036306) q[2];
sx q[2];
rz(-1.3653501) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.55351174) q[1];
sx q[1];
rz(-1.7953824) q[1];
sx q[1];
rz(2.5976546) q[1];
x q[2];
rz(0.065683059) q[3];
sx q[3];
rz(-2.5816397) q[3];
sx q[3];
rz(2.5064859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.078538744) q[2];
sx q[2];
rz(-2.9324014) q[2];
sx q[2];
rz(0.83855808) q[2];
rz(-0.36455425) q[3];
sx q[3];
rz(-1.3804599) q[3];
sx q[3];
rz(2.3013039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7446328) q[0];
sx q[0];
rz(-1.3752221) q[0];
sx q[0];
rz(-2.3378085) q[0];
rz(1.0546168) q[1];
sx q[1];
rz(-2.5024253) q[1];
sx q[1];
rz(-1.1484336) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.046612) q[0];
sx q[0];
rz(-1.5946832) q[0];
sx q[0];
rz(-1.4987962) q[0];
rz(-pi) q[1];
rz(-2.7788413) q[2];
sx q[2];
rz(-1.0748378) q[2];
sx q[2];
rz(0.19778684) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.941274) q[1];
sx q[1];
rz(-0.66220821) q[1];
sx q[1];
rz(2.3108285) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5356482) q[3];
sx q[3];
rz(-1.9225549) q[3];
sx q[3];
rz(0.37946057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0917197) q[2];
sx q[2];
rz(-2.1240978) q[2];
sx q[2];
rz(-0.77825528) q[2];
rz(2.9459279) q[3];
sx q[3];
rz(-1.0396495) q[3];
sx q[3];
rz(2.1197317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8675999) q[0];
sx q[0];
rz(-0.99853981) q[0];
sx q[0];
rz(-0.25319779) q[0];
rz(-2.4018438) q[1];
sx q[1];
rz(-2.4052129) q[1];
sx q[1];
rz(-0.69828066) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88903504) q[0];
sx q[0];
rz(-1.6425942) q[0];
sx q[0];
rz(0.30307146) q[0];
rz(-pi) q[1];
x q[1];
rz(0.024784879) q[2];
sx q[2];
rz(-1.4966655) q[2];
sx q[2];
rz(-2.5075846) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0638949) q[1];
sx q[1];
rz(-0.5118891) q[1];
sx q[1];
rz(-1.1179395) q[1];
rz(-pi) q[2];
rz(-2.9000557) q[3];
sx q[3];
rz(-1.0165443) q[3];
sx q[3];
rz(-0.17061558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.11761052) q[2];
sx q[2];
rz(-0.88279804) q[2];
sx q[2];
rz(2.3013766) q[2];
rz(-2.501287) q[3];
sx q[3];
rz(-0.87602031) q[3];
sx q[3];
rz(-1.1673814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2587851) q[0];
sx q[0];
rz(-0.85492815) q[0];
sx q[0];
rz(1.7229765) q[0];
rz(0.029126833) q[1];
sx q[1];
rz(-0.11321414) q[1];
sx q[1];
rz(-1.3197457) q[1];
rz(-2.5393456) q[2];
sx q[2];
rz(-1.2802274) q[2];
sx q[2];
rz(1.0614492) q[2];
rz(-2.2610353) q[3];
sx q[3];
rz(-1.6869443) q[3];
sx q[3];
rz(-1.501367) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];