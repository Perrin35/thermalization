OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3553319) q[0];
sx q[0];
rz(-3.0769899) q[0];
sx q[0];
rz(-0.021615418) q[0];
rz(2.1463483) q[1];
sx q[1];
rz(-1.8145476) q[1];
sx q[1];
rz(-1.8099161) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0776805) q[0];
sx q[0];
rz(-0.67718107) q[0];
sx q[0];
rz(1.022524) q[0];
rz(-pi) q[1];
rz(1.1515491) q[2];
sx q[2];
rz(-0.24818072) q[2];
sx q[2];
rz(-0.56564769) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.057030765) q[1];
sx q[1];
rz(-1.7006526) q[1];
sx q[1];
rz(-2.3741541) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33820037) q[3];
sx q[3];
rz(-0.97465289) q[3];
sx q[3];
rz(-2.3232834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8109479) q[2];
sx q[2];
rz(-1.4928879) q[2];
sx q[2];
rz(-2.5374106) q[2];
rz(2.1172681) q[3];
sx q[3];
rz(-1.9842792) q[3];
sx q[3];
rz(2.0143051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8169096) q[0];
sx q[0];
rz(-3.1284101) q[0];
sx q[0];
rz(-2.0781793) q[0];
rz(0.88513199) q[1];
sx q[1];
rz(-1.5567895) q[1];
sx q[1];
rz(0.0016454776) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20365276) q[0];
sx q[0];
rz(-1.5814039) q[0];
sx q[0];
rz(1.525735) q[0];
x q[1];
rz(-2.4992141) q[2];
sx q[2];
rz(-1.7787691) q[2];
sx q[2];
rz(3.1285398) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.26047036) q[1];
sx q[1];
rz(-2.3350041) q[1];
sx q[1];
rz(-1.2299728) q[1];
x q[2];
rz(-1.315829) q[3];
sx q[3];
rz(-1.4062738) q[3];
sx q[3];
rz(-2.4428575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1559747) q[2];
sx q[2];
rz(-1.6321471) q[2];
sx q[2];
rz(-2.4334811) q[2];
rz(1.0937141) q[3];
sx q[3];
rz(-1.3392859) q[3];
sx q[3];
rz(0.99350199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22096069) q[0];
sx q[0];
rz(-1.7376124) q[0];
sx q[0];
rz(2.3143342) q[0];
rz(-0.0050841252) q[1];
sx q[1];
rz(-1.2132443) q[1];
sx q[1];
rz(-2.0522096) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0076865772) q[0];
sx q[0];
rz(-1.0431164) q[0];
sx q[0];
rz(-2.6469127) q[0];
rz(-pi) q[1];
rz(-1.1459648) q[2];
sx q[2];
rz(-2.1263188) q[2];
sx q[2];
rz(2.5512763) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3315862) q[1];
sx q[1];
rz(-0.3948822) q[1];
sx q[1];
rz(1.7942795) q[1];
rz(-pi) q[2];
rz(-0.79077625) q[3];
sx q[3];
rz(-1.1547935) q[3];
sx q[3];
rz(0.63000597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0777145) q[2];
sx q[2];
rz(-2.2641247) q[2];
sx q[2];
rz(0.88469488) q[2];
rz(1.9583154) q[3];
sx q[3];
rz(-1.3180472) q[3];
sx q[3];
rz(1.4900835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5692212) q[0];
sx q[0];
rz(-2.5771993) q[0];
sx q[0];
rz(2.4147721) q[0];
rz(2.3379393) q[1];
sx q[1];
rz(-2.0539961) q[1];
sx q[1];
rz(-2.7817536) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026304631) q[0];
sx q[0];
rz(-1.3867154) q[0];
sx q[0];
rz(2.227965) q[0];
x q[1];
rz(-1.6845409) q[2];
sx q[2];
rz(-2.2188088) q[2];
sx q[2];
rz(-1.5930454) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6105146) q[1];
sx q[1];
rz(-0.58137608) q[1];
sx q[1];
rz(-0.73552144) q[1];
rz(1.9197649) q[3];
sx q[3];
rz(-1.8604606) q[3];
sx q[3];
rz(1.3834013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5741253) q[2];
sx q[2];
rz(-1.1541157) q[2];
sx q[2];
rz(-0.7652258) q[2];
rz(-0.75677538) q[3];
sx q[3];
rz(-0.59195834) q[3];
sx q[3];
rz(-2.9005907) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9969479) q[0];
sx q[0];
rz(-0.47676555) q[0];
sx q[0];
rz(1.416052) q[0];
rz(-0.36711806) q[1];
sx q[1];
rz(-1.3857625) q[1];
sx q[1];
rz(2.1062772) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9177592) q[0];
sx q[0];
rz(-0.98063722) q[0];
sx q[0];
rz(0.53442861) q[0];
rz(-1.5315227) q[2];
sx q[2];
rz(-2.0619259) q[2];
sx q[2];
rz(-0.66205762) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1977735) q[1];
sx q[1];
rz(-3.0316331) q[1];
sx q[1];
rz(0.52093671) q[1];
rz(-pi) q[2];
rz(1.7807998) q[3];
sx q[3];
rz(-1.3011419) q[3];
sx q[3];
rz(0.3013914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7845903) q[2];
sx q[2];
rz(-1.720153) q[2];
sx q[2];
rz(-2.3727097) q[2];
rz(-0.33603493) q[3];
sx q[3];
rz(-0.78032812) q[3];
sx q[3];
rz(-1.4679573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1795905) q[0];
sx q[0];
rz(-2.5080894) q[0];
sx q[0];
rz(1.1451716) q[0];
rz(-1.1046474) q[1];
sx q[1];
rz(-1.9112174) q[1];
sx q[1];
rz(-2.9343658) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80400318) q[0];
sx q[0];
rz(-2.1875256) q[0];
sx q[0];
rz(-1.7417275) q[0];
rz(-2.5404262) q[2];
sx q[2];
rz(-1.1784369) q[2];
sx q[2];
rz(0.62077921) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8689649) q[1];
sx q[1];
rz(-0.3513063) q[1];
sx q[1];
rz(2.3428194) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81760041) q[3];
sx q[3];
rz(-1.451965) q[3];
sx q[3];
rz(-0.51018754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.8032288) q[2];
sx q[2];
rz(-2.2613566) q[2];
sx q[2];
rz(2.4198789) q[2];
rz(1.8668113) q[3];
sx q[3];
rz(-1.5721679) q[3];
sx q[3];
rz(0.26003626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8969144) q[0];
sx q[0];
rz(-1.5748064) q[0];
sx q[0];
rz(-0.73202837) q[0];
rz(-0.0094982068) q[1];
sx q[1];
rz(-2.5932725) q[1];
sx q[1];
rz(0.2917372) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1222629) q[0];
sx q[0];
rz(-2.7195192) q[0];
sx q[0];
rz(-0.12515573) q[0];
rz(-pi) q[1];
rz(-1.2417492) q[2];
sx q[2];
rz(-2.8223158) q[2];
sx q[2];
rz(0.74908756) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8768423) q[1];
sx q[1];
rz(-1.2683588) q[1];
sx q[1];
rz(1.7251863) q[1];
rz(-pi) q[2];
rz(-0.50645701) q[3];
sx q[3];
rz(-2.5589057) q[3];
sx q[3];
rz(1.9401693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.26414028) q[2];
sx q[2];
rz(-1.143127) q[2];
sx q[2];
rz(-2.3042802) q[2];
rz(-1.1710179) q[3];
sx q[3];
rz(-1.5191017) q[3];
sx q[3];
rz(3.0795735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0983122) q[0];
sx q[0];
rz(-2.4536112) q[0];
sx q[0];
rz(0.6638546) q[0];
rz(0.10617667) q[1];
sx q[1];
rz(-2.5352434) q[1];
sx q[1];
rz(2.1829139) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7230969) q[0];
sx q[0];
rz(-2.3765058) q[0];
sx q[0];
rz(-0.78406783) q[0];
rz(-pi) q[1];
rz(2.8616222) q[2];
sx q[2];
rz(-2.0915871) q[2];
sx q[2];
rz(2.1998646) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9194591) q[1];
sx q[1];
rz(-2.6047915) q[1];
sx q[1];
rz(1.7807351) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.18387353) q[3];
sx q[3];
rz(-0.5920147) q[3];
sx q[3];
rz(1.5873991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.133193) q[2];
sx q[2];
rz(-1.2352751) q[2];
sx q[2];
rz(0.91910249) q[2];
rz(1.3778) q[3];
sx q[3];
rz(-1.2104687) q[3];
sx q[3];
rz(1.9581883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8063426) q[0];
sx q[0];
rz(-2.2882473) q[0];
sx q[0];
rz(0.43689716) q[0];
rz(2.4412952) q[1];
sx q[1];
rz(-1.7179787) q[1];
sx q[1];
rz(-0.46554309) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1818905) q[0];
sx q[0];
rz(-0.86599106) q[0];
sx q[0];
rz(-1.6824791) q[0];
rz(-1.9829468) q[2];
sx q[2];
rz(-1.4530384) q[2];
sx q[2];
rz(2.462537) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0034345) q[1];
sx q[1];
rz(-0.58193602) q[1];
sx q[1];
rz(-0.65852965) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25119987) q[3];
sx q[3];
rz(-1.2217055) q[3];
sx q[3];
rz(2.420345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7312701) q[2];
sx q[2];
rz(-0.53933829) q[2];
sx q[2];
rz(-1.643606) q[2];
rz(2.9368029) q[3];
sx q[3];
rz(-1.0054761) q[3];
sx q[3];
rz(0.27206102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64514226) q[0];
sx q[0];
rz(-2.5800939) q[0];
sx q[0];
rz(-1.1219332) q[0];
rz(2.3902068) q[1];
sx q[1];
rz(-1.6151927) q[1];
sx q[1];
rz(2.5591992) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1727027) q[0];
sx q[0];
rz(-0.70969289) q[0];
sx q[0];
rz(-1.2560647) q[0];
x q[1];
rz(-2.0612129) q[2];
sx q[2];
rz(-1.4368125) q[2];
sx q[2];
rz(1.4659363) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.45422428) q[1];
sx q[1];
rz(-0.2174938) q[1];
sx q[1];
rz(0.99517676) q[1];
x q[2];
rz(-0.21721812) q[3];
sx q[3];
rz(-1.5981042) q[3];
sx q[3];
rz(0.56009968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1311243) q[2];
sx q[2];
rz(-2.1605587) q[2];
sx q[2];
rz(-2.3790512) q[2];
rz(1.7307581) q[3];
sx q[3];
rz(-0.91791955) q[3];
sx q[3];
rz(-1.5677174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.0653771) q[0];
sx q[0];
rz(-2.1552754) q[0];
sx q[0];
rz(-1.7393204) q[0];
rz(-1.8021884) q[1];
sx q[1];
rz(-1.4504455) q[1];
sx q[1];
rz(1.6557678) q[1];
rz(-2.8600678) q[2];
sx q[2];
rz(-0.83419656) q[2];
sx q[2];
rz(0.4783334) q[2];
rz(-1.0675666) q[3];
sx q[3];
rz(-1.9964841) q[3];
sx q[3];
rz(-1.869429) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
