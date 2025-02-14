OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.96149093) q[0];
sx q[0];
rz(-2.5749126) q[0];
sx q[0];
rz(0.79071796) q[0];
rz(2.2136731) q[1];
sx q[1];
rz(3.1766422) q[1];
sx q[1];
rz(9.1215134) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.431031) q[0];
sx q[0];
rz(-2.2236784) q[0];
sx q[0];
rz(-1.2827355) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9531115) q[2];
sx q[2];
rz(-1.9046724) q[2];
sx q[2];
rz(-0.44313639) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11285148) q[1];
sx q[1];
rz(-1.0489221) q[1];
sx q[1];
rz(2.584409) q[1];
rz(-pi) q[2];
rz(0.73533006) q[3];
sx q[3];
rz(-1.5092675) q[3];
sx q[3];
rz(1.3593591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7683893) q[2];
sx q[2];
rz(-0.36036569) q[2];
sx q[2];
rz(2.4599794) q[2];
rz(2.5810589) q[3];
sx q[3];
rz(-2.3571641) q[3];
sx q[3];
rz(2.4973629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19175567) q[0];
sx q[0];
rz(-2.2332709) q[0];
sx q[0];
rz(-2.7857696) q[0];
rz(2.8882354) q[1];
sx q[1];
rz(-2.2919877) q[1];
sx q[1];
rz(1.2454978) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1863128) q[0];
sx q[0];
rz(-1.3471502) q[0];
sx q[0];
rz(2.6564084) q[0];
rz(-pi) q[1];
rz(2.5077057) q[2];
sx q[2];
rz(-1.1340464) q[2];
sx q[2];
rz(-0.89981198) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.624234) q[1];
sx q[1];
rz(-1.1453724) q[1];
sx q[1];
rz(2.5555771) q[1];
rz(-2.6368876) q[3];
sx q[3];
rz(-2.0619806) q[3];
sx q[3];
rz(2.1757656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7443098) q[2];
sx q[2];
rz(-2.3645568) q[2];
sx q[2];
rz(3.0139319) q[2];
rz(0.67576659) q[3];
sx q[3];
rz(-2.6830169) q[3];
sx q[3];
rz(2.7202386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0449988) q[0];
sx q[0];
rz(-2.8969722) q[0];
sx q[0];
rz(-2.007572) q[0];
rz(0.89316142) q[1];
sx q[1];
rz(-0.90407073) q[1];
sx q[1];
rz(1.8697416) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9006827) q[0];
sx q[0];
rz(-2.5782452) q[0];
sx q[0];
rz(1.6202089) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8142305) q[2];
sx q[2];
rz(-2.8231502) q[2];
sx q[2];
rz(-1.8093579) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0772897) q[1];
sx q[1];
rz(-1.3659119) q[1];
sx q[1];
rz(-2.2868392) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2168733) q[3];
sx q[3];
rz(-2.8689403) q[3];
sx q[3];
rz(-0.70875185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9280055) q[2];
sx q[2];
rz(-0.86557937) q[2];
sx q[2];
rz(-1.340284) q[2];
rz(1.853893) q[3];
sx q[3];
rz(-1.6940593) q[3];
sx q[3];
rz(-2.652216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61557788) q[0];
sx q[0];
rz(-2.5647793) q[0];
sx q[0];
rz(-2.4082129) q[0];
rz(-2.3461657) q[1];
sx q[1];
rz(-2.9454102) q[1];
sx q[1];
rz(-1.1682074) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0709658) q[0];
sx q[0];
rz(-0.42558181) q[0];
sx q[0];
rz(-1.312567) q[0];
rz(-pi) q[1];
x q[1];
rz(2.204037) q[2];
sx q[2];
rz(-2.0481234) q[2];
sx q[2];
rz(0.027743159) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4609264) q[1];
sx q[1];
rz(-3.089) q[1];
sx q[1];
rz(-0.61477964) q[1];
x q[2];
rz(0.73428728) q[3];
sx q[3];
rz(-2.3288615) q[3];
sx q[3];
rz(-1.6664291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1802133) q[2];
sx q[2];
rz(-0.81575477) q[2];
sx q[2];
rz(0.95480314) q[2];
rz(2.09156) q[3];
sx q[3];
rz(-0.78767109) q[3];
sx q[3];
rz(-2.5900456) q[3];
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
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10057218) q[0];
sx q[0];
rz(-0.53871483) q[0];
sx q[0];
rz(2.4464497) q[0];
rz(-1.9255385) q[1];
sx q[1];
rz(-2.1247517) q[1];
sx q[1];
rz(-3.1207808) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7835158) q[0];
sx q[0];
rz(-2.4672271) q[0];
sx q[0];
rz(-1.800531) q[0];
x q[1];
rz(-1.2131507) q[2];
sx q[2];
rz(-0.76947533) q[2];
sx q[2];
rz(-2.3703252) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4592308) q[1];
sx q[1];
rz(-1.8098847) q[1];
sx q[1];
rz(0.011929529) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6324051) q[3];
sx q[3];
rz(-1.6635259) q[3];
sx q[3];
rz(-1.9237067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1702599) q[2];
sx q[2];
rz(-2.8016165) q[2];
sx q[2];
rz(-0.60410947) q[2];
rz(-2.4938834) q[3];
sx q[3];
rz(-2.6205687) q[3];
sx q[3];
rz(0.46085301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0172417) q[0];
sx q[0];
rz(-2.0483973) q[0];
sx q[0];
rz(-2.7830615) q[0];
rz(-0.56309807) q[1];
sx q[1];
rz(-0.9022572) q[1];
sx q[1];
rz(0.30555746) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3747923) q[0];
sx q[0];
rz(-1.2106441) q[0];
sx q[0];
rz(2.4712579) q[0];
x q[1];
rz(1.7420058) q[2];
sx q[2];
rz(-1.3911026) q[2];
sx q[2];
rz(3.0269738) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8628775) q[1];
sx q[1];
rz(-2.3722665) q[1];
sx q[1];
rz(0.74689405) q[1];
rz(-pi) q[2];
x q[2];
rz(0.077353954) q[3];
sx q[3];
rz(-1.6408444) q[3];
sx q[3];
rz(-1.182568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.87865889) q[2];
sx q[2];
rz(-0.96076751) q[2];
sx q[2];
rz(-0.052138694) q[2];
rz(-0.42929286) q[3];
sx q[3];
rz(-1.41058) q[3];
sx q[3];
rz(-0.2521635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4572064) q[0];
sx q[0];
rz(-2.6173213) q[0];
sx q[0];
rz(0.23505178) q[0];
rz(2.5382407) q[1];
sx q[1];
rz(-2.3124606) q[1];
sx q[1];
rz(2.5650909) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7991195) q[0];
sx q[0];
rz(-0.89460374) q[0];
sx q[0];
rz(-0.66838092) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2320249) q[2];
sx q[2];
rz(-0.82451754) q[2];
sx q[2];
rz(0.19245779) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.41771463) q[1];
sx q[1];
rz(-0.52199328) q[1];
sx q[1];
rz(-0.88845171) q[1];
x q[2];
rz(2.0527538) q[3];
sx q[3];
rz(-1.6982984) q[3];
sx q[3];
rz(-1.8021999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.06360402) q[2];
sx q[2];
rz(-2.3006907) q[2];
sx q[2];
rz(-1.3464751) q[2];
rz(-0.50955647) q[3];
sx q[3];
rz(-1.8815123) q[3];
sx q[3];
rz(-0.30663651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52617306) q[0];
sx q[0];
rz(-2.058448) q[0];
sx q[0];
rz(2.7857067) q[0];
rz(-2.4399759) q[1];
sx q[1];
rz(-2.9569148) q[1];
sx q[1];
rz(-2.1772749) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6414757) q[0];
sx q[0];
rz(-1.7044412) q[0];
sx q[0];
rz(2.0290613) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2779916) q[2];
sx q[2];
rz(-1.5902012) q[2];
sx q[2];
rz(0.524117) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1279932) q[1];
sx q[1];
rz(-1.0759872) q[1];
sx q[1];
rz(1.593099) q[1];
rz(2.8300222) q[3];
sx q[3];
rz(-1.6487247) q[3];
sx q[3];
rz(-0.0078709906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.648916) q[2];
sx q[2];
rz(-2.4339088) q[2];
sx q[2];
rz(0.81563449) q[2];
rz(0.91551578) q[3];
sx q[3];
rz(-2.2736277) q[3];
sx q[3];
rz(-2.9149132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8409519) q[0];
sx q[0];
rz(-0.17423593) q[0];
sx q[0];
rz(1.9006282) q[0];
rz(-3.0366483) q[1];
sx q[1];
rz(-2.8009156) q[1];
sx q[1];
rz(2.6822283) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1007538) q[0];
sx q[0];
rz(-1.213165) q[0];
sx q[0];
rz(0.25672428) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2392339) q[2];
sx q[2];
rz(-1.3156097) q[2];
sx q[2];
rz(-0.58411264) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1314754) q[1];
sx q[1];
rz(-1.7820211) q[1];
sx q[1];
rz(-0.44472163) q[1];
rz(-pi) q[2];
rz(2.9271803) q[3];
sx q[3];
rz(-0.58480703) q[3];
sx q[3];
rz(2.6406276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0910054) q[2];
sx q[2];
rz(-2.4531187) q[2];
sx q[2];
rz(-0.090593226) q[2];
rz(2.374384) q[3];
sx q[3];
rz(-1.4827261) q[3];
sx q[3];
rz(-1.4513133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7200658) q[0];
sx q[0];
rz(-0.92362112) q[0];
sx q[0];
rz(-1.5371171) q[0];
rz(-0.9556669) q[1];
sx q[1];
rz(-1.6803398) q[1];
sx q[1];
rz(2.59424) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5621278) q[0];
sx q[0];
rz(-1.5766931) q[0];
sx q[0];
rz(1.4910758) q[0];
rz(2.5432406) q[2];
sx q[2];
rz(-0.6171591) q[2];
sx q[2];
rz(2.0218818) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.66710752) q[1];
sx q[1];
rz(-0.82208668) q[1];
sx q[1];
rz(-2.8182823) q[1];
rz(1.0676974) q[3];
sx q[3];
rz(-0.92313952) q[3];
sx q[3];
rz(1.5474418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0305816) q[2];
sx q[2];
rz(-2.9647398) q[2];
sx q[2];
rz(-0.64876968) q[2];
rz(-2.9099921) q[3];
sx q[3];
rz(-0.88445556) q[3];
sx q[3];
rz(-0.086656682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34407525) q[0];
sx q[0];
rz(-1.9301013) q[0];
sx q[0];
rz(1.7420266) q[0];
rz(-0.65921417) q[1];
sx q[1];
rz(-1.8623687) q[1];
sx q[1];
rz(1.6946793) q[1];
rz(-1.7991015) q[2];
sx q[2];
rz(-1.5226428) q[2];
sx q[2];
rz(-1.5626283) q[2];
rz(1.4249887) q[3];
sx q[3];
rz(-1.6248228) q[3];
sx q[3];
rz(1.1820033) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
