OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7944827) q[0];
sx q[0];
rz(-1.4248983) q[0];
sx q[0];
rz(0.27118924) q[0];
rz(0.36692342) q[1];
sx q[1];
rz(-2.38382) q[1];
sx q[1];
rz(1.4581207) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9957033) q[0];
sx q[0];
rz(-1.567739) q[0];
sx q[0];
rz(-1.5778592) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.053399) q[2];
sx q[2];
rz(-1.6763684) q[2];
sx q[2];
rz(-2.2029049) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7768727) q[1];
sx q[1];
rz(-2.128771) q[1];
sx q[1];
rz(-0.67529894) q[1];
rz(-pi) q[2];
rz(0.062964418) q[3];
sx q[3];
rz(-2.8573572) q[3];
sx q[3];
rz(-0.96719826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1272507) q[2];
sx q[2];
rz(-2.47561) q[2];
sx q[2];
rz(-1.8951529) q[2];
rz(-2.1261101) q[3];
sx q[3];
rz(-0.4963488) q[3];
sx q[3];
rz(0.61771667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.604973) q[0];
sx q[0];
rz(-3.0687357) q[0];
sx q[0];
rz(-0.68742043) q[0];
rz(2.5303326) q[1];
sx q[1];
rz(-1.6300853) q[1];
sx q[1];
rz(0.90868178) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9729823) q[0];
sx q[0];
rz(-1.6618722) q[0];
sx q[0];
rz(0.9206207) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15678997) q[2];
sx q[2];
rz(-2.1794381) q[2];
sx q[2];
rz(-1.5896686) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7707053) q[1];
sx q[1];
rz(-1.0418278) q[1];
sx q[1];
rz(2.4338399) q[1];
x q[2];
rz(0.55697316) q[3];
sx q[3];
rz(-1.6786298) q[3];
sx q[3];
rz(1.432991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0574135) q[2];
sx q[2];
rz(-0.97687352) q[2];
sx q[2];
rz(1.9834391) q[2];
rz(-0.17364764) q[3];
sx q[3];
rz(-1.4519139) q[3];
sx q[3];
rz(2.4524073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77370206) q[0];
sx q[0];
rz(-2.9359718) q[0];
sx q[0];
rz(2.6389627) q[0];
rz(1.3357119) q[1];
sx q[1];
rz(-0.10919658) q[1];
sx q[1];
rz(1.4044382) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40622845) q[0];
sx q[0];
rz(-1.9695008) q[0];
sx q[0];
rz(0.47946341) q[0];
rz(1.7638172) q[2];
sx q[2];
rz(-1.9380331) q[2];
sx q[2];
rz(2.6835416) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3223649) q[1];
sx q[1];
rz(-0.33999264) q[1];
sx q[1];
rz(-2.875706) q[1];
rz(-0.59111528) q[3];
sx q[3];
rz(-2.0264668) q[3];
sx q[3];
rz(-2.299813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8414843) q[2];
sx q[2];
rz(-2.4418162) q[2];
sx q[2];
rz(2.4669199) q[2];
rz(0.35215968) q[3];
sx q[3];
rz(-1.1511185) q[3];
sx q[3];
rz(-1.5153511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.801055) q[0];
sx q[0];
rz(-0.041943701) q[0];
sx q[0];
rz(-1.9757784) q[0];
rz(0.55533987) q[1];
sx q[1];
rz(-0.97709877) q[1];
sx q[1];
rz(1.4110483) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4395907) q[0];
sx q[0];
rz(-1.5541663) q[0];
sx q[0];
rz(-2.8921333) q[0];
x q[1];
rz(1.246366) q[2];
sx q[2];
rz(-2.590605) q[2];
sx q[2];
rz(1.18206) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1852604) q[1];
sx q[1];
rz(-1.0721551) q[1];
sx q[1];
rz(0.39875984) q[1];
rz(-2.9054303) q[3];
sx q[3];
rz(-0.8311701) q[3];
sx q[3];
rz(-0.96280747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0393151) q[2];
sx q[2];
rz(-2.6199876) q[2];
sx q[2];
rz(-0.48307854) q[2];
rz(0.77110243) q[3];
sx q[3];
rz(-1.5811788) q[3];
sx q[3];
rz(-1.34351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1191206) q[0];
sx q[0];
rz(-2.6799057) q[0];
sx q[0];
rz(-2.9126677) q[0];
rz(-1.4643033) q[1];
sx q[1];
rz(-2.5720282) q[1];
sx q[1];
rz(-2.6615108) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10413607) q[0];
sx q[0];
rz(-2.7711282) q[0];
sx q[0];
rz(-2.3568826) q[0];
rz(-pi) q[1];
rz(1.8304906) q[2];
sx q[2];
rz(-1.9995481) q[2];
sx q[2];
rz(-1.7130121) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6704178) q[1];
sx q[1];
rz(-1.5199317) q[1];
sx q[1];
rz(-2.7689217) q[1];
rz(-pi) q[2];
rz(2.1958417) q[3];
sx q[3];
rz(-1.4212928) q[3];
sx q[3];
rz(-0.76555071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0698801) q[2];
sx q[2];
rz(-2.0166848) q[2];
sx q[2];
rz(1.9055535) q[2];
rz(1.49336) q[3];
sx q[3];
rz(-2.3107216) q[3];
sx q[3];
rz(-1.6911471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1278095) q[0];
sx q[0];
rz(-2.4007128) q[0];
sx q[0];
rz(2.2770449) q[0];
rz(-2.3132482) q[1];
sx q[1];
rz(-1.8048077) q[1];
sx q[1];
rz(0.72992647) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8720388) q[0];
sx q[0];
rz(-2.3667496) q[0];
sx q[0];
rz(-2.274735) q[0];
rz(-pi) q[1];
rz(0.24940074) q[2];
sx q[2];
rz(-0.6847544) q[2];
sx q[2];
rz(2.8114787) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.79086458) q[1];
sx q[1];
rz(-2.5525188) q[1];
sx q[1];
rz(1.1228611) q[1];
x q[2];
rz(0.72336332) q[3];
sx q[3];
rz(-1.5302684) q[3];
sx q[3];
rz(1.877829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.56968969) q[2];
sx q[2];
rz(-0.14293417) q[2];
sx q[2];
rz(-1.8611056) q[2];
rz(1.5661543) q[3];
sx q[3];
rz(-1.0404634) q[3];
sx q[3];
rz(2.2787826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5363619) q[0];
sx q[0];
rz(-1.0030168) q[0];
sx q[0];
rz(0.59148106) q[0];
rz(-0.65762562) q[1];
sx q[1];
rz(-1.769519) q[1];
sx q[1];
rz(-0.11810158) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5786835) q[0];
sx q[0];
rz(-1.2244891) q[0];
sx q[0];
rz(1.3265394) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0131091) q[2];
sx q[2];
rz(-0.5196577) q[2];
sx q[2];
rz(-2.3574587) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76712045) q[1];
sx q[1];
rz(-0.53009696) q[1];
sx q[1];
rz(-2.9494826) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2198417) q[3];
sx q[3];
rz(-2.8650682) q[3];
sx q[3];
rz(3.0365503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20146951) q[2];
sx q[2];
rz(-2.757299) q[2];
sx q[2];
rz(-1.8817792) q[2];
rz(-1.2093557) q[3];
sx q[3];
rz(-1.3426547) q[3];
sx q[3];
rz(0.847009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2131293) q[0];
sx q[0];
rz(-2.6452112) q[0];
sx q[0];
rz(-1.553836) q[0];
rz(-1.3262879) q[1];
sx q[1];
rz(-2.1554048) q[1];
sx q[1];
rz(-0.80002588) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6094506) q[0];
sx q[0];
rz(-0.47663078) q[0];
sx q[0];
rz(2.5146288) q[0];
rz(-pi) q[1];
rz(-2.3800092) q[2];
sx q[2];
rz(-1.933145) q[2];
sx q[2];
rz(-1.5630388) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.043083755) q[1];
sx q[1];
rz(-1.0810163) q[1];
sx q[1];
rz(-0.54651641) q[1];
x q[2];
rz(0.991865) q[3];
sx q[3];
rz(-0.63922113) q[3];
sx q[3];
rz(-1.0238436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.942261) q[2];
sx q[2];
rz(-2.3776725) q[2];
sx q[2];
rz(-1.8755111) q[2];
rz(2.074504) q[3];
sx q[3];
rz(-1.7231924) q[3];
sx q[3];
rz(-3.0838695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98457321) q[0];
sx q[0];
rz(-2.3047801) q[0];
sx q[0];
rz(2.9465604) q[0];
rz(0.86206478) q[1];
sx q[1];
rz(-1.74086) q[1];
sx q[1];
rz(-0.21924266) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1006323) q[0];
sx q[0];
rz(-2.1233243) q[0];
sx q[0];
rz(2.841903) q[0];
x q[1];
rz(0.17941189) q[2];
sx q[2];
rz(-2.0127014) q[2];
sx q[2];
rz(0.68540689) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5089607) q[1];
sx q[1];
rz(-2.1034985) q[1];
sx q[1];
rz(-1.6091634) q[1];
rz(-1.0107993) q[3];
sx q[3];
rz(-2.3084497) q[3];
sx q[3];
rz(3.0411947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7361043) q[2];
sx q[2];
rz(-2.2443503) q[2];
sx q[2];
rz(-1.415095) q[2];
rz(1.803558) q[3];
sx q[3];
rz(-1.3058563) q[3];
sx q[3];
rz(3.0726748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5445484) q[0];
sx q[0];
rz(-0.91240779) q[0];
sx q[0];
rz(-1.9507116) q[0];
rz(-1.4471588) q[1];
sx q[1];
rz(-1.2971327) q[1];
sx q[1];
rz(-1.4318633) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3668684) q[0];
sx q[0];
rz(-0.69126832) q[0];
sx q[0];
rz(0.33095215) q[0];
x q[1];
rz(2.55673) q[2];
sx q[2];
rz(-1.9205127) q[2];
sx q[2];
rz(-0.51285267) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9976342) q[1];
sx q[1];
rz(-1.9737509) q[1];
sx q[1];
rz(2.2657299) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1095554) q[3];
sx q[3];
rz(-0.83870974) q[3];
sx q[3];
rz(0.17994954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.65149629) q[2];
sx q[2];
rz(-2.8054674) q[2];
sx q[2];
rz(0.87745848) q[2];
rz(-1.6089926) q[3];
sx q[3];
rz(-1.420615) q[3];
sx q[3];
rz(2.1144313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.0090016) q[0];
sx q[0];
rz(-1.4132211) q[0];
sx q[0];
rz(2.6268517) q[0];
rz(0.71350907) q[1];
sx q[1];
rz(-2.3206354) q[1];
sx q[1];
rz(-1.6101507) q[1];
rz(1.2825251) q[2];
sx q[2];
rz(-1.2119515) q[2];
sx q[2];
rz(-1.8744946) q[2];
rz(1.3446829) q[3];
sx q[3];
rz(-1.825708) q[3];
sx q[3];
rz(0.59151266) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
