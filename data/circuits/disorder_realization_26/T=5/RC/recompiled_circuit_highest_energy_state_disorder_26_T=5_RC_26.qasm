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
rz(0.17231365) q[0];
sx q[0];
rz(3.3495164) q[0];
sx q[0];
rz(10.715545) q[0];
rz(0.15200226) q[1];
sx q[1];
rz(3.7842964) q[1];
sx q[1];
rz(9.8531129) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5458459) q[0];
sx q[0];
rz(-1.4662678) q[0];
sx q[0];
rz(-1.680611) q[0];
rz(-pi) q[1];
rz(2.7784155) q[2];
sx q[2];
rz(-2.4360949) q[2];
sx q[2];
rz(0.35340912) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.92378919) q[1];
sx q[1];
rz(-0.83309595) q[1];
sx q[1];
rz(1.6290725) q[1];
rz(-pi) q[2];
rz(0.74083565) q[3];
sx q[3];
rz(-0.29262283) q[3];
sx q[3];
rz(-0.13340852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9316445) q[2];
sx q[2];
rz(-0.95982176) q[2];
sx q[2];
rz(1.9317365) q[2];
rz(0.079553902) q[3];
sx q[3];
rz(-2.7456561) q[3];
sx q[3];
rz(-0.52566665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.2570268) q[0];
sx q[0];
rz(-0.50978065) q[0];
sx q[0];
rz(2.7575745) q[0];
rz(2.2513023) q[1];
sx q[1];
rz(-1.236981) q[1];
sx q[1];
rz(-0.30337897) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7913032) q[0];
sx q[0];
rz(-1.7614189) q[0];
sx q[0];
rz(-1.0269357) q[0];
x q[1];
rz(-2.0069507) q[2];
sx q[2];
rz(-2.580603) q[2];
sx q[2];
rz(-2.2962388) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0133923) q[1];
sx q[1];
rz(-1.6329995) q[1];
sx q[1];
rz(-1.661646) q[1];
x q[2];
rz(2.5796765) q[3];
sx q[3];
rz(-0.73560916) q[3];
sx q[3];
rz(0.058678415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.38637912) q[2];
sx q[2];
rz(-1.1073802) q[2];
sx q[2];
rz(0.73491043) q[2];
rz(0.84878659) q[3];
sx q[3];
rz(-0.74257344) q[3];
sx q[3];
rz(-2.7809704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3291149) q[0];
sx q[0];
rz(-2.766093) q[0];
sx q[0];
rz(0.078711674) q[0];
rz(-0.82376897) q[1];
sx q[1];
rz(-2.0562101) q[1];
sx q[1];
rz(-1.8471921) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5536907) q[0];
sx q[0];
rz(-1.6312092) q[0];
sx q[0];
rz(-1.6325476) q[0];
rz(-pi) q[1];
rz(-2.2357777) q[2];
sx q[2];
rz(-0.73630263) q[2];
sx q[2];
rz(1.5018963) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6523167) q[1];
sx q[1];
rz(-1.6879663) q[1];
sx q[1];
rz(0.28276171) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9031676) q[3];
sx q[3];
rz(-2.2608071) q[3];
sx q[3];
rz(0.7499786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.89982975) q[2];
sx q[2];
rz(-2.7624942) q[2];
sx q[2];
rz(-0.81643334) q[2];
rz(-1.624931) q[3];
sx q[3];
rz(-2.1818325) q[3];
sx q[3];
rz(-2.4412156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65275943) q[0];
sx q[0];
rz(-0.82010287) q[0];
sx q[0];
rz(0.56831992) q[0];
rz(-1.142451) q[1];
sx q[1];
rz(-1.3520974) q[1];
sx q[1];
rz(-1.6894587) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88599151) q[0];
sx q[0];
rz(-1.3491677) q[0];
sx q[0];
rz(-0.63943531) q[0];
x q[1];
rz(0.4770437) q[2];
sx q[2];
rz(-2.9397268) q[2];
sx q[2];
rz(-0.22722009) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4751079) q[1];
sx q[1];
rz(-2.1748073) q[1];
sx q[1];
rz(2.758065) q[1];
rz(2.9630757) q[3];
sx q[3];
rz(-1.2367289) q[3];
sx q[3];
rz(2.8802383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6292754) q[2];
sx q[2];
rz(-0.53048152) q[2];
sx q[2];
rz(3.0647035) q[2];
rz(-2.7422089) q[3];
sx q[3];
rz(-2.2279584) q[3];
sx q[3];
rz(-0.54401773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9870616) q[0];
sx q[0];
rz(-0.35683826) q[0];
sx q[0];
rz(2.8416908) q[0];
rz(-1.9918282) q[1];
sx q[1];
rz(-1.0118142) q[1];
sx q[1];
rz(2.6531175) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91575275) q[0];
sx q[0];
rz(-1.7560648) q[0];
sx q[0];
rz(-3.1082694) q[0];
rz(2.5484094) q[2];
sx q[2];
rz(-1.1399684) q[2];
sx q[2];
rz(-0.54774988) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8676973) q[1];
sx q[1];
rz(-0.58389837) q[1];
sx q[1];
rz(-0.03429596) q[1];
rz(-pi) q[2];
rz(1.3018872) q[3];
sx q[3];
rz(-0.9465512) q[3];
sx q[3];
rz(1.3866803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7601295) q[2];
sx q[2];
rz(-3.0035786) q[2];
sx q[2];
rz(1.4786973) q[2];
rz(-1.1100769) q[3];
sx q[3];
rz(-2.1434982) q[3];
sx q[3];
rz(-0.64322513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4019302) q[0];
sx q[0];
rz(-1.1818385) q[0];
sx q[0];
rz(-0.31841835) q[0];
rz(-0.034612522) q[1];
sx q[1];
rz(-2.5739659) q[1];
sx q[1];
rz(2.6908223) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0552669) q[0];
sx q[0];
rz(-2.0811715) q[0];
sx q[0];
rz(-0.016252131) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66025205) q[2];
sx q[2];
rz(-2.1532032) q[2];
sx q[2];
rz(-2.9803515) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.602421) q[1];
sx q[1];
rz(-0.97214593) q[1];
sx q[1];
rz(-0.4954229) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9616021) q[3];
sx q[3];
rz(-2.3473661) q[3];
sx q[3];
rz(-0.022217928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.047711756) q[2];
sx q[2];
rz(-0.1897976) q[2];
sx q[2];
rz(0.06165687) q[2];
rz(-0.098585248) q[3];
sx q[3];
rz(-2.3969789) q[3];
sx q[3];
rz(-2.4390167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75673574) q[0];
sx q[0];
rz(-1.1133794) q[0];
sx q[0];
rz(3.1016438) q[0];
rz(-2.3710251) q[1];
sx q[1];
rz(-0.71208411) q[1];
sx q[1];
rz(-0.22824731) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21439274) q[0];
sx q[0];
rz(-1.6480371) q[0];
sx q[0];
rz(-1.4820251) q[0];
x q[1];
rz(1.575861) q[2];
sx q[2];
rz(-0.3061848) q[2];
sx q[2];
rz(-1.926946) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.7668263) q[1];
sx q[1];
rz(-3.0438381) q[1];
sx q[1];
rz(-2.5663816) q[1];
rz(-pi) q[2];
rz(-2.1114717) q[3];
sx q[3];
rz(-2.432593) q[3];
sx q[3];
rz(-1.0741155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0780636) q[2];
sx q[2];
rz(-2.0827796) q[2];
sx q[2];
rz(-0.25827363) q[2];
rz(-0.20049788) q[3];
sx q[3];
rz(-2.343488) q[3];
sx q[3];
rz(0.18331461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9789326) q[0];
sx q[0];
rz(-1.3717835) q[0];
sx q[0];
rz(0.3072511) q[0];
rz(1.2112674) q[1];
sx q[1];
rz(-0.4937506) q[1];
sx q[1];
rz(2.5603851) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6224339) q[0];
sx q[0];
rz(-0.40589505) q[0];
sx q[0];
rz(-1.4916363) q[0];
rz(-pi) q[1];
rz(1.7343592) q[2];
sx q[2];
rz(-2.2738761) q[2];
sx q[2];
rz(-0.66649461) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0710268) q[1];
sx q[1];
rz(-0.80447996) q[1];
sx q[1];
rz(2.8623373) q[1];
rz(-2.5804306) q[3];
sx q[3];
rz(-2.5303839) q[3];
sx q[3];
rz(-1.8451913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4622978) q[2];
sx q[2];
rz(-2.0079948) q[2];
sx q[2];
rz(-3.1104258) q[2];
rz(-2.9160685) q[3];
sx q[3];
rz(-1.7799107) q[3];
sx q[3];
rz(2.1112736) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0533326) q[0];
sx q[0];
rz(-0.044476155) q[0];
sx q[0];
rz(0.70892507) q[0];
rz(-2.9290579) q[1];
sx q[1];
rz(-2.1268851) q[1];
sx q[1];
rz(-0.85740352) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014985059) q[0];
sx q[0];
rz(-1.7005655) q[0];
sx q[0];
rz(0.0051104498) q[0];
x q[1];
rz(0.53376824) q[2];
sx q[2];
rz(-1.2883175) q[2];
sx q[2];
rz(0.9953645) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.58536968) q[1];
sx q[1];
rz(-1.1818647) q[1];
sx q[1];
rz(-0.57805581) q[1];
rz(0.69714947) q[3];
sx q[3];
rz(-1.7350041) q[3];
sx q[3];
rz(-2.177161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7015486) q[2];
sx q[2];
rz(-2.7893119) q[2];
sx q[2];
rz(-2.3360543) q[2];
rz(2.7696179) q[3];
sx q[3];
rz(-1.6476846) q[3];
sx q[3];
rz(-0.39723799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1086248) q[0];
sx q[0];
rz(-1.2297577) q[0];
sx q[0];
rz(-2.1627872) q[0];
rz(0.72273123) q[1];
sx q[1];
rz(-1.1284072) q[1];
sx q[1];
rz(-0.61789787) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.582983) q[0];
sx q[0];
rz(-1.2961868) q[0];
sx q[0];
rz(-2.3586169) q[0];
x q[1];
rz(0.65872569) q[2];
sx q[2];
rz(-1.1470231) q[2];
sx q[2];
rz(1.9383333) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4010716) q[1];
sx q[1];
rz(-0.24722543) q[1];
sx q[1];
rz(0.81584357) q[1];
x q[2];
rz(-1.6553819) q[3];
sx q[3];
rz(-1.654155) q[3];
sx q[3];
rz(0.74970923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2021947) q[2];
sx q[2];
rz(-1.3556182) q[2];
sx q[2];
rz(3.1414269) q[2];
rz(2.557142) q[3];
sx q[3];
rz(-1.0060468) q[3];
sx q[3];
rz(-2.5463026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38681876) q[0];
sx q[0];
rz(-1.5712354) q[0];
sx q[0];
rz(-1.5729217) q[0];
rz(-1.8008925) q[1];
sx q[1];
rz(-2.0410213) q[1];
sx q[1];
rz(-1.6165728) q[1];
rz(-0.93201119) q[2];
sx q[2];
rz(-1.1975653) q[2];
sx q[2];
rz(-1.4690659) q[2];
rz(-1.9289005) q[3];
sx q[3];
rz(-1.9869589) q[3];
sx q[3];
rz(-0.12369894) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
