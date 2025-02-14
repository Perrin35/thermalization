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
rz(-2.4988889) q[1];
sx q[1];
rz(0.42833498) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1550395) q[0];
sx q[0];
rz(-1.4615834) q[0];
sx q[0];
rz(0.10515736) q[0];
rz(0.36317715) q[2];
sx q[2];
rz(-2.4360949) q[2];
sx q[2];
rz(-0.35340912) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1312771) q[1];
sx q[1];
rz(-2.4020264) q[1];
sx q[1];
rz(-0.063994813) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.400757) q[3];
sx q[3];
rz(-0.29262283) q[3];
sx q[3];
rz(3.0081841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.20994818) q[2];
sx q[2];
rz(-0.95982176) q[2];
sx q[2];
rz(1.9317365) q[2];
rz(-3.0620388) q[3];
sx q[3];
rz(-0.39593655) q[3];
sx q[3];
rz(-2.615926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2570268) q[0];
sx q[0];
rz(-0.50978065) q[0];
sx q[0];
rz(0.38401815) q[0];
rz(-2.2513023) q[1];
sx q[1];
rz(-1.236981) q[1];
sx q[1];
rz(0.30337897) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3502894) q[0];
sx q[0];
rz(-1.3801738) q[0];
sx q[0];
rz(1.0269357) q[0];
rz(1.0531002) q[2];
sx q[2];
rz(-1.3441022) q[2];
sx q[2];
rz(-2.0402619) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5785256) q[1];
sx q[1];
rz(-1.4801229) q[1];
sx q[1];
rz(-3.0791326) q[1];
rz(2.5796765) q[3];
sx q[3];
rz(-0.73560916) q[3];
sx q[3];
rz(0.058678415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7552135) q[2];
sx q[2];
rz(-1.1073802) q[2];
sx q[2];
rz(-2.4066822) q[2];
rz(-2.2928061) q[3];
sx q[3];
rz(-0.74257344) q[3];
sx q[3];
rz(0.36062226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3291149) q[0];
sx q[0];
rz(-0.37549967) q[0];
sx q[0];
rz(0.078711674) q[0];
rz(0.82376897) q[1];
sx q[1];
rz(-1.0853826) q[1];
sx q[1];
rz(1.2944006) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58790197) q[0];
sx q[0];
rz(-1.6312092) q[0];
sx q[0];
rz(-1.6325476) q[0];
x q[1];
rz(-0.90581494) q[2];
sx q[2];
rz(-0.73630263) q[2];
sx q[2];
rz(-1.5018963) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11547281) q[1];
sx q[1];
rz(-1.2900262) q[1];
sx q[1];
rz(-1.4488279) q[1];
x q[2];
rz(2.9031676) q[3];
sx q[3];
rz(-2.2608071) q[3];
sx q[3];
rz(2.391614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2417629) q[2];
sx q[2];
rz(-0.37909847) q[2];
sx q[2];
rz(0.81643334) q[2];
rz(-1.5166616) q[3];
sx q[3];
rz(-2.1818325) q[3];
sx q[3];
rz(2.4412156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4888332) q[0];
sx q[0];
rz(-0.82010287) q[0];
sx q[0];
rz(0.56831992) q[0];
rz(1.142451) q[1];
sx q[1];
rz(-1.3520974) q[1];
sx q[1];
rz(-1.4521339) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88599151) q[0];
sx q[0];
rz(-1.7924249) q[0];
sx q[0];
rz(2.5021573) q[0];
rz(-pi) q[1];
rz(-2.9617519) q[2];
sx q[2];
rz(-1.6629863) q[2];
sx q[2];
rz(1.8123019) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6664847) q[1];
sx q[1];
rz(-2.1748073) q[1];
sx q[1];
rz(-0.38352769) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9098572) q[3];
sx q[3];
rz(-1.7393469) q[3];
sx q[3];
rz(-1.3685365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51231724) q[2];
sx q[2];
rz(-2.6111111) q[2];
sx q[2];
rz(-3.0647035) q[2];
rz(2.7422089) q[3];
sx q[3];
rz(-2.2279584) q[3];
sx q[3];
rz(0.54401773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15453108) q[0];
sx q[0];
rz(-0.35683826) q[0];
sx q[0];
rz(2.8416908) q[0];
rz(-1.9918282) q[1];
sx q[1];
rz(-1.0118142) q[1];
sx q[1];
rz(2.6531175) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73672527) q[0];
sx q[0];
rz(-0.18820764) q[0];
sx q[0];
rz(1.7467278) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59318329) q[2];
sx q[2];
rz(-1.1399684) q[2];
sx q[2];
rz(0.54774988) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8160745) q[1];
sx q[1];
rz(-1.5897004) q[1];
sx q[1];
rz(0.58362785) q[1];
rz(-pi) q[2];
rz(0.3533479) q[3];
sx q[3];
rz(-0.67253695) q[3];
sx q[3];
rz(1.8272682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38146314) q[2];
sx q[2];
rz(-3.0035786) q[2];
sx q[2];
rz(1.6628954) q[2];
rz(-1.1100769) q[3];
sx q[3];
rz(-0.9980945) q[3];
sx q[3];
rz(-2.4983675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7396624) q[0];
sx q[0];
rz(-1.1818385) q[0];
sx q[0];
rz(-2.8231743) q[0];
rz(0.034612522) q[1];
sx q[1];
rz(-2.5739659) q[1];
sx q[1];
rz(0.45077032) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0530659) q[0];
sx q[0];
rz(-2.6309816) q[0];
sx q[0];
rz(1.5998163) q[0];
x q[1];
rz(-0.87574739) q[2];
sx q[2];
rz(-2.108547) q[2];
sx q[2];
rz(1.328383) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.83799441) q[1];
sx q[1];
rz(-2.384409) q[1];
sx q[1];
rz(-0.96214575) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9616021) q[3];
sx q[3];
rz(-2.3473661) q[3];
sx q[3];
rz(-3.1193747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.047711756) q[2];
sx q[2];
rz(-0.1897976) q[2];
sx q[2];
rz(3.0799358) q[2];
rz(-3.0430074) q[3];
sx q[3];
rz(-0.74461377) q[3];
sx q[3];
rz(-2.4390167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.3848569) q[0];
sx q[0];
rz(-1.1133794) q[0];
sx q[0];
rz(0.039948832) q[0];
rz(-2.3710251) q[1];
sx q[1];
rz(-2.4295085) q[1];
sx q[1];
rz(0.22824731) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4994482) q[0];
sx q[0];
rz(-3.0239883) q[0];
sx q[0];
rz(2.2884877) q[0];
rz(-pi) q[1];
rz(3.1399916) q[2];
sx q[2];
rz(-1.2646156) q[2];
sx q[2];
rz(1.209335) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.7668263) q[1];
sx q[1];
rz(-0.09775459) q[1];
sx q[1];
rz(-0.575211) q[1];
rz(1.030121) q[3];
sx q[3];
rz(-0.70899963) q[3];
sx q[3];
rz(-2.0674771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0780636) q[2];
sx q[2];
rz(-1.0588131) q[2];
sx q[2];
rz(-0.25827363) q[2];
rz(0.20049788) q[3];
sx q[3];
rz(-2.343488) q[3];
sx q[3];
rz(-0.18331461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.9789326) q[0];
sx q[0];
rz(-1.7698092) q[0];
sx q[0];
rz(-0.3072511) q[0];
rz(-1.2112674) q[1];
sx q[1];
rz(-2.6478421) q[1];
sx q[1];
rz(2.5603851) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2659822) q[0];
sx q[0];
rz(-1.6020244) q[0];
sx q[0];
rz(-1.1660378) q[0];
rz(-pi) q[1];
rz(0.70971428) q[2];
sx q[2];
rz(-1.4462556) q[2];
sx q[2];
rz(-2.1309851) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.678486) q[1];
sx q[1];
rz(-2.3358279) q[1];
sx q[1];
rz(1.2918949) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5804306) q[3];
sx q[3];
rz(-2.5303839) q[3];
sx q[3];
rz(-1.8451913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4622978) q[2];
sx q[2];
rz(-1.1335979) q[2];
sx q[2];
rz(-3.1104258) q[2];
rz(-2.9160685) q[3];
sx q[3];
rz(-1.7799107) q[3];
sx q[3];
rz(-1.0303191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0533326) q[0];
sx q[0];
rz(-3.0971165) q[0];
sx q[0];
rz(-0.70892507) q[0];
rz(0.21253474) q[1];
sx q[1];
rz(-1.0147076) q[1];
sx q[1];
rz(0.85740352) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5564726) q[0];
sx q[0];
rz(-1.5758638) q[0];
sx q[0];
rz(1.4410254) q[0];
x q[1];
rz(1.8959664) q[2];
sx q[2];
rz(-2.0813) q[2];
sx q[2];
rz(2.4028962) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6826806) q[1];
sx q[1];
rz(-2.4574728) q[1];
sx q[1];
rz(-0.64351179) q[1];
x q[2];
rz(1.7836501) q[3];
sx q[3];
rz(-0.88485938) q[3];
sx q[3];
rz(2.6712772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7015486) q[2];
sx q[2];
rz(-0.35228071) q[2];
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
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.032967903) q[0];
sx q[0];
rz(-1.2297577) q[0];
sx q[0];
rz(-2.1627872) q[0];
rz(2.4188614) q[1];
sx q[1];
rz(-1.1284072) q[1];
sx q[1];
rz(0.61789787) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5586097) q[0];
sx q[0];
rz(-1.2961868) q[0];
sx q[0];
rz(-0.78297575) q[0];
rz(-1.0523791) q[2];
sx q[2];
rz(-0.978865) q[2];
sx q[2];
rz(-0.059378864) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.56983518) q[1];
sx q[1];
rz(-1.7392842) q[1];
sx q[1];
rz(-1.7525826) q[1];
rz(-pi) q[2];
x q[2];
rz(0.79094751) q[3];
sx q[3];
rz(-3.0229048) q[3];
sx q[3];
rz(3.0968248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.939398) q[2];
sx q[2];
rz(-1.7859744) q[2];
sx q[2];
rz(-0.00016577684) q[2];
rz(2.557142) q[3];
sx q[3];
rz(-2.1355459) q[3];
sx q[3];
rz(2.5463026) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38681876) q[0];
sx q[0];
rz(-1.5703572) q[0];
sx q[0];
rz(1.5686709) q[0];
rz(1.3407002) q[1];
sx q[1];
rz(-2.0410213) q[1];
sx q[1];
rz(-1.6165728) q[1];
rz(2.2095815) q[2];
sx q[2];
rz(-1.1975653) q[2];
sx q[2];
rz(-1.4690659) q[2];
rz(1.9289005) q[3];
sx q[3];
rz(-1.1546338) q[3];
sx q[3];
rz(3.0178937) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
