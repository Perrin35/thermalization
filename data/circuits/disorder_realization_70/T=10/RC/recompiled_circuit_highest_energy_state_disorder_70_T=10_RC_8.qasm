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
rz(-1.3938067) q[0];
sx q[0];
rz(-0.90016794) q[0];
sx q[0];
rz(-1.409344) q[0];
rz(-2.4318168) q[1];
sx q[1];
rz(-0.33325279) q[1];
sx q[1];
rz(2.7594653) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2173806) q[0];
sx q[0];
rz(-1.4367625) q[0];
sx q[0];
rz(-0.73705518) q[0];
x q[1];
rz(1.8252188) q[2];
sx q[2];
rz(-1.7013936) q[2];
sx q[2];
rz(-2.8859464) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14146261) q[1];
sx q[1];
rz(-1.9381818) q[1];
sx q[1];
rz(-2.7482583) q[1];
x q[2];
rz(2.2629991) q[3];
sx q[3];
rz(-1.4784364) q[3];
sx q[3];
rz(2.0372932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0673151) q[2];
sx q[2];
rz(-1.6396319) q[2];
sx q[2];
rz(0.074946694) q[2];
rz(1.2429169) q[3];
sx q[3];
rz(-0.26151812) q[3];
sx q[3];
rz(2.7459131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0089834) q[0];
sx q[0];
rz(-0.41985303) q[0];
sx q[0];
rz(-0.62927759) q[0];
rz(2.3043326) q[1];
sx q[1];
rz(-2.3910797) q[1];
sx q[1];
rz(-0.57964051) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64855498) q[0];
sx q[0];
rz(-0.84024197) q[0];
sx q[0];
rz(-1.8454946) q[0];
rz(-pi) q[1];
rz(0.39548611) q[2];
sx q[2];
rz(-1.8137534) q[2];
sx q[2];
rz(-2.0045351) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.77387735) q[1];
sx q[1];
rz(-1.9019777) q[1];
sx q[1];
rz(0.41172387) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9778538) q[3];
sx q[3];
rz(-1.5904015) q[3];
sx q[3];
rz(-0.99359182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2418182) q[2];
sx q[2];
rz(-2.9782229) q[2];
sx q[2];
rz(0.74431288) q[2];
rz(-2.7742079) q[3];
sx q[3];
rz(-1.8495879) q[3];
sx q[3];
rz(1.415409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.431417) q[0];
sx q[0];
rz(-0.95040584) q[0];
sx q[0];
rz(1.0071734) q[0];
rz(-3.1305283) q[1];
sx q[1];
rz(-2.8304351) q[1];
sx q[1];
rz(0.99753582) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38356203) q[0];
sx q[0];
rz(-0.13995384) q[0];
sx q[0];
rz(1.7566232) q[0];
rz(-2.3332525) q[2];
sx q[2];
rz(-1.6534717) q[2];
sx q[2];
rz(-2.9660564) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5845452) q[1];
sx q[1];
rz(-1.5633094) q[1];
sx q[1];
rz(-2.6232031) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5854534) q[3];
sx q[3];
rz(-2.08257) q[3];
sx q[3];
rz(1.1095593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.02115383) q[2];
sx q[2];
rz(-2.8595371) q[2];
sx q[2];
rz(0.8417449) q[2];
rz(0.65819955) q[3];
sx q[3];
rz(-2.2704312) q[3];
sx q[3];
rz(3.1200718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72092527) q[0];
sx q[0];
rz(-3.0166716) q[0];
sx q[0];
rz(2.4125873) q[0];
rz(0.76338243) q[1];
sx q[1];
rz(-0.63012505) q[1];
sx q[1];
rz(-2.7785832) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.027452) q[0];
sx q[0];
rz(-1.318256) q[0];
sx q[0];
rz(2.8414498) q[0];
rz(-0.98539008) q[2];
sx q[2];
rz(-2.3040758) q[2];
sx q[2];
rz(2.458141) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2263086) q[1];
sx q[1];
rz(-1.4292418) q[1];
sx q[1];
rz(2.841921) q[1];
rz(-pi) q[2];
rz(-1.6239802) q[3];
sx q[3];
rz(-1.6973572) q[3];
sx q[3];
rz(2.3317331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.18569854) q[2];
sx q[2];
rz(-0.82315767) q[2];
sx q[2];
rz(1.6419179) q[2];
rz(-0.58756346) q[3];
sx q[3];
rz(-2.1524119) q[3];
sx q[3];
rz(-0.51830083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8374306) q[0];
sx q[0];
rz(-0.70697933) q[0];
sx q[0];
rz(2.8564603) q[0];
rz(-0.25310165) q[1];
sx q[1];
rz(-2.1123835) q[1];
sx q[1];
rz(-2.0957799) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3110712) q[0];
sx q[0];
rz(-1.1806628) q[0];
sx q[0];
rz(-2.5907787) q[0];
x q[1];
rz(-1.990647) q[2];
sx q[2];
rz(-1.1641181) q[2];
sx q[2];
rz(2.5620154) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1872028) q[1];
sx q[1];
rz(-0.39785255) q[1];
sx q[1];
rz(-2.9955762) q[1];
rz(-pi) q[2];
x q[2];
rz(1.771403) q[3];
sx q[3];
rz(-0.90292519) q[3];
sx q[3];
rz(-2.843586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.73695856) q[2];
sx q[2];
rz(-1.0598695) q[2];
sx q[2];
rz(-0.57445478) q[2];
rz(-2.5308841) q[3];
sx q[3];
rz(-0.52053958) q[3];
sx q[3];
rz(-2.0846539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8059998) q[0];
sx q[0];
rz(-1.7570423) q[0];
sx q[0];
rz(-0.35032508) q[0];
rz(-2.8490745) q[1];
sx q[1];
rz(-0.12648335) q[1];
sx q[1];
rz(0.77936053) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8838046) q[0];
sx q[0];
rz(-1.5124) q[0];
sx q[0];
rz(1.6898872) q[0];
rz(-1.3931403) q[2];
sx q[2];
rz(-2.0442932) q[2];
sx q[2];
rz(0.23261468) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8694589) q[1];
sx q[1];
rz(-1.3682406) q[1];
sx q[1];
rz(0.32702469) q[1];
rz(-pi) q[2];
rz(2.8661714) q[3];
sx q[3];
rz(-0.53715992) q[3];
sx q[3];
rz(-2.8583683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.23955762) q[2];
sx q[2];
rz(-1.5152479) q[2];
sx q[2];
rz(2.1774192) q[2];
rz(0.17298175) q[3];
sx q[3];
rz(-1.0123342) q[3];
sx q[3];
rz(0.13121901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5445589) q[0];
sx q[0];
rz(-1.1976765) q[0];
sx q[0];
rz(3.0721989) q[0];
rz(1.3736877) q[1];
sx q[1];
rz(-1.2488139) q[1];
sx q[1];
rz(-0.51838851) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98250472) q[0];
sx q[0];
rz(-2.2347274) q[0];
sx q[0];
rz(2.9256665) q[0];
rz(2.3396427) q[2];
sx q[2];
rz(-0.63849245) q[2];
sx q[2];
rz(-1.0494159) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6024149) q[1];
sx q[1];
rz(-1.5449734) q[1];
sx q[1];
rz(-3.0033545) q[1];
rz(-pi) q[2];
rz(-0.33643548) q[3];
sx q[3];
rz(-2.9671768) q[3];
sx q[3];
rz(2.1788696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9531276) q[2];
sx q[2];
rz(-1.943482) q[2];
sx q[2];
rz(0.24448621) q[2];
rz(1.9237579) q[3];
sx q[3];
rz(-0.094450258) q[3];
sx q[3];
rz(-2.794246) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.090488) q[0];
sx q[0];
rz(-0.29260391) q[0];
sx q[0];
rz(-2.4421413) q[0];
rz(-0.72169101) q[1];
sx q[1];
rz(-1.3612008) q[1];
sx q[1];
rz(-1.1915421) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22868294) q[0];
sx q[0];
rz(-1.7828724) q[0];
sx q[0];
rz(-1.6985083) q[0];
rz(2.4999111) q[2];
sx q[2];
rz(-2.8770718) q[2];
sx q[2];
rz(-1.7873639) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4556132) q[1];
sx q[1];
rz(-1.0334823) q[1];
sx q[1];
rz(-0.098453589) q[1];
rz(0.64719836) q[3];
sx q[3];
rz(-2.5997926) q[3];
sx q[3];
rz(1.7447646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9222074) q[2];
sx q[2];
rz(-0.81874138) q[2];
sx q[2];
rz(1.7259664) q[2];
rz(2.7042232) q[3];
sx q[3];
rz(-0.24961095) q[3];
sx q[3];
rz(-2.6387446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9472083) q[0];
sx q[0];
rz(-1.9358862) q[0];
sx q[0];
rz(0.44596392) q[0];
rz(-1.3392316) q[1];
sx q[1];
rz(-2.4868592) q[1];
sx q[1];
rz(1.9816678) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9438749) q[0];
sx q[0];
rz(-0.76043429) q[0];
sx q[0];
rz(0.51087733) q[0];
rz(-pi) q[1];
rz(-1.8322938) q[2];
sx q[2];
rz(-1.6435677) q[2];
sx q[2];
rz(2.6852599) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7801108) q[1];
sx q[1];
rz(-2.105753) q[1];
sx q[1];
rz(-1.2992211) q[1];
rz(-pi) q[2];
rz(2.5612513) q[3];
sx q[3];
rz(-2.146462) q[3];
sx q[3];
rz(2.2656296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7356073) q[2];
sx q[2];
rz(-2.1129825) q[2];
sx q[2];
rz(-2.410991) q[2];
rz(-2.2405911) q[3];
sx q[3];
rz(-0.42947072) q[3];
sx q[3];
rz(3.1072531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63448298) q[0];
sx q[0];
rz(-2.6431838) q[0];
sx q[0];
rz(0.46257567) q[0];
rz(-1.0628465) q[1];
sx q[1];
rz(-1.7741508) q[1];
sx q[1];
rz(-0.06632334) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1076814) q[0];
sx q[0];
rz(-2.4723791) q[0];
sx q[0];
rz(0.50743703) q[0];
rz(-pi) q[1];
rz(1.9059876) q[2];
sx q[2];
rz(-1.4045713) q[2];
sx q[2];
rz(0.34632296) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6219225) q[1];
sx q[1];
rz(-2.1035476) q[1];
sx q[1];
rz(2.4581535) q[1];
x q[2];
rz(1.0995501) q[3];
sx q[3];
rz(-2.6726035) q[3];
sx q[3];
rz(-0.8595621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1041717) q[2];
sx q[2];
rz(-2.6080242) q[2];
sx q[2];
rz(1.2332234) q[2];
rz(-0.45796606) q[3];
sx q[3];
rz(-0.27126867) q[3];
sx q[3];
rz(0.39877322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0290699) q[0];
sx q[0];
rz(-1.7081013) q[0];
sx q[0];
rz(-1.3992455) q[0];
rz(-0.15432547) q[1];
sx q[1];
rz(-1.4230774) q[1];
sx q[1];
rz(1.9565061) q[1];
rz(2.828601) q[2];
sx q[2];
rz(-1.9421158) q[2];
sx q[2];
rz(-2.4696642) q[2];
rz(0.48404398) q[3];
sx q[3];
rz(-2.0917907) q[3];
sx q[3];
rz(1.2863822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
