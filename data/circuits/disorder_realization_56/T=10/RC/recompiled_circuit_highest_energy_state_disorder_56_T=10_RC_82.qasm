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
rz(2.7409878) q[0];
sx q[0];
rz(2.7317943) q[0];
sx q[0];
rz(8.3538342) q[0];
rz(0.61869705) q[1];
sx q[1];
rz(3.7096042) q[1];
sx q[1];
rz(8.5811442) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5045568) q[0];
sx q[0];
rz(-1.8352274) q[0];
sx q[0];
rz(0.030585551) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29542342) q[2];
sx q[2];
rz(-1.43338) q[2];
sx q[2];
rz(2.5538921) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.31402097) q[1];
sx q[1];
rz(-2.2306475) q[1];
sx q[1];
rz(0.41484264) q[1];
x q[2];
rz(2.6671404) q[3];
sx q[3];
rz(-2.0182924) q[3];
sx q[3];
rz(1.3407269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3598651) q[2];
sx q[2];
rz(-2.3024776) q[2];
sx q[2];
rz(0.016999379) q[2];
rz(-1.4012236) q[3];
sx q[3];
rz(-0.56178105) q[3];
sx q[3];
rz(-1.2167654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5725937) q[0];
sx q[0];
rz(-1.3324791) q[0];
sx q[0];
rz(-0.78729415) q[0];
rz(-2.9632118) q[1];
sx q[1];
rz(-1.9352103) q[1];
sx q[1];
rz(-1.23752) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2923861) q[0];
sx q[0];
rz(-0.9885177) q[0];
sx q[0];
rz(-1.4968064) q[0];
x q[1];
rz(-0.95713116) q[2];
sx q[2];
rz(-1.3295914) q[2];
sx q[2];
rz(-1.4184191) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5893277) q[1];
sx q[1];
rz(-2.5922212) q[1];
sx q[1];
rz(-2.2087847) q[1];
rz(-pi) q[2];
rz(-2.5982598) q[3];
sx q[3];
rz(-2.6381734) q[3];
sx q[3];
rz(-1.0375298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7201207) q[2];
sx q[2];
rz(-1.6933491) q[2];
sx q[2];
rz(0.063701542) q[2];
rz(1.2740159) q[3];
sx q[3];
rz(-2.2985022) q[3];
sx q[3];
rz(-1.5773076) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.098123) q[0];
sx q[0];
rz(-1.1935357) q[0];
sx q[0];
rz(-2.3980339) q[0];
rz(2.6877563) q[1];
sx q[1];
rz(-1.3498787) q[1];
sx q[1];
rz(0.41890621) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.009368816) q[0];
sx q[0];
rz(-2.0133349) q[0];
sx q[0];
rz(-2.924463) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9100283) q[2];
sx q[2];
rz(-0.76398173) q[2];
sx q[2];
rz(0.52896777) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.53330227) q[1];
sx q[1];
rz(-2.4455297) q[1];
sx q[1];
rz(1.1185557) q[1];
rz(-2.2499372) q[3];
sx q[3];
rz(-1.4640995) q[3];
sx q[3];
rz(-1.5859491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0539187) q[2];
sx q[2];
rz(-2.3232465) q[2];
sx q[2];
rz(1.9107001) q[2];
rz(-2.6026717) q[3];
sx q[3];
rz(-1.6659707) q[3];
sx q[3];
rz(-1.6787136) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26381275) q[0];
sx q[0];
rz(-0.75436622) q[0];
sx q[0];
rz(-2.5394649) q[0];
rz(1.6944132) q[1];
sx q[1];
rz(-2.4592631) q[1];
sx q[1];
rz(0.86300659) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2857638) q[0];
sx q[0];
rz(-1.6842972) q[0];
sx q[0];
rz(2.4185926) q[0];
rz(2.8406937) q[2];
sx q[2];
rz(-1.2105807) q[2];
sx q[2];
rz(0.27733251) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0095686) q[1];
sx q[1];
rz(-1.3088262) q[1];
sx q[1];
rz(-0.97673194) q[1];
x q[2];
rz(1.9163777) q[3];
sx q[3];
rz(-2.3732819) q[3];
sx q[3];
rz(-2.5490724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0756695) q[2];
sx q[2];
rz(-1.1752335) q[2];
sx q[2];
rz(-2.2389331) q[2];
rz(-2.3265808) q[3];
sx q[3];
rz(-0.60324001) q[3];
sx q[3];
rz(2.7284315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0478504) q[0];
sx q[0];
rz(-0.047690064) q[0];
sx q[0];
rz(-0.93511859) q[0];
rz(1.6085666) q[1];
sx q[1];
rz(-0.32564274) q[1];
sx q[1];
rz(0.26587048) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4971127) q[0];
sx q[0];
rz(-1.2605259) q[0];
sx q[0];
rz(-3.029986) q[0];
rz(-2.8433617) q[2];
sx q[2];
rz(-2.3687393) q[2];
sx q[2];
rz(1.7333584) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2468202) q[1];
sx q[1];
rz(-2.7547464) q[1];
sx q[1];
rz(-0.43718613) q[1];
rz(-pi) q[2];
rz(-1.3060027) q[3];
sx q[3];
rz(-1.6700498) q[3];
sx q[3];
rz(2.3276727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.461146) q[2];
sx q[2];
rz(-1.5004044) q[2];
sx q[2];
rz(-2.6908596) q[2];
rz(1.9510673) q[3];
sx q[3];
rz(-0.7014941) q[3];
sx q[3];
rz(-0.43959555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9615237) q[0];
sx q[0];
rz(-2.838205) q[0];
sx q[0];
rz(-3.0971089) q[0];
rz(-0.26875177) q[1];
sx q[1];
rz(-0.93695295) q[1];
sx q[1];
rz(-0.43089795) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.160153) q[0];
sx q[0];
rz(-2.4593076) q[0];
sx q[0];
rz(-0.42295608) q[0];
rz(-pi) q[1];
rz(2.1304279) q[2];
sx q[2];
rz(-0.75087386) q[2];
sx q[2];
rz(-2.4023285) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0468586) q[1];
sx q[1];
rz(-2.2206535) q[1];
sx q[1];
rz(-2.6296494) q[1];
rz(-pi) q[2];
rz(-0.66029064) q[3];
sx q[3];
rz(-1.5662282) q[3];
sx q[3];
rz(-0.82483236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0253133) q[2];
sx q[2];
rz(-2.7677324) q[2];
sx q[2];
rz(-2.5171793) q[2];
rz(-2.0404909) q[3];
sx q[3];
rz(-2.3291984) q[3];
sx q[3];
rz(-0.98744923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2295912) q[0];
sx q[0];
rz(-2.0822552) q[0];
sx q[0];
rz(-0.68369317) q[0];
rz(-1.4423485) q[1];
sx q[1];
rz(-0.78386274) q[1];
sx q[1];
rz(-0.20194617) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5324704) q[0];
sx q[0];
rz(-1.8989855) q[0];
sx q[0];
rz(-0.0013602982) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6228637) q[2];
sx q[2];
rz(-0.74801842) q[2];
sx q[2];
rz(1.5128653) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4836639) q[1];
sx q[1];
rz(-0.52573181) q[1];
sx q[1];
rz(0.13327285) q[1];
rz(0.45777623) q[3];
sx q[3];
rz(-1.5321863) q[3];
sx q[3];
rz(-0.82070551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.086229) q[2];
sx q[2];
rz(-1.7217041) q[2];
sx q[2];
rz(-1.281338) q[2];
rz(2.4751439) q[3];
sx q[3];
rz(-1.0969578) q[3];
sx q[3];
rz(-0.54653978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1619038) q[0];
sx q[0];
rz(-0.55210102) q[0];
sx q[0];
rz(-0.51396489) q[0];
rz(-0.95298302) q[1];
sx q[1];
rz(-2.516808) q[1];
sx q[1];
rz(2.0228588) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0054487451) q[0];
sx q[0];
rz(-1.6523673) q[0];
sx q[0];
rz(2.7177627) q[0];
rz(2.7114026) q[2];
sx q[2];
rz(-2.3025844) q[2];
sx q[2];
rz(2.6553287) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9671849) q[1];
sx q[1];
rz(-0.46665114) q[1];
sx q[1];
rz(1.0789167) q[1];
x q[2];
rz(-2.5932157) q[3];
sx q[3];
rz(-1.8781239) q[3];
sx q[3];
rz(-2.4509094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3861367) q[2];
sx q[2];
rz(-2.7894207) q[2];
sx q[2];
rz(-0.90292162) q[2];
rz(-1.4593982) q[3];
sx q[3];
rz(-1.3420339) q[3];
sx q[3];
rz(2.3193147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.4626386) q[0];
sx q[0];
rz(-2.8215388) q[0];
sx q[0];
rz(1.1902887) q[0];
rz(-2.8207488) q[1];
sx q[1];
rz(-1.4535934) q[1];
sx q[1];
rz(1.2215337) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0319043) q[0];
sx q[0];
rz(-1.3656989) q[0];
sx q[0];
rz(1.431365) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8655445) q[2];
sx q[2];
rz(-0.85535565) q[2];
sx q[2];
rz(-0.1705585) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1561191) q[1];
sx q[1];
rz(-1.3669984) q[1];
sx q[1];
rz(0.76038313) q[1];
x q[2];
rz(2.9469195) q[3];
sx q[3];
rz(-1.750573) q[3];
sx q[3];
rz(-2.468022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.90423501) q[2];
sx q[2];
rz(-0.15494896) q[2];
sx q[2];
rz(-0.3332738) q[2];
rz(-2.1244369) q[3];
sx q[3];
rz(-2.1867496) q[3];
sx q[3];
rz(-2.8216383) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3453813) q[0];
sx q[0];
rz(-2.6295202) q[0];
sx q[0];
rz(0.91878015) q[0];
rz(0.19110075) q[1];
sx q[1];
rz(-0.42593503) q[1];
sx q[1];
rz(0.79413116) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9730703) q[0];
sx q[0];
rz(-2.2380073) q[0];
sx q[0];
rz(2.4486827) q[0];
rz(-pi) q[1];
rz(-2.625611) q[2];
sx q[2];
rz(-0.40922874) q[2];
sx q[2];
rz(-0.87142631) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0865143) q[1];
sx q[1];
rz(-2.3740413) q[1];
sx q[1];
rz(-1.1805736) q[1];
rz(-pi) q[2];
rz(2.8201032) q[3];
sx q[3];
rz(-1.3703385) q[3];
sx q[3];
rz(0.050762477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0309151) q[2];
sx q[2];
rz(-2.5258749) q[2];
sx q[2];
rz(0.75882971) q[2];
rz(2.2714254) q[3];
sx q[3];
rz(-2.3535959) q[3];
sx q[3];
rz(-1.0677392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1398685) q[0];
sx q[0];
rz(-1.682946) q[0];
sx q[0];
rz(-1.2276822) q[0];
rz(0.43791804) q[1];
sx q[1];
rz(-1.4358078) q[1];
sx q[1];
rz(-1.5461071) q[1];
rz(-2.6545637) q[2];
sx q[2];
rz(-1.9430046) q[2];
sx q[2];
rz(-0.76443048) q[2];
rz(-2.1049166) q[3];
sx q[3];
rz(-0.86363367) q[3];
sx q[3];
rz(0.96994079) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
