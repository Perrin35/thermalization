OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.70541731) q[0];
sx q[0];
rz(-2.5751312) q[0];
sx q[0];
rz(2.9705272) q[0];
rz(1.45362) q[1];
sx q[1];
rz(-0.34314081) q[1];
sx q[1];
rz(1.8106102) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0199466) q[0];
sx q[0];
rz(-2.2197147) q[0];
sx q[0];
rz(1.6365504) q[0];
rz(-pi) q[1];
rz(-1.5491484) q[2];
sx q[2];
rz(-2.0865555) q[2];
sx q[2];
rz(1.6025008) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6273856) q[1];
sx q[1];
rz(-1.5537019) q[1];
sx q[1];
rz(-1.8219201) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64294502) q[3];
sx q[3];
rz(-1.7503947) q[3];
sx q[3];
rz(2.5760866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.77582899) q[2];
sx q[2];
rz(-0.77314955) q[2];
sx q[2];
rz(1.3295056) q[2];
rz(1.7261516) q[3];
sx q[3];
rz(-1.5751782) q[3];
sx q[3];
rz(-0.20234385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.1502894) q[0];
sx q[0];
rz(-2.5973899) q[0];
sx q[0];
rz(1.6888899) q[0];
rz(0.20092043) q[1];
sx q[1];
rz(-2.0566437) q[1];
sx q[1];
rz(-2.8413049) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18314221) q[0];
sx q[0];
rz(-2.169325) q[0];
sx q[0];
rz(1.4248225) q[0];
rz(0.47313182) q[2];
sx q[2];
rz(-2.5169249) q[2];
sx q[2];
rz(-1.2166789) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3631565) q[1];
sx q[1];
rz(-1.2769715) q[1];
sx q[1];
rz(-2.9086962) q[1];
rz(-pi) q[2];
rz(-2.2829451) q[3];
sx q[3];
rz(-0.24752366) q[3];
sx q[3];
rz(1.1502707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.369027) q[2];
sx q[2];
rz(-2.7894661) q[2];
sx q[2];
rz(2.1035813) q[2];
rz(0.84233061) q[3];
sx q[3];
rz(-1.7269644) q[3];
sx q[3];
rz(0.37030181) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5623986) q[0];
sx q[0];
rz(-2.6694522) q[0];
sx q[0];
rz(-2.753479) q[0];
rz(0.072470486) q[1];
sx q[1];
rz(-1.4265172) q[1];
sx q[1];
rz(0.31633502) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65773327) q[0];
sx q[0];
rz(-1.3257926) q[0];
sx q[0];
rz(-0.0018951456) q[0];
rz(1.146831) q[2];
sx q[2];
rz(-1.0397054) q[2];
sx q[2];
rz(0.46734992) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9005147) q[1];
sx q[1];
rz(-1.1763651) q[1];
sx q[1];
rz(-2.6532252) q[1];
rz(-0.86795904) q[3];
sx q[3];
rz(-2.7083525) q[3];
sx q[3];
rz(0.058335282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0324273) q[2];
sx q[2];
rz(-0.18523231) q[2];
sx q[2];
rz(0.16080984) q[2];
rz(0.12604776) q[3];
sx q[3];
rz(-1.7536609) q[3];
sx q[3];
rz(-2.1179312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2407103) q[0];
sx q[0];
rz(-2.4375589) q[0];
sx q[0];
rz(2.3098992) q[0];
rz(1.7968934) q[1];
sx q[1];
rz(-0.99935499) q[1];
sx q[1];
rz(1.5136738) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3095703) q[0];
sx q[0];
rz(-0.67414588) q[0];
sx q[0];
rz(3.0679697) q[0];
rz(-pi) q[1];
rz(3.0604826) q[2];
sx q[2];
rz(-2.1334279) q[2];
sx q[2];
rz(-0.62603355) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.44805376) q[1];
sx q[1];
rz(-1.5539907) q[1];
sx q[1];
rz(2.0190713) q[1];
rz(-pi) q[2];
rz(-0.45659153) q[3];
sx q[3];
rz(-2.3559104) q[3];
sx q[3];
rz(-2.8031138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4251129) q[2];
sx q[2];
rz(-1.0093062) q[2];
sx q[2];
rz(1.224219) q[2];
rz(-0.41695693) q[3];
sx q[3];
rz(-1.0427534) q[3];
sx q[3];
rz(-0.52880374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.058218) q[0];
sx q[0];
rz(-2.871802) q[0];
sx q[0];
rz(-1.303724) q[0];
rz(-2.6858792) q[1];
sx q[1];
rz(-2.8790751) q[1];
sx q[1];
rz(3.0900893) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52841016) q[0];
sx q[0];
rz(-0.15325704) q[0];
sx q[0];
rz(0.86092237) q[0];
rz(-pi) q[1];
rz(1.0837469) q[2];
sx q[2];
rz(-0.32099989) q[2];
sx q[2];
rz(0.1333065) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.81899535) q[1];
sx q[1];
rz(-1.7261191) q[1];
sx q[1];
rz(2.0218693) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73370917) q[3];
sx q[3];
rz(-1.6624311) q[3];
sx q[3];
rz(2.0165781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6872528) q[2];
sx q[2];
rz(-1.3663102) q[2];
sx q[2];
rz(0.59801897) q[2];
rz(-0.55001843) q[3];
sx q[3];
rz(-0.23967448) q[3];
sx q[3];
rz(-1.1365183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0705868) q[0];
sx q[0];
rz(-3.0650009) q[0];
sx q[0];
rz(2.9396074) q[0];
rz(2.1760991) q[1];
sx q[1];
rz(-2.1172724) q[1];
sx q[1];
rz(-3.0583256) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4250454) q[0];
sx q[0];
rz(-1.3971359) q[0];
sx q[0];
rz(2.9437149) q[0];
rz(-pi) q[1];
rz(-1.303057) q[2];
sx q[2];
rz(-1.8362852) q[2];
sx q[2];
rz(2.7051085) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6492918) q[1];
sx q[1];
rz(-1.9363083) q[1];
sx q[1];
rz(-0.11458061) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5842651) q[3];
sx q[3];
rz(-1.8263655) q[3];
sx q[3];
rz(2.717358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0319556) q[2];
sx q[2];
rz(-1.6836616) q[2];
sx q[2];
rz(1.3827682) q[2];
rz(-0.2494732) q[3];
sx q[3];
rz(-1.8981257) q[3];
sx q[3];
rz(-1.680254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4270585) q[0];
sx q[0];
rz(-2.336851) q[0];
sx q[0];
rz(0.89685857) q[0];
rz(-2.8915021) q[1];
sx q[1];
rz(-0.18083328) q[1];
sx q[1];
rz(2.0239963) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8539124) q[0];
sx q[0];
rz(-2.2932862) q[0];
sx q[0];
rz(1.7476837) q[0];
rz(-0.859185) q[2];
sx q[2];
rz(-1.4118328) q[2];
sx q[2];
rz(-0.32260103) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1114137) q[1];
sx q[1];
rz(-1.3707268) q[1];
sx q[1];
rz(-2.5740037) q[1];
rz(-pi) q[2];
rz(1.8459122) q[3];
sx q[3];
rz(-1.0529622) q[3];
sx q[3];
rz(2.8003611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0573132) q[2];
sx q[2];
rz(-1.3998312) q[2];
sx q[2];
rz(-0.21952595) q[2];
rz(-0.38671842) q[3];
sx q[3];
rz(-2.3733449) q[3];
sx q[3];
rz(-0.99532551) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0893843) q[0];
sx q[0];
rz(-1.5445671) q[0];
sx q[0];
rz(-2.9242933) q[0];
rz(-3.1106588) q[1];
sx q[1];
rz(-0.63806454) q[1];
sx q[1];
rz(-2.4826179) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9650759) q[0];
sx q[0];
rz(-0.79916164) q[0];
sx q[0];
rz(1.6060711) q[0];
x q[1];
rz(2.4023513) q[2];
sx q[2];
rz(-0.87202245) q[2];
sx q[2];
rz(3.0657363) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.79170376) q[1];
sx q[1];
rz(-2.3096497) q[1];
sx q[1];
rz(1.8658584) q[1];
x q[2];
rz(-1.6625255) q[3];
sx q[3];
rz(-2.496521) q[3];
sx q[3];
rz(2.0446387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7028246) q[2];
sx q[2];
rz(-1.3805026) q[2];
sx q[2];
rz(0.32021114) q[2];
rz(-1.6163588) q[3];
sx q[3];
rz(-1.1956918) q[3];
sx q[3];
rz(-1.2493791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3808688) q[0];
sx q[0];
rz(-2.9601233) q[0];
sx q[0];
rz(2.4587801) q[0];
rz(1.0702417) q[1];
sx q[1];
rz(-1.2425334) q[1];
sx q[1];
rz(-1.210093) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72737981) q[0];
sx q[0];
rz(-2.7276037) q[0];
sx q[0];
rz(-1.0340286) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2145043) q[2];
sx q[2];
rz(-0.63535832) q[2];
sx q[2];
rz(-2.2941342) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6283419) q[1];
sx q[1];
rz(-1.8645789) q[1];
sx q[1];
rz(2.252584) q[1];
x q[2];
rz(-1.9425415) q[3];
sx q[3];
rz(-2.1575655) q[3];
sx q[3];
rz(1.1128807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5658297) q[2];
sx q[2];
rz(-2.1255707) q[2];
sx q[2];
rz(-0.56662095) q[2];
rz(-2.2120655) q[3];
sx q[3];
rz(-2.1598787) q[3];
sx q[3];
rz(1.367759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41314769) q[0];
sx q[0];
rz(-2.8840273) q[0];
sx q[0];
rz(1.7096827) q[0];
rz(-0.52629772) q[1];
sx q[1];
rz(-2.6142575) q[1];
sx q[1];
rz(-0.79968232) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73197094) q[0];
sx q[0];
rz(-2.7831166) q[0];
sx q[0];
rz(1.8147857) q[0];
rz(-1.4893555) q[2];
sx q[2];
rz(-1.4105721) q[2];
sx q[2];
rz(2.3022431) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.848222) q[1];
sx q[1];
rz(-1.9105517) q[1];
sx q[1];
rz(1.0287702) q[1];
rz(1.4950072) q[3];
sx q[3];
rz(-1.4255376) q[3];
sx q[3];
rz(-2.4491572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2852823) q[2];
sx q[2];
rz(-2.6976863) q[2];
sx q[2];
rz(0.69520673) q[2];
rz(0.55784145) q[3];
sx q[3];
rz(-1.7772243) q[3];
sx q[3];
rz(0.69303304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0859062) q[0];
sx q[0];
rz(-1.1924556) q[0];
sx q[0];
rz(-2.6299155) q[0];
rz(1.862539) q[1];
sx q[1];
rz(-2.3976354) q[1];
sx q[1];
rz(2.4639113) q[1];
rz(-0.002481133) q[2];
sx q[2];
rz(-2.8056792) q[2];
sx q[2];
rz(0.16028595) q[2];
rz(-0.79808509) q[3];
sx q[3];
rz(-1.6860262) q[3];
sx q[3];
rz(-1.891914) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];