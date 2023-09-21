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
rz(-0.17106549) q[0];
rz(1.45362) q[1];
sx q[1];
rz(-0.34314081) q[1];
sx q[1];
rz(-1.3309825) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1216461) q[0];
sx q[0];
rz(-2.2197147) q[0];
sx q[0];
rz(-1.5050423) q[0];
rz(-pi) q[1];
rz(1.5491484) q[2];
sx q[2];
rz(-2.0865555) q[2];
sx q[2];
rz(-1.6025008) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.060974412) q[1];
sx q[1];
rz(-1.3197101) q[1];
sx q[1];
rz(3.1239448) q[1];
x q[2];
rz(0.64294502) q[3];
sx q[3];
rz(-1.7503947) q[3];
sx q[3];
rz(-0.56550607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77582899) q[2];
sx q[2];
rz(-2.3684431) q[2];
sx q[2];
rz(-1.8120871) q[2];
rz(1.7261516) q[3];
sx q[3];
rz(-1.5664145) q[3];
sx q[3];
rz(0.20234385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1502894) q[0];
sx q[0];
rz(-2.5973899) q[0];
sx q[0];
rz(-1.4527028) q[0];
rz(-0.20092043) q[1];
sx q[1];
rz(-2.0566437) q[1];
sx q[1];
rz(-0.30028775) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6712924) q[0];
sx q[0];
rz(-1.4503345) q[0];
sx q[0];
rz(0.60351535) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2533721) q[2];
sx q[2];
rz(-2.1183287) q[2];
sx q[2];
rz(-2.4878793) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3631565) q[1];
sx q[1];
rz(-1.8646212) q[1];
sx q[1];
rz(-2.9086962) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16365666) q[3];
sx q[3];
rz(-1.3842584) q[3];
sx q[3];
rz(-1.2638307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.369027) q[2];
sx q[2];
rz(-0.3521266) q[2];
sx q[2];
rz(-2.1035813) q[2];
rz(0.84233061) q[3];
sx q[3];
rz(-1.4146283) q[3];
sx q[3];
rz(-0.37030181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.5623986) q[0];
sx q[0];
rz(-2.6694522) q[0];
sx q[0];
rz(2.753479) q[0];
rz(-0.072470486) q[1];
sx q[1];
rz(-1.7150755) q[1];
sx q[1];
rz(0.31633502) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4916723) q[0];
sx q[0];
rz(-0.24501093) q[0];
sx q[0];
rz(-1.578376) q[0];
rz(1.146831) q[2];
sx q[2];
rz(-2.1018873) q[2];
sx q[2];
rz(2.6742427) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1856826) q[1];
sx q[1];
rz(-2.5240286) q[1];
sx q[1];
rz(-2.4159141) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8510677) q[3];
sx q[3];
rz(-1.2447262) q[3];
sx q[3];
rz(2.3323004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1091653) q[2];
sx q[2];
rz(-0.18523231) q[2];
sx q[2];
rz(-0.16080984) q[2];
rz(3.0155449) q[3];
sx q[3];
rz(-1.7536609) q[3];
sx q[3];
rz(2.1179312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9008824) q[0];
sx q[0];
rz(-2.4375589) q[0];
sx q[0];
rz(-2.3098992) q[0];
rz(1.7968934) q[1];
sx q[1];
rz(-0.99935499) q[1];
sx q[1];
rz(1.5136738) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8228089) q[0];
sx q[0];
rz(-1.6167287) q[0];
sx q[0];
rz(-0.6728234) q[0];
x q[1];
rz(0.081110031) q[2];
sx q[2];
rz(-2.1334279) q[2];
sx q[2];
rz(-2.5155591) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0107683) q[1];
sx q[1];
rz(-2.0190034) q[1];
sx q[1];
rz(-3.1229449) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6850011) q[3];
sx q[3];
rz(-2.3559104) q[3];
sx q[3];
rz(0.33847889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7164798) q[2];
sx q[2];
rz(-2.1322865) q[2];
sx q[2];
rz(-1.224219) q[2];
rz(-0.41695693) q[3];
sx q[3];
rz(-2.0988393) q[3];
sx q[3];
rz(-2.6127889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.058218) q[0];
sx q[0];
rz(-0.26979065) q[0];
sx q[0];
rz(1.303724) q[0];
rz(-2.6858792) q[1];
sx q[1];
rz(-2.8790751) q[1];
sx q[1];
rz(3.0900893) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52841016) q[0];
sx q[0];
rz(-2.9883356) q[0];
sx q[0];
rz(-2.2806703) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2850045) q[2];
sx q[2];
rz(-1.7190061) q[2];
sx q[2];
rz(0.9718026) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4645849) q[1];
sx q[1];
rz(-2.0160463) q[1];
sx q[1];
rz(-0.17226179) q[1];
rz(2.4078835) q[3];
sx q[3];
rz(-1.6624311) q[3];
sx q[3];
rz(-2.0165781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.45433989) q[2];
sx q[2];
rz(-1.3663102) q[2];
sx q[2];
rz(-0.59801897) q[2];
rz(0.55001843) q[3];
sx q[3];
rz(-0.23967448) q[3];
sx q[3];
rz(1.1365183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0710058) q[0];
sx q[0];
rz(-0.07659176) q[0];
sx q[0];
rz(-2.9396074) q[0];
rz(0.96549353) q[1];
sx q[1];
rz(-2.1172724) q[1];
sx q[1];
rz(3.0583256) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.252713) q[0];
sx q[0];
rz(-1.3759334) q[0];
sx q[0];
rz(-1.3937508) q[0];
rz(-pi) q[1];
rz(-1.8385356) q[2];
sx q[2];
rz(-1.8362852) q[2];
sx q[2];
rz(-2.7051085) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.80379936) q[1];
sx q[1];
rz(-0.38227907) q[1];
sx q[1];
rz(1.2804968) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2721328) q[3];
sx q[3];
rz(-2.1080058) q[3];
sx q[3];
rz(-0.99029535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0319556) q[2];
sx q[2];
rz(-1.6836616) q[2];
sx q[2];
rz(1.3827682) q[2];
rz(-2.8921195) q[3];
sx q[3];
rz(-1.243467) q[3];
sx q[3];
rz(1.4613387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4270585) q[0];
sx q[0];
rz(-0.80474168) q[0];
sx q[0];
rz(2.2447341) q[0];
rz(-0.25009051) q[1];
sx q[1];
rz(-0.18083328) q[1];
sx q[1];
rz(-2.0239963) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8539124) q[0];
sx q[0];
rz(-0.84830647) q[0];
sx q[0];
rz(-1.7476837) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2824077) q[2];
sx q[2];
rz(-1.7297598) q[2];
sx q[2];
rz(-0.32260103) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.84278216) q[1];
sx q[1];
rz(-0.59814765) q[1];
sx q[1];
rz(-2.7808955) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8459122) q[3];
sx q[3];
rz(-2.0886305) q[3];
sx q[3];
rz(0.3412316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0573132) q[2];
sx q[2];
rz(-1.3998312) q[2];
sx q[2];
rz(-0.21952595) q[2];
rz(-2.7548742) q[3];
sx q[3];
rz(-0.76824776) q[3];
sx q[3];
rz(2.1462671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0893843) q[0];
sx q[0];
rz(-1.5970255) q[0];
sx q[0];
rz(-2.9242933) q[0];
rz(0.030933881) q[1];
sx q[1];
rz(-2.5035281) q[1];
sx q[1];
rz(-0.6589748) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0156408) q[0];
sx q[0];
rz(-2.3693187) q[0];
sx q[0];
rz(3.1053567) q[0];
rz(-pi) q[1];
rz(-2.4023513) q[2];
sx q[2];
rz(-0.87202245) q[2];
sx q[2];
rz(-3.0657363) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.36775667) q[1];
sx q[1];
rz(-0.78513297) q[1];
sx q[1];
rz(2.8326041) q[1];
rz(-pi) q[2];
rz(-1.6625255) q[3];
sx q[3];
rz(-0.64507161) q[3];
sx q[3];
rz(1.0969539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7028246) q[2];
sx q[2];
rz(-1.3805026) q[2];
sx q[2];
rz(-0.32021114) q[2];
rz(1.6163588) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(1.8922136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76072389) q[0];
sx q[0];
rz(-2.9601233) q[0];
sx q[0];
rz(-0.6828126) q[0];
rz(-2.071351) q[1];
sx q[1];
rz(-1.8990592) q[1];
sx q[1];
rz(1.210093) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4142128) q[0];
sx q[0];
rz(-0.41398898) q[0];
sx q[0];
rz(-1.0340286) q[0];
rz(1.9270883) q[2];
sx q[2];
rz(-2.5062343) q[2];
sx q[2];
rz(2.2941342) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4267537) q[1];
sx q[1];
rz(-2.4086082) q[1];
sx q[1];
rz(1.1232125) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50001486) q[3];
sx q[3];
rz(-2.4588636) q[3];
sx q[3];
rz(-0.4993712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5658297) q[2];
sx q[2];
rz(-1.016022) q[2];
sx q[2];
rz(0.56662095) q[2];
rz(-0.9295272) q[3];
sx q[3];
rz(-0.98171392) q[3];
sx q[3];
rz(1.367759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.41314769) q[0];
sx q[0];
rz(-0.25756535) q[0];
sx q[0];
rz(-1.7096827) q[0];
rz(-2.6152949) q[1];
sx q[1];
rz(-2.6142575) q[1];
sx q[1];
rz(0.79968232) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.66946) q[0];
sx q[0];
rz(-1.9181983) q[0];
sx q[0];
rz(3.0513289) q[0];
x q[1];
rz(-0.16074796) q[2];
sx q[2];
rz(-1.6511917) q[2];
sx q[2];
rz(-0.71842566) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4754776) q[1];
sx q[1];
rz(-1.0628504) q[1];
sx q[1];
rz(-2.75027) q[1];
rz(-1.4950072) q[3];
sx q[3];
rz(-1.7160551) q[3];
sx q[3];
rz(-2.4491572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2852823) q[2];
sx q[2];
rz(-2.6976863) q[2];
sx q[2];
rz(2.4463859) q[2];
rz(-2.5837512) q[3];
sx q[3];
rz(-1.3643684) q[3];
sx q[3];
rz(-0.69303304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0859062) q[0];
sx q[0];
rz(-1.9491371) q[0];
sx q[0];
rz(0.51167713) q[0];
rz(1.862539) q[1];
sx q[1];
rz(-2.3976354) q[1];
sx q[1];
rz(2.4639113) q[1];
rz(-3.1391115) q[2];
sx q[2];
rz(-0.33591349) q[2];
sx q[2];
rz(-2.9813067) q[2];
rz(1.7351032) q[3];
sx q[3];
rz(-2.3621029) q[3];
sx q[3];
rz(2.7030871) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];