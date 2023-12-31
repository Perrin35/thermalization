OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4361753) q[0];
sx q[0];
rz(-0.56646148) q[0];
sx q[0];
rz(0.17106549) q[0];
rz(1.45362) q[1];
sx q[1];
rz(-0.34314081) q[1];
sx q[1];
rz(1.8106102) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91141191) q[0];
sx q[0];
rz(-0.6517621) q[0];
sx q[0];
rz(3.0551811) q[0];
rz(0.51585977) q[2];
sx q[2];
rz(-1.5519648) q[2];
sx q[2];
rz(-3.0992103) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1316489) q[1];
sx q[1];
rz(-2.8899) q[1];
sx q[1];
rz(-1.5021055) q[1];
rz(2.8475548) q[3];
sx q[3];
rz(-0.66411823) q[3];
sx q[3];
rz(0.77120632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.77582899) q[2];
sx q[2];
rz(-2.3684431) q[2];
sx q[2];
rz(1.3295056) q[2];
rz(-1.4154411) q[3];
sx q[3];
rz(-1.5751782) q[3];
sx q[3];
rz(2.9392488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99130327) q[0];
sx q[0];
rz(-2.5973899) q[0];
sx q[0];
rz(-1.4527028) q[0];
rz(0.20092043) q[1];
sx q[1];
rz(-2.0566437) q[1];
sx q[1];
rz(0.30028775) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9584504) q[0];
sx q[0];
rz(-2.169325) q[0];
sx q[0];
rz(1.4248225) q[0];
x q[1];
rz(2.5710201) q[2];
sx q[2];
rz(-1.3010446) q[2];
sx q[2];
rz(-2.3938993) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.86094942) q[1];
sx q[1];
rz(-1.3480535) q[1];
sx q[1];
rz(1.8722948) q[1];
rz(-2.2829451) q[3];
sx q[3];
rz(-0.24752366) q[3];
sx q[3];
rz(1.1502707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.77256569) q[2];
sx q[2];
rz(-2.7894661) q[2];
sx q[2];
rz(-1.0380113) q[2];
rz(-2.299262) q[3];
sx q[3];
rz(-1.7269644) q[3];
sx q[3];
rz(-2.7712908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-1.5623986) q[0];
sx q[0];
rz(-2.6694522) q[0];
sx q[0];
rz(2.753479) q[0];
rz(3.0691222) q[1];
sx q[1];
rz(-1.4265172) q[1];
sx q[1];
rz(-0.31633502) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2280699) q[0];
sx q[0];
rz(-1.5689578) q[0];
sx q[0];
rz(1.8158005) q[0];
rz(-pi) q[1];
rz(0.57245589) q[2];
sx q[2];
rz(-1.2081895) q[2];
sx q[2];
rz(2.2628757) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12831941) q[1];
sx q[1];
rz(-1.1228021) q[1];
sx q[1];
rz(-1.1303348) q[1];
rz(-2.2736336) q[3];
sx q[3];
rz(-0.43324019) q[3];
sx q[3];
rz(-3.0832574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1091653) q[2];
sx q[2];
rz(-0.18523231) q[2];
sx q[2];
rz(-2.9807828) q[2];
rz(-3.0155449) q[3];
sx q[3];
rz(-1.7536609) q[3];
sx q[3];
rz(-2.1179312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9008824) q[0];
sx q[0];
rz(-2.4375589) q[0];
sx q[0];
rz(-2.3098992) q[0];
rz(-1.7968934) q[1];
sx q[1];
rz(-0.99935499) q[1];
sx q[1];
rz(-1.5136738) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2154402) q[0];
sx q[0];
rz(-0.8988131) q[0];
sx q[0];
rz(-1.5120904) q[0];
rz(-0.081110031) q[2];
sx q[2];
rz(-1.0081648) q[2];
sx q[2];
rz(0.62603355) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.44805376) q[1];
sx q[1];
rz(-1.5876019) q[1];
sx q[1];
rz(1.1225213) q[1];
rz(-pi) q[2];
rz(1.1553331) q[3];
sx q[3];
rz(-2.2586125) q[3];
sx q[3];
rz(0.26879877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4251129) q[2];
sx q[2];
rz(-2.1322865) q[2];
sx q[2];
rz(1.9173737) q[2];
rz(0.41695693) q[3];
sx q[3];
rz(-2.0988393) q[3];
sx q[3];
rz(-0.52880374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-0.2625176) q[1];
sx q[1];
rz(-3.0900893) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9542959) q[0];
sx q[0];
rz(-1.4547537) q[0];
sx q[0];
rz(0.10033484) q[0];
x q[1];
rz(-1.2850045) q[2];
sx q[2];
rz(-1.7190061) q[2];
sx q[2];
rz(-0.9718026) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0806549) q[1];
sx q[1];
rz(-2.6662711) q[1];
sx q[1];
rz(1.9156384) q[1];
x q[2];
rz(0.73370917) q[3];
sx q[3];
rz(-1.6624311) q[3];
sx q[3];
rz(2.0165781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45433989) q[2];
sx q[2];
rz(-1.7752825) q[2];
sx q[2];
rz(-0.59801897) q[2];
rz(-0.55001843) q[3];
sx q[3];
rz(-0.23967448) q[3];
sx q[3];
rz(2.0050744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0710058) q[0];
sx q[0];
rz(-3.0650009) q[0];
sx q[0];
rz(-2.9396074) q[0];
rz(2.1760991) q[1];
sx q[1];
rz(-1.0243203) q[1];
sx q[1];
rz(-0.083267033) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88887963) q[0];
sx q[0];
rz(-1.3759334) q[0];
sx q[0];
rz(-1.3937508) q[0];
rz(-0.7716878) q[2];
sx q[2];
rz(-0.37479127) q[2];
sx q[2];
rz(-0.37116606) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0373842) q[1];
sx q[1];
rz(-1.4638149) q[1];
sx q[1];
rz(-1.9385098) q[1];
rz(-pi) q[2];
rz(-1.2721328) q[3];
sx q[3];
rz(-1.0335869) q[3];
sx q[3];
rz(-0.99029535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.10963708) q[2];
sx q[2];
rz(-1.457931) q[2];
sx q[2];
rz(-1.7588245) q[2];
rz(0.2494732) q[3];
sx q[3];
rz(-1.8981257) q[3];
sx q[3];
rz(-1.4613387) q[3];
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
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7145342) q[0];
sx q[0];
rz(-0.80474168) q[0];
sx q[0];
rz(-2.2447341) q[0];
rz(-2.8915021) q[1];
sx q[1];
rz(-2.9607594) q[1];
sx q[1];
rz(1.1175964) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9761336) q[0];
sx q[0];
rz(-1.7031859) q[0];
sx q[0];
rz(0.73029851) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3300702) q[2];
sx q[2];
rz(-0.72611085) q[2];
sx q[2];
rz(1.7118529) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.41457957) q[1];
sx q[1];
rz(-2.1257183) q[1];
sx q[1];
rz(1.3347866) q[1];
rz(2.696633) q[3];
sx q[3];
rz(-0.58044725) q[3];
sx q[3];
rz(2.2821033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0573132) q[2];
sx q[2];
rz(-1.3998312) q[2];
sx q[2];
rz(0.21952595) q[2];
rz(2.7548742) q[3];
sx q[3];
rz(-0.76824776) q[3];
sx q[3];
rz(0.99532551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052208386) q[0];
sx q[0];
rz(-1.5970255) q[0];
sx q[0];
rz(-2.9242933) q[0];
rz(3.1106588) q[1];
sx q[1];
rz(-2.5035281) q[1];
sx q[1];
rz(0.6589748) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4188822) q[0];
sx q[0];
rz(-1.5455149) q[0];
sx q[0];
rz(2.3696469) q[0];
rz(0.72139229) q[2];
sx q[2];
rz(-2.112769) q[2];
sx q[2];
rz(-1.1162356) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.79170376) q[1];
sx q[1];
rz(-2.3096497) q[1];
sx q[1];
rz(1.2757343) q[1];
x q[2];
rz(-1.6625255) q[3];
sx q[3];
rz(-2.496521) q[3];
sx q[3];
rz(-1.0969539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4387681) q[2];
sx q[2];
rz(-1.7610901) q[2];
sx q[2];
rz(-0.32021114) q[2];
rz(-1.6163588) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(-1.8922136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3808688) q[0];
sx q[0];
rz(-0.18146935) q[0];
sx q[0];
rz(-2.4587801) q[0];
rz(-1.0702417) q[1];
sx q[1];
rz(-1.8990592) q[1];
sx q[1];
rz(-1.210093) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72737981) q[0];
sx q[0];
rz(-0.41398898) q[0];
sx q[0];
rz(-1.0340286) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2145043) q[2];
sx q[2];
rz(-2.5062343) q[2];
sx q[2];
rz(-2.2941342) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6283419) q[1];
sx q[1];
rz(-1.2770137) q[1];
sx q[1];
rz(-2.252584) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50001486) q[3];
sx q[3];
rz(-2.4588636) q[3];
sx q[3];
rz(2.6422215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5658297) q[2];
sx q[2];
rz(-2.1255707) q[2];
sx q[2];
rz(0.56662095) q[2];
rz(-0.9295272) q[3];
sx q[3];
rz(-2.1598787) q[3];
sx q[3];
rz(-1.367759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41314769) q[0];
sx q[0];
rz(-0.25756535) q[0];
sx q[0];
rz(1.7096827) q[0];
rz(-0.52629772) q[1];
sx q[1];
rz(-0.52733517) q[1];
sx q[1];
rz(-2.3419103) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47213263) q[0];
sx q[0];
rz(-1.2233943) q[0];
sx q[0];
rz(3.0513289) q[0];
rz(-pi) q[1];
rz(-1.4893555) q[2];
sx q[2];
rz(-1.4105721) q[2];
sx q[2];
rz(2.3022431) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.29337063) q[1];
sx q[1];
rz(-1.231041) q[1];
sx q[1];
rz(1.0287702) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47761376) q[3];
sx q[3];
rz(-2.9778746) q[3];
sx q[3];
rz(1.9660266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2852823) q[2];
sx q[2];
rz(-2.6976863) q[2];
sx q[2];
rz(0.69520673) q[2];
rz(-0.55784145) q[3];
sx q[3];
rz(-1.7772243) q[3];
sx q[3];
rz(-0.69303304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0556864) q[0];
sx q[0];
rz(-1.9491371) q[0];
sx q[0];
rz(0.51167713) q[0];
rz(1.862539) q[1];
sx q[1];
rz(-2.3976354) q[1];
sx q[1];
rz(2.4639113) q[1];
rz(-2.8056801) q[2];
sx q[2];
rz(-1.5699785) q[2];
sx q[2];
rz(1.7287398) q[2];
rz(-2.3435076) q[3];
sx q[3];
rz(-1.4555664) q[3];
sx q[3];
rz(1.2496787) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
