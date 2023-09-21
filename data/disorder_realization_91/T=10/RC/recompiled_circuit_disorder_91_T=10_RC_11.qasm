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
rz(-1.3309825) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0199466) q[0];
sx q[0];
rz(-2.2197147) q[0];
sx q[0];
rz(1.5050423) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1034307) q[2];
sx q[2];
rz(-2.6254203) q[2];
sx q[2];
rz(1.49522) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.0099437873) q[1];
sx q[1];
rz(-0.25169262) q[1];
sx q[1];
rz(-1.6394872) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29403789) q[3];
sx q[3];
rz(-2.4774744) q[3];
sx q[3];
rz(-0.77120632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3657637) q[2];
sx q[2];
rz(-0.77314955) q[2];
sx q[2];
rz(1.3295056) q[2];
rz(1.4154411) q[3];
sx q[3];
rz(-1.5751782) q[3];
sx q[3];
rz(0.20234385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.8413049) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0694885) q[0];
sx q[0];
rz(-0.61394962) q[0];
sx q[0];
rz(2.9314562) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8882206) q[2];
sx q[2];
rz(-2.1183287) q[2];
sx q[2];
rz(-2.4878793) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.7784361) q[1];
sx q[1];
rz(-1.8646212) q[1];
sx q[1];
rz(-0.2328965) q[1];
x q[2];
rz(-2.2829451) q[3];
sx q[3];
rz(-2.894069) q[3];
sx q[3];
rz(1.991322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.369027) q[2];
sx q[2];
rz(-0.3521266) q[2];
sx q[2];
rz(2.1035813) q[2];
rz(2.299262) q[3];
sx q[3];
rz(-1.4146283) q[3];
sx q[3];
rz(-2.7712908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5623986) q[0];
sx q[0];
rz(-0.47214046) q[0];
sx q[0];
rz(2.753479) q[0];
rz(0.072470486) q[1];
sx q[1];
rz(-1.7150755) q[1];
sx q[1];
rz(-0.31633502) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2280699) q[0];
sx q[0];
rz(-1.5726349) q[0];
sx q[0];
rz(-1.8158005) q[0];
rz(-pi) q[1];
rz(-1.9947617) q[2];
sx q[2];
rz(-2.1018873) q[2];
sx q[2];
rz(-0.46734992) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2410779) q[1];
sx q[1];
rz(-1.1763651) q[1];
sx q[1];
rz(0.48836744) q[1];
rz(-pi) q[2];
rz(-1.9100788) q[3];
sx q[3];
rz(-1.8456036) q[3];
sx q[3];
rz(0.85698444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1091653) q[2];
sx q[2];
rz(-2.9563603) q[2];
sx q[2];
rz(2.9807828) q[2];
rz(-0.12604776) q[3];
sx q[3];
rz(-1.3879317) q[3];
sx q[3];
rz(-2.1179312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2407103) q[0];
sx q[0];
rz(-0.70403376) q[0];
sx q[0];
rz(-2.3098992) q[0];
rz(1.7968934) q[1];
sx q[1];
rz(-0.99935499) q[1];
sx q[1];
rz(-1.6279189) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2154402) q[0];
sx q[0];
rz(-0.8988131) q[0];
sx q[0];
rz(1.5120904) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0604826) q[2];
sx q[2];
rz(-2.1334279) q[2];
sx q[2];
rz(-0.62603355) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0537783) q[1];
sx q[1];
rz(-0.44856854) q[1];
sx q[1];
rz(1.5320369) q[1];
rz(1.9862595) q[3];
sx q[3];
rz(-2.2586125) q[3];
sx q[3];
rz(2.8727939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4251129) q[2];
sx q[2];
rz(-1.0093062) q[2];
sx q[2];
rz(-1.9173737) q[2];
rz(2.7246357) q[3];
sx q[3];
rz(-1.0427534) q[3];
sx q[3];
rz(-0.52880374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.083374627) q[0];
sx q[0];
rz(-0.26979065) q[0];
sx q[0];
rz(1.303724) q[0];
rz(-2.6858792) q[1];
sx q[1];
rz(-2.8790751) q[1];
sx q[1];
rz(-0.051503332) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3951552) q[0];
sx q[0];
rz(-1.4711385) q[0];
sx q[0];
rz(-1.4541724) q[0];
rz(-1.8565882) q[2];
sx q[2];
rz(-1.4225866) q[2];
sx q[2];
rz(-0.9718026) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.67700779) q[1];
sx q[1];
rz(-1.1255463) q[1];
sx q[1];
rz(-2.9693309) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4476942) q[3];
sx q[3];
rz(-0.8408635) q[3];
sx q[3];
rz(2.6134932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.45433989) q[2];
sx q[2];
rz(-1.7752825) q[2];
sx q[2];
rz(-2.5435737) q[2];
rz(-2.5915742) q[3];
sx q[3];
rz(-2.9019182) q[3];
sx q[3];
rz(-1.1365183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0705868) q[0];
sx q[0];
rz(-3.0650009) q[0];
sx q[0];
rz(0.20198527) q[0];
rz(-0.96549353) q[1];
sx q[1];
rz(-2.1172724) q[1];
sx q[1];
rz(0.083267033) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.252713) q[0];
sx q[0];
rz(-1.7656592) q[0];
sx q[0];
rz(1.7478419) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27481885) q[2];
sx q[2];
rz(-1.3126557) q[2];
sx q[2];
rz(1.2061662) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.80379936) q[1];
sx q[1];
rz(-2.7593136) q[1];
sx q[1];
rz(1.2804968) q[1];
rz(0.45883026) q[3];
sx q[3];
rz(-2.5341431) q[3];
sx q[3];
rz(1.6096887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
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
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7145342) q[0];
sx q[0];
rz(-2.336851) q[0];
sx q[0];
rz(-0.89685857) q[0];
rz(-0.25009051) q[1];
sx q[1];
rz(-0.18083328) q[1];
sx q[1];
rz(1.1175964) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9761336) q[0];
sx q[0];
rz(-1.7031859) q[0];
sx q[0];
rz(2.4112941) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2824077) q[2];
sx q[2];
rz(-1.7297598) q[2];
sx q[2];
rz(-0.32260103) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0301789) q[1];
sx q[1];
rz(-1.7708659) q[1];
sx q[1];
rz(-0.56758893) q[1];
rz(-1.8459122) q[3];
sx q[3];
rz(-1.0529622) q[3];
sx q[3];
rz(-2.8003611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0573132) q[2];
sx q[2];
rz(-1.7417615) q[2];
sx q[2];
rz(2.9220667) q[2];
rz(0.38671842) q[3];
sx q[3];
rz(-2.3733449) q[3];
sx q[3];
rz(0.99532551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0893843) q[0];
sx q[0];
rz(-1.5445671) q[0];
sx q[0];
rz(-2.9242933) q[0];
rz(3.1106588) q[1];
sx q[1];
rz(-2.5035281) q[1];
sx q[1];
rz(0.6589748) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12595183) q[0];
sx q[0];
rz(-0.77227393) q[0];
sx q[0];
rz(-0.036235972) q[0];
x q[1];
rz(0.72139229) q[2];
sx q[2];
rz(-2.112769) q[2];
sx q[2];
rz(2.0253571) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.98098552) q[1];
sx q[1];
rz(-1.7874582) q[1];
sx q[1];
rz(-0.76088455) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92774763) q[3];
sx q[3];
rz(-1.5156931) q[3];
sx q[3];
rz(-2.7411214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4387681) q[2];
sx q[2];
rz(-1.3805026) q[2];
sx q[2];
rz(-2.8213815) q[2];
rz(-1.6163588) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(-1.8922136) q[3];
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
rz(pi/2) q[3];
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
rz(2.3808688) q[0];
sx q[0];
rz(-0.18146935) q[0];
sx q[0];
rz(-2.4587801) q[0];
rz(-2.071351) q[1];
sx q[1];
rz(-1.2425334) q[1];
sx q[1];
rz(1.9314996) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72737981) q[0];
sx q[0];
rz(-0.41398898) q[0];
sx q[0];
rz(2.1075641) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96610539) q[2];
sx q[2];
rz(-1.7793057) q[2];
sx q[2];
rz(2.1272121) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6283419) q[1];
sx q[1];
rz(-1.2770137) q[1];
sx q[1];
rz(0.88900868) q[1];
x q[2];
rz(-2.5217767) q[3];
sx q[3];
rz(-1.8780939) q[3];
sx q[3];
rz(-0.67051552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5658297) q[2];
sx q[2];
rz(-2.1255707) q[2];
sx q[2];
rz(2.5749717) q[2];
rz(-2.2120655) q[3];
sx q[3];
rz(-2.1598787) q[3];
sx q[3];
rz(1.367759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.728445) q[0];
sx q[0];
rz(-2.8840273) q[0];
sx q[0];
rz(1.7096827) q[0];
rz(-0.52629772) q[1];
sx q[1];
rz(-0.52733517) q[1];
sx q[1];
rz(-2.3419103) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4096217) q[0];
sx q[0];
rz(-2.7831166) q[0];
sx q[0];
rz(1.326807) q[0];
rz(-pi) q[1];
rz(-1.4893555) q[2];
sx q[2];
rz(-1.4105721) q[2];
sx q[2];
rz(-0.83934957) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4754776) q[1];
sx q[1];
rz(-2.0787422) q[1];
sx q[1];
rz(2.75027) q[1];
x q[2];
rz(-1.6465854) q[3];
sx q[3];
rz(-1.4255376) q[3];
sx q[3];
rz(0.69243542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2852823) q[2];
sx q[2];
rz(-2.6976863) q[2];
sx q[2];
rz(2.4463859) q[2];
rz(0.55784145) q[3];
sx q[3];
rz(-1.7772243) q[3];
sx q[3];
rz(-2.4485596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(3.1391115) q[2];
sx q[2];
rz(-2.8056792) q[2];
sx q[2];
rz(0.16028595) q[2];
rz(-2.9813319) q[3];
sx q[3];
rz(-0.80453034) q[3];
sx q[3];
rz(2.9321032) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
