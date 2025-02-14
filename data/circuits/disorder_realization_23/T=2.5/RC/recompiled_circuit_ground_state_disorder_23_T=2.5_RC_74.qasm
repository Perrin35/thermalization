OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5601114) q[0];
sx q[0];
rz(-1.2553296) q[0];
sx q[0];
rz(-2.6502967) q[0];
rz(-3.4499912) q[1];
sx q[1];
rz(1.8412794) q[1];
sx q[1];
rz(12.507764) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1442382) q[0];
sx q[0];
rz(-0.73737007) q[0];
sx q[0];
rz(-2.2585346) q[0];
x q[1];
rz(0.7841447) q[2];
sx q[2];
rz(-1.2375187) q[2];
sx q[2];
rz(3.0212457) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9323001) q[1];
sx q[1];
rz(-1.2697311) q[1];
sx q[1];
rz(-1.4026227) q[1];
rz(-1.2514312) q[3];
sx q[3];
rz(-1.6465558) q[3];
sx q[3];
rz(-3.0204242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0356174) q[2];
sx q[2];
rz(-2.3507037) q[2];
sx q[2];
rz(-0.38113511) q[2];
rz(-0.032958802) q[3];
sx q[3];
rz(-1.4573174) q[3];
sx q[3];
rz(1.0491252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9829262) q[0];
sx q[0];
rz(-2.6148112) q[0];
sx q[0];
rz(-0.43989936) q[0];
rz(0.061821763) q[1];
sx q[1];
rz(-1.6114176) q[1];
sx q[1];
rz(0.27438146) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2915115) q[0];
sx q[0];
rz(-0.78118229) q[0];
sx q[0];
rz(2.2456332) q[0];
x q[1];
rz(2.0437061) q[2];
sx q[2];
rz(-0.47240546) q[2];
sx q[2];
rz(2.1771569) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.134717) q[1];
sx q[1];
rz(-1.1235079) q[1];
sx q[1];
rz(3.0577117) q[1];
x q[2];
rz(0.54631375) q[3];
sx q[3];
rz(-1.3689181) q[3];
sx q[3];
rz(-0.047919433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2939833) q[2];
sx q[2];
rz(-3.1276438) q[2];
sx q[2];
rz(0.070276109) q[2];
rz(2.516732) q[3];
sx q[3];
rz(-1.8127541) q[3];
sx q[3];
rz(-3.0666053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4570419) q[0];
sx q[0];
rz(-1.6383189) q[0];
sx q[0];
rz(0.3918089) q[0];
rz(-2.1458972) q[1];
sx q[1];
rz(-0.83637801) q[1];
sx q[1];
rz(3.0973184) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56681117) q[0];
sx q[0];
rz(-1.0100845) q[0];
sx q[0];
rz(-0.36578481) q[0];
rz(-pi) q[1];
rz(2.2703553) q[2];
sx q[2];
rz(-0.92513621) q[2];
sx q[2];
rz(-1.439656) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.6926319) q[1];
sx q[1];
rz(-1.7794183) q[1];
sx q[1];
rz(1.7600928) q[1];
x q[2];
rz(-0.1416154) q[3];
sx q[3];
rz(-0.60185233) q[3];
sx q[3];
rz(0.54094368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5135045) q[2];
sx q[2];
rz(-0.91614437) q[2];
sx q[2];
rz(1.2543031) q[2];
rz(0.11145505) q[3];
sx q[3];
rz(-0.97193757) q[3];
sx q[3];
rz(2.286262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3067538) q[0];
sx q[0];
rz(-0.67973891) q[0];
sx q[0];
rz(-2.2656031) q[0];
rz(2.7442979) q[1];
sx q[1];
rz(-1.5715503) q[1];
sx q[1];
rz(-1.809459) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5529163) q[0];
sx q[0];
rz(-1.1015) q[0];
sx q[0];
rz(-1.3568715) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1213673) q[2];
sx q[2];
rz(-0.75006333) q[2];
sx q[2];
rz(-1.5062065) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6369824) q[1];
sx q[1];
rz(-2.9688997) q[1];
sx q[1];
rz(2.6673698) q[1];
rz(1.9537266) q[3];
sx q[3];
rz(-0.57102244) q[3];
sx q[3];
rz(-1.8117222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1490271) q[2];
sx q[2];
rz(-1.6529447) q[2];
sx q[2];
rz(-1.3943025) q[2];
rz(-2.1483138) q[3];
sx q[3];
rz(-1.1634049) q[3];
sx q[3];
rz(-1.8786028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.29436) q[0];
sx q[0];
rz(-0.17437409) q[0];
sx q[0];
rz(-2.2827523) q[0];
rz(2.7186588) q[1];
sx q[1];
rz(-0.49086389) q[1];
sx q[1];
rz(-1.4609969) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5678755) q[0];
sx q[0];
rz(-1.4481312) q[0];
sx q[0];
rz(-0.50267641) q[0];
rz(-1.5761574) q[2];
sx q[2];
rz(-2.6963384) q[2];
sx q[2];
rz(1.3878979) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.40845937) q[1];
sx q[1];
rz(-0.28318048) q[1];
sx q[1];
rz(2.3896809) q[1];
rz(-2.0764396) q[3];
sx q[3];
rz(-1.6958456) q[3];
sx q[3];
rz(1.7495383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.6964263) q[2];
sx q[2];
rz(-0.40881613) q[2];
sx q[2];
rz(-1.6430631) q[2];
rz(-2.0131352) q[3];
sx q[3];
rz(-2.0339637) q[3];
sx q[3];
rz(-2.4988153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4437101) q[0];
sx q[0];
rz(-1.3297798) q[0];
sx q[0];
rz(0.71769303) q[0];
rz(0.36092654) q[1];
sx q[1];
rz(-2.2310427) q[1];
sx q[1];
rz(-2.8275729) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0180394) q[0];
sx q[0];
rz(-1.529487) q[0];
sx q[0];
rz(2.6459751) q[0];
rz(1.6144361) q[2];
sx q[2];
rz(-2.4167762) q[2];
sx q[2];
rz(-2.5284655) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.18923387) q[1];
sx q[1];
rz(-1.2554607) q[1];
sx q[1];
rz(1.0811514) q[1];
x q[2];
rz(-1.3889503) q[3];
sx q[3];
rz(-0.80917984) q[3];
sx q[3];
rz(-1.2004747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6641984) q[2];
sx q[2];
rz(-1.3496642) q[2];
sx q[2];
rz(2.9273709) q[2];
rz(0.91841206) q[3];
sx q[3];
rz(-2.1954506) q[3];
sx q[3];
rz(1.4001747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65900954) q[0];
sx q[0];
rz(-2.5795586) q[0];
sx q[0];
rz(2.5157978) q[0];
rz(1.91045) q[1];
sx q[1];
rz(-1.2673255) q[1];
sx q[1];
rz(2.321718) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17894402) q[0];
sx q[0];
rz(-1.8838619) q[0];
sx q[0];
rz(0.59525916) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.49687044) q[2];
sx q[2];
rz(-2.3665603) q[2];
sx q[2];
rz(0.97737038) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.956683) q[1];
sx q[1];
rz(-1.239983) q[1];
sx q[1];
rz(-0.079903101) q[1];
rz(-pi) q[2];
rz(0.094535703) q[3];
sx q[3];
rz(-1.5889639) q[3];
sx q[3];
rz(-1.2430735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.95255533) q[2];
sx q[2];
rz(-1.8841212) q[2];
sx q[2];
rz(2.7541568) q[2];
rz(0.45346013) q[3];
sx q[3];
rz(-2.267434) q[3];
sx q[3];
rz(0.15997729) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3979935) q[0];
sx q[0];
rz(-2.0026119) q[0];
sx q[0];
rz(-2.0177662) q[0];
rz(-0.75346142) q[1];
sx q[1];
rz(-1.9517784) q[1];
sx q[1];
rz(-1.168728) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74334621) q[0];
sx q[0];
rz(-0.69730824) q[0];
sx q[0];
rz(-1.4278317) q[0];
rz(-1.1134603) q[2];
sx q[2];
rz(-1.5761318) q[2];
sx q[2];
rz(-2.920546) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1064042) q[1];
sx q[1];
rz(-0.6995844) q[1];
sx q[1];
rz(-2.0856218) q[1];
rz(0.91442502) q[3];
sx q[3];
rz(-1.1946988) q[3];
sx q[3];
rz(0.25901435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0231102) q[2];
sx q[2];
rz(-2.1239943) q[2];
sx q[2];
rz(-0.47038308) q[2];
rz(-1.6810301) q[3];
sx q[3];
rz(-1.2602256) q[3];
sx q[3];
rz(0.98106658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39101446) q[0];
sx q[0];
rz(-1.7711201) q[0];
sx q[0];
rz(-1.4960666) q[0];
rz(-2.0059313) q[1];
sx q[1];
rz(-1.2874425) q[1];
sx q[1];
rz(-0.48887238) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2100567) q[0];
sx q[0];
rz(-1.6029412) q[0];
sx q[0];
rz(-1.4308306) q[0];
rz(-pi) q[1];
rz(-3.0453072) q[2];
sx q[2];
rz(-2.0803148) q[2];
sx q[2];
rz(-2.4063039) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.012163398) q[1];
sx q[1];
rz(-1.5930854) q[1];
sx q[1];
rz(0.049456923) q[1];
x q[2];
rz(-3.0472894) q[3];
sx q[3];
rz(-0.62178388) q[3];
sx q[3];
rz(1.0434542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.44498542) q[2];
sx q[2];
rz(-0.61034909) q[2];
sx q[2];
rz(-2.0390017) q[2];
rz(1.6164814) q[3];
sx q[3];
rz(-2.3307255) q[3];
sx q[3];
rz(-1.471224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.1333604) q[0];
sx q[0];
rz(-2.6120549) q[0];
sx q[0];
rz(-1.3326921) q[0];
rz(-0.38707271) q[1];
sx q[1];
rz(-1.8862855) q[1];
sx q[1];
rz(-2.0852087) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1176193) q[0];
sx q[0];
rz(-1.4314382) q[0];
sx q[0];
rz(-1.3246714) q[0];
rz(2.8661997) q[2];
sx q[2];
rz(-0.43245855) q[2];
sx q[2];
rz(1.5713009) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3247128) q[1];
sx q[1];
rz(-1.9512647) q[1];
sx q[1];
rz(2.3808384) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4293618) q[3];
sx q[3];
rz(-0.81347403) q[3];
sx q[3];
rz(1.939919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1308088) q[2];
sx q[2];
rz(-0.84785145) q[2];
sx q[2];
rz(-1.7424142) q[2];
rz(-2.0684659) q[3];
sx q[3];
rz(-1.224204) q[3];
sx q[3];
rz(-2.6789902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90049967) q[0];
sx q[0];
rz(-1.4971965) q[0];
sx q[0];
rz(1.5370488) q[0];
rz(2.4201139) q[1];
sx q[1];
rz(-2.571048) q[1];
sx q[1];
rz(-1.2068988) q[1];
rz(-1.9611364) q[2];
sx q[2];
rz(-1.9705638) q[2];
sx q[2];
rz(0.86314291) q[2];
rz(-0.10151699) q[3];
sx q[3];
rz(-2.3675167) q[3];
sx q[3];
rz(1.7270796) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
