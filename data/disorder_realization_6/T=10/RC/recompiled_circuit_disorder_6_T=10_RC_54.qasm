OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52580994) q[0];
sx q[0];
rz(-1.7237741) q[0];
sx q[0];
rz(-0.56086993) q[0];
rz(-2.0286735) q[1];
sx q[1];
rz(-1.3781883) q[1];
sx q[1];
rz(-1.2150432) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5117447) q[0];
sx q[0];
rz(-1.6748322) q[0];
sx q[0];
rz(1.758979) q[0];
x q[1];
rz(0.36429976) q[2];
sx q[2];
rz(-2.4476353) q[2];
sx q[2];
rz(2.7147164) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4814264) q[1];
sx q[1];
rz(-0.62454849) q[1];
sx q[1];
rz(2.156483) q[1];
rz(-pi) q[2];
rz(1.3052985) q[3];
sx q[3];
rz(-1.7146829) q[3];
sx q[3];
rz(0.063751566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3540196) q[2];
sx q[2];
rz(-0.95280567) q[2];
sx q[2];
rz(2.9585178) q[2];
rz(-2.7637774) q[3];
sx q[3];
rz(-2.0928045) q[3];
sx q[3];
rz(-2.8474076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29782444) q[0];
sx q[0];
rz(-2.4968708) q[0];
sx q[0];
rz(-3.0644754) q[0];
rz(0.33879694) q[1];
sx q[1];
rz(-1.1145376) q[1];
sx q[1];
rz(-1.6024626) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3829271) q[0];
sx q[0];
rz(-2.5490952) q[0];
sx q[0];
rz(-1.4325607) q[0];
rz(-2.9039731) q[2];
sx q[2];
rz(-0.75287205) q[2];
sx q[2];
rz(-1.876229) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6290366) q[1];
sx q[1];
rz(-1.3965544) q[1];
sx q[1];
rz(1.0838572) q[1];
x q[2];
rz(1.8061403) q[3];
sx q[3];
rz(-1.1574405) q[3];
sx q[3];
rz(-1.8542765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.845528) q[2];
sx q[2];
rz(-1.277593) q[2];
sx q[2];
rz(0.65845931) q[2];
rz(2.9902839) q[3];
sx q[3];
rz(-2.1189809) q[3];
sx q[3];
rz(-0.69491274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028458683) q[0];
sx q[0];
rz(-0.78335339) q[0];
sx q[0];
rz(2.7084896) q[0];
rz(-1.9494879) q[1];
sx q[1];
rz(-1.2116218) q[1];
sx q[1];
rz(2.5862397) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1486737) q[0];
sx q[0];
rz(-2.2745471) q[0];
sx q[0];
rz(0.91627319) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0981512) q[2];
sx q[2];
rz(-0.61908365) q[2];
sx q[2];
rz(1.4738136) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8086116) q[1];
sx q[1];
rz(-2.4858027) q[1];
sx q[1];
rz(1.3036742) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2461353) q[3];
sx q[3];
rz(-0.56926308) q[3];
sx q[3];
rz(0.07490052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4042523) q[2];
sx q[2];
rz(-0.78812391) q[2];
sx q[2];
rz(-1.2505442) q[2];
rz(-0.2441497) q[3];
sx q[3];
rz(-1.282225) q[3];
sx q[3];
rz(-1.4499433) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26043949) q[0];
sx q[0];
rz(-0.45757159) q[0];
sx q[0];
rz(-2.326791) q[0];
rz(1.3793777) q[1];
sx q[1];
rz(-2.791399) q[1];
sx q[1];
rz(0.25517685) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6738621) q[0];
sx q[0];
rz(-2.6669589) q[0];
sx q[0];
rz(2.4509096) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63919477) q[2];
sx q[2];
rz(-1.8393469) q[2];
sx q[2];
rz(3.064379) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7123588) q[1];
sx q[1];
rz(-0.36839596) q[1];
sx q[1];
rz(-0.952094) q[1];
rz(-0.051024036) q[3];
sx q[3];
rz(-2.0912366) q[3];
sx q[3];
rz(-0.40363064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2531551) q[2];
sx q[2];
rz(-1.5909114) q[2];
sx q[2];
rz(-2.9684084) q[2];
rz(-0.52982461) q[3];
sx q[3];
rz(-2.9960222) q[3];
sx q[3];
rz(-0.1023275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2816876) q[0];
sx q[0];
rz(-1.4929993) q[0];
sx q[0];
rz(-1.3758855) q[0];
rz(-1.8638523) q[1];
sx q[1];
rz(-0.81218305) q[1];
sx q[1];
rz(3.0854991) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2984021) q[0];
sx q[0];
rz(-2.139233) q[0];
sx q[0];
rz(1.1355023) q[0];
rz(-pi) q[1];
rz(0.9390097) q[2];
sx q[2];
rz(-2.5917705) q[2];
sx q[2];
rz(2.3914571) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4896506) q[1];
sx q[1];
rz(-0.30311668) q[1];
sx q[1];
rz(-3.0928844) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.089245307) q[3];
sx q[3];
rz(-1.0122932) q[3];
sx q[3];
rz(-2.9104779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.00099480199) q[2];
sx q[2];
rz(-0.91789118) q[2];
sx q[2];
rz(0.43219217) q[2];
rz(-2.2473992) q[3];
sx q[3];
rz(-1.0995355) q[3];
sx q[3];
rz(1.4661219) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11319259) q[0];
sx q[0];
rz(-0.88554651) q[0];
sx q[0];
rz(-2.4940441) q[0];
rz(1.2619069) q[1];
sx q[1];
rz(-1.4636661) q[1];
sx q[1];
rz(-2.1870959) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.944825) q[0];
sx q[0];
rz(-1.4155354) q[0];
sx q[0];
rz(-1.0374271) q[0];
x q[1];
rz(-0.034770413) q[2];
sx q[2];
rz(-2.2010942) q[2];
sx q[2];
rz(0.45331732) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8050025) q[1];
sx q[1];
rz(-2.046236) q[1];
sx q[1];
rz(-0.57979433) q[1];
rz(-pi) q[2];
rz(-1.2797221) q[3];
sx q[3];
rz(-1.4143412) q[3];
sx q[3];
rz(-2.1978956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.59297562) q[2];
sx q[2];
rz(-1.9081215) q[2];
sx q[2];
rz(-2.0992289) q[2];
rz(-2.7029165) q[3];
sx q[3];
rz(-2.091566) q[3];
sx q[3];
rz(-1.3180805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.8577268) q[0];
sx q[0];
rz(-2.9086869) q[0];
sx q[0];
rz(-0.74321157) q[0];
rz(1.5076393) q[1];
sx q[1];
rz(-0.71989027) q[1];
sx q[1];
rz(0.61002237) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3726495) q[0];
sx q[0];
rz(-2.5805051) q[0];
sx q[0];
rz(-2.6555496) q[0];
rz(-1.4201944) q[2];
sx q[2];
rz(-0.65956958) q[2];
sx q[2];
rz(1.9285551) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2513189) q[1];
sx q[1];
rz(-0.90369019) q[1];
sx q[1];
rz(-3.1097079) q[1];
x q[2];
rz(2.0537297) q[3];
sx q[3];
rz(-0.42290877) q[3];
sx q[3];
rz(-1.9382697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8053749) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(0.91840333) q[2];
rz(-1.5911128) q[3];
sx q[3];
rz(-0.95033002) q[3];
sx q[3];
rz(2.7526855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36088762) q[0];
sx q[0];
rz(-2.4724859) q[0];
sx q[0];
rz(1.5135182) q[0];
rz(2.6121415) q[1];
sx q[1];
rz(-2.0748731) q[1];
sx q[1];
rz(2.4050074) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37305957) q[0];
sx q[0];
rz(-3.1096418) q[0];
sx q[0];
rz(1.91747) q[0];
rz(-1.2484776) q[2];
sx q[2];
rz(-1.6741447) q[2];
sx q[2];
rz(-1.5600187) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8383933) q[1];
sx q[1];
rz(-1.1078849) q[1];
sx q[1];
rz(-1.6846912) q[1];
x q[2];
rz(0.96111091) q[3];
sx q[3];
rz(-0.52934066) q[3];
sx q[3];
rz(1.5101658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1071876) q[2];
sx q[2];
rz(-1.2074869) q[2];
sx q[2];
rz(0.68391189) q[2];
rz(-1.2290139) q[3];
sx q[3];
rz(-1.3701655) q[3];
sx q[3];
rz(-1.7470523) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9443611) q[0];
sx q[0];
rz(-1.5988388) q[0];
sx q[0];
rz(2.9558682) q[0];
rz(-2.1445403) q[1];
sx q[1];
rz(-1.8763708) q[1];
sx q[1];
rz(-2.396778) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54430994) q[0];
sx q[0];
rz(-2.3047857) q[0];
sx q[0];
rz(2.701176) q[0];
rz(0.79297519) q[2];
sx q[2];
rz(-1.8599531) q[2];
sx q[2];
rz(-2.5536429) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3639431) q[1];
sx q[1];
rz(-1.793982) q[1];
sx q[1];
rz(1.2121483) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7471077) q[3];
sx q[3];
rz(-0.82291616) q[3];
sx q[3];
rz(-2.8701973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1054489) q[2];
sx q[2];
rz(-2.3858586) q[2];
sx q[2];
rz(1.194681) q[2];
rz(0.99669325) q[3];
sx q[3];
rz(-1.2160622) q[3];
sx q[3];
rz(-2.1452346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3982518) q[0];
sx q[0];
rz(-0.78813362) q[0];
sx q[0];
rz(-2.7375896) q[0];
rz(0.031127302) q[1];
sx q[1];
rz(-1.6571836) q[1];
sx q[1];
rz(-1.1709447) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084328018) q[0];
sx q[0];
rz(-1.143976) q[0];
sx q[0];
rz(1.5469993) q[0];
rz(-2.3659336) q[2];
sx q[2];
rz(-1.7663029) q[2];
sx q[2];
rz(2.069371) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3552637) q[1];
sx q[1];
rz(-1.5621645) q[1];
sx q[1];
rz(-1.5934056) q[1];
rz(-pi) q[2];
rz(2.714614) q[3];
sx q[3];
rz(-1.8378165) q[3];
sx q[3];
rz(3.0468575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6955473) q[2];
sx q[2];
rz(-1.7798767) q[2];
sx q[2];
rz(0.5919624) q[2];
rz(2.5752318) q[3];
sx q[3];
rz(-0.16470328) q[3];
sx q[3];
rz(-1.5238354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82407172) q[0];
sx q[0];
rz(-2.1614647) q[0];
sx q[0];
rz(1.9807057) q[0];
rz(-3.042165) q[1];
sx q[1];
rz(-1.8933404) q[1];
sx q[1];
rz(1.0642687) q[1];
rz(-0.912491) q[2];
sx q[2];
rz(-2.1952663) q[2];
sx q[2];
rz(-1.2637539) q[2];
rz(0.016146544) q[3];
sx q[3];
rz(-1.9042249) q[3];
sx q[3];
rz(-0.92845542) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
