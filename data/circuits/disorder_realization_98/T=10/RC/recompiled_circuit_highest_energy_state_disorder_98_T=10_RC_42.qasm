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
rz(0.75818169) q[0];
sx q[0];
rz(3.5092847) q[0];
sx q[0];
rz(7.379471) q[0];
rz(0.81049377) q[1];
sx q[1];
rz(-0.23263045) q[1];
sx q[1];
rz(-2.0356324) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61845022) q[0];
sx q[0];
rz(-1.179427) q[0];
sx q[0];
rz(-2.6871293) q[0];
rz(-pi) q[1];
rz(-2.8739086) q[2];
sx q[2];
rz(-1.6639198) q[2];
sx q[2];
rz(-2.6665319) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.65484337) q[1];
sx q[1];
rz(-0.037529201) q[1];
sx q[1];
rz(1.4629745) q[1];
rz(-1.5777784) q[3];
sx q[3];
rz(-1.8551702) q[3];
sx q[3];
rz(0.46753903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7652863) q[2];
sx q[2];
rz(-0.47069612) q[2];
sx q[2];
rz(-2.7619696) q[2];
rz(1.5073353) q[3];
sx q[3];
rz(-1.0999271) q[3];
sx q[3];
rz(0.96989337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64604243) q[0];
sx q[0];
rz(-0.33691418) q[0];
sx q[0];
rz(-0.90091339) q[0];
rz(3.0100929) q[1];
sx q[1];
rz(-1.9410746) q[1];
sx q[1];
rz(0.44201717) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7666047) q[0];
sx q[0];
rz(-2.8061013) q[0];
sx q[0];
rz(-2.8906402) q[0];
rz(1.5698264) q[2];
sx q[2];
rz(-1.5661777) q[2];
sx q[2];
rz(-0.49170353) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17168853) q[1];
sx q[1];
rz(-2.0398047) q[1];
sx q[1];
rz(-0.63290872) q[1];
x q[2];
rz(1.6824746) q[3];
sx q[3];
rz(-1.1726716) q[3];
sx q[3];
rz(-2.3071837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.76571959) q[2];
sx q[2];
rz(-1.7123545) q[2];
sx q[2];
rz(-2.5326552) q[2];
rz(0.40306148) q[3];
sx q[3];
rz(-1.1761222) q[3];
sx q[3];
rz(-0.5016996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8159863) q[0];
sx q[0];
rz(-0.41929647) q[0];
sx q[0];
rz(-2.0035279) q[0];
rz(-1.5886547) q[1];
sx q[1];
rz(-0.76858968) q[1];
sx q[1];
rz(2.3345711) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0332542) q[0];
sx q[0];
rz(-1.6699808) q[0];
sx q[0];
rz(1.5434815) q[0];
x q[1];
rz(-2.8093178) q[2];
sx q[2];
rz(-2.6746779) q[2];
sx q[2];
rz(0.98540598) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.623329) q[1];
sx q[1];
rz(-2.9179472) q[1];
sx q[1];
rz(-2.5751051) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.85523) q[3];
sx q[3];
rz(-1.2060646) q[3];
sx q[3];
rz(1.3160365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1300065) q[2];
sx q[2];
rz(-0.61379543) q[2];
sx q[2];
rz(1.9278795) q[2];
rz(3.0774434) q[3];
sx q[3];
rz(-1.6545273) q[3];
sx q[3];
rz(-0.00042644342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59378687) q[0];
sx q[0];
rz(-1.8812027) q[0];
sx q[0];
rz(2.2547145) q[0];
rz(-0.15874323) q[1];
sx q[1];
rz(-2.6491149) q[1];
sx q[1];
rz(1.8633206) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0510382) q[0];
sx q[0];
rz(-2.0060385) q[0];
sx q[0];
rz(-0.79751517) q[0];
x q[1];
rz(2.0416987) q[2];
sx q[2];
rz(-2.0713901) q[2];
sx q[2];
rz(-0.72906993) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.63250971) q[1];
sx q[1];
rz(-2.4466568) q[1];
sx q[1];
rz(-1.664723) q[1];
rz(-pi) q[2];
rz(-2.9904891) q[3];
sx q[3];
rz(-2.1204815) q[3];
sx q[3];
rz(-1.0719887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0025803) q[2];
sx q[2];
rz(-0.85550344) q[2];
sx q[2];
rz(-0.95376897) q[2];
rz(1.194713) q[3];
sx q[3];
rz(-2.298893) q[3];
sx q[3];
rz(2.6660624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3139528) q[0];
sx q[0];
rz(-2.4251745) q[0];
sx q[0];
rz(0.59999505) q[0];
rz(-2.4586239) q[1];
sx q[1];
rz(-2.3536847) q[1];
sx q[1];
rz(-0.33040985) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1106772) q[0];
sx q[0];
rz(-0.85665032) q[0];
sx q[0];
rz(2.899709) q[0];
rz(-pi) q[1];
rz(2.1938964) q[2];
sx q[2];
rz(-2.8665339) q[2];
sx q[2];
rz(0.81240053) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.40881316) q[1];
sx q[1];
rz(-2.6169852) q[1];
sx q[1];
rz(1.9636159) q[1];
rz(-pi) q[2];
rz(-2.7269269) q[3];
sx q[3];
rz(-2.2288481) q[3];
sx q[3];
rz(-1.0533223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7052475) q[2];
sx q[2];
rz(-2.0630431) q[2];
sx q[2];
rz(-0.62502965) q[2];
rz(2.446512) q[3];
sx q[3];
rz(-1.2806226) q[3];
sx q[3];
rz(0.052791031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92796749) q[0];
sx q[0];
rz(-1.1518421) q[0];
sx q[0];
rz(-1.7684162) q[0];
rz(-1.7534509) q[1];
sx q[1];
rz(-0.87166798) q[1];
sx q[1];
rz(-1.6067827) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2508188) q[0];
sx q[0];
rz(-1.5355331) q[0];
sx q[0];
rz(1.4232487) q[0];
rz(-pi) q[1];
rz(2.0892145) q[2];
sx q[2];
rz(-0.3568584) q[2];
sx q[2];
rz(-1.269358) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.36027137) q[1];
sx q[1];
rz(-2.2151673) q[1];
sx q[1];
rz(2.086198) q[1];
rz(-1.8880879) q[3];
sx q[3];
rz(-0.41966885) q[3];
sx q[3];
rz(-2.9515226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9390949) q[2];
sx q[2];
rz(-1.9715344) q[2];
sx q[2];
rz(-1.6004174) q[2];
rz(2.1206858) q[3];
sx q[3];
rz(-1.7424135) q[3];
sx q[3];
rz(-0.30103621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0114667) q[0];
sx q[0];
rz(-3.0045894) q[0];
sx q[0];
rz(-2.2504508) q[0];
rz(-2.4961684) q[1];
sx q[1];
rz(-1.3963457) q[1];
sx q[1];
rz(-0.83271629) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37207212) q[0];
sx q[0];
rz(-1.4345048) q[0];
sx q[0];
rz(2.6762677) q[0];
rz(1.4786167) q[2];
sx q[2];
rz(-2.5155009) q[2];
sx q[2];
rz(2.7220059) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.558185) q[1];
sx q[1];
rz(-1.3761569) q[1];
sx q[1];
rz(-2.6411112) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54033684) q[3];
sx q[3];
rz(-1.9129921) q[3];
sx q[3];
rz(-1.0409174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5523395) q[2];
sx q[2];
rz(-0.66698843) q[2];
sx q[2];
rz(-1.5480631) q[2];
rz(-2.7175236) q[3];
sx q[3];
rz(-1.721902) q[3];
sx q[3];
rz(-2.8467395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9645204) q[0];
sx q[0];
rz(-1.1966713) q[0];
sx q[0];
rz(-0.82343423) q[0];
rz(1.1973165) q[1];
sx q[1];
rz(-1.4875393) q[1];
sx q[1];
rz(-1.6808602) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0076559) q[0];
sx q[0];
rz(-1.888935) q[0];
sx q[0];
rz(0.11548059) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97779556) q[2];
sx q[2];
rz(-1.7109181) q[2];
sx q[2];
rz(-1.8210379) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4464691) q[1];
sx q[1];
rz(-0.90988509) q[1];
sx q[1];
rz(1.5170694) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4183059) q[3];
sx q[3];
rz(-2.6761665) q[3];
sx q[3];
rz(1.775072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0921649) q[2];
sx q[2];
rz(-2.4804513) q[2];
sx q[2];
rz(0.1235505) q[2];
rz(2.3030247) q[3];
sx q[3];
rz(-2.0288012) q[3];
sx q[3];
rz(0.53028321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-3.0274444) q[0];
sx q[0];
rz(-1.0076948) q[0];
sx q[0];
rz(-1.4134407) q[0];
rz(2.24276) q[1];
sx q[1];
rz(-0.90679589) q[1];
sx q[1];
rz(2.5487505) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7031738) q[0];
sx q[0];
rz(-0.98383622) q[0];
sx q[0];
rz(-2.274375) q[0];
rz(2.7024621) q[2];
sx q[2];
rz(-2.9408717) q[2];
sx q[2];
rz(1.8619271) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.71708116) q[1];
sx q[1];
rz(-2.6047978) q[1];
sx q[1];
rz(-2.2473141) q[1];
rz(-0.96372202) q[3];
sx q[3];
rz(-0.98120171) q[3];
sx q[3];
rz(-2.251925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8574519) q[2];
sx q[2];
rz(-2.393674) q[2];
sx q[2];
rz(-0.22135529) q[2];
rz(2.0981233) q[3];
sx q[3];
rz(-2.2011493) q[3];
sx q[3];
rz(2.7531457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2986044) q[0];
sx q[0];
rz(-2.8589111) q[0];
sx q[0];
rz(-2.2931732) q[0];
rz(-0.09659718) q[1];
sx q[1];
rz(-1.074147) q[1];
sx q[1];
rz(0.75327795) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71614186) q[0];
sx q[0];
rz(-1.4708859) q[0];
sx q[0];
rz(-1.3425171) q[0];
rz(-pi) q[1];
rz(3.1084849) q[2];
sx q[2];
rz(-0.82947846) q[2];
sx q[2];
rz(-0.25698369) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2172474) q[1];
sx q[1];
rz(-2.0567523) q[1];
sx q[1];
rz(2.7107581) q[1];
x q[2];
rz(2.0263591) q[3];
sx q[3];
rz(-0.90901819) q[3];
sx q[3];
rz(1.7152804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0775371) q[2];
sx q[2];
rz(-0.83751837) q[2];
sx q[2];
rz(0.45292863) q[2];
rz(-2.5405267) q[3];
sx q[3];
rz(-1.6744924) q[3];
sx q[3];
rz(-2.1452904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8375028) q[0];
sx q[0];
rz(-1.6927728) q[0];
sx q[0];
rz(0.46020831) q[0];
rz(-2.9651463) q[1];
sx q[1];
rz(-0.23527589) q[1];
sx q[1];
rz(-1.489524) q[1];
rz(1.5701736) q[2];
sx q[2];
rz(-1.9194308) q[2];
sx q[2];
rz(0.44375833) q[2];
rz(-2.413977) q[3];
sx q[3];
rz(-2.03491) q[3];
sx q[3];
rz(2.2274957) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
