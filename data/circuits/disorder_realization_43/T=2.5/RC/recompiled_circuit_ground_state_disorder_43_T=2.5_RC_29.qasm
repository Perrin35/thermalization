OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3620152) q[0];
sx q[0];
rz(0.22688046) q[0];
sx q[0];
rz(14.332834) q[0];
rz(3.0472164) q[1];
sx q[1];
rz(-0.91369319) q[1];
sx q[1];
rz(0.28092608) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37051113) q[0];
sx q[0];
rz(-0.91508741) q[0];
sx q[0];
rz(-2.9408668) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0269555) q[2];
sx q[2];
rz(-1.2817941) q[2];
sx q[2];
rz(-2.8665045) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0795796) q[1];
sx q[1];
rz(-1.307319) q[1];
sx q[1];
rz(2.0017712) q[1];
rz(2.6219764) q[3];
sx q[3];
rz(-2.5027788) q[3];
sx q[3];
rz(2.6878302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1137587) q[2];
sx q[2];
rz(-1.4049302) q[2];
sx q[2];
rz(-3.0244381) q[2];
rz(2.8095918) q[3];
sx q[3];
rz(-0.74688512) q[3];
sx q[3];
rz(-2.3565256) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4431045) q[0];
sx q[0];
rz(-2.3790058) q[0];
sx q[0];
rz(2.7681328) q[0];
rz(0.39730486) q[1];
sx q[1];
rz(-1.1445069) q[1];
sx q[1];
rz(2.4041046) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0076548) q[0];
sx q[0];
rz(-2.2077401) q[0];
sx q[0];
rz(0.68151125) q[0];
rz(-pi) q[1];
rz(0.55464427) q[2];
sx q[2];
rz(-2.0337542) q[2];
sx q[2];
rz(-2.2210768) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8792845) q[1];
sx q[1];
rz(-1.4170987) q[1];
sx q[1];
rz(-2.5889695) q[1];
rz(0.37009671) q[3];
sx q[3];
rz(-0.79269275) q[3];
sx q[3];
rz(-1.1350138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.071659) q[2];
sx q[2];
rz(-1.5295014) q[2];
sx q[2];
rz(-2.7640479) q[2];
rz(-2.8318882) q[3];
sx q[3];
rz(-1.0891424) q[3];
sx q[3];
rz(1.496605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.51652235) q[0];
sx q[0];
rz(-2.3023038) q[0];
sx q[0];
rz(1.2155493) q[0];
rz(1.0141605) q[1];
sx q[1];
rz(-1.4233669) q[1];
sx q[1];
rz(-0.97022143) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.559186) q[0];
sx q[0];
rz(-1.0716039) q[0];
sx q[0];
rz(3.0954719) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2283111) q[2];
sx q[2];
rz(-2.0184419) q[2];
sx q[2];
rz(-0.44128445) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.41170317) q[1];
sx q[1];
rz(-1.3940934) q[1];
sx q[1];
rz(1.2248301) q[1];
rz(-pi) q[2];
rz(-0.30556314) q[3];
sx q[3];
rz(-2.773958) q[3];
sx q[3];
rz(0.43722269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0479451) q[2];
sx q[2];
rz(-2.2990172) q[2];
sx q[2];
rz(-2.688431) q[2];
rz(1.9122745) q[3];
sx q[3];
rz(-1.2776351) q[3];
sx q[3];
rz(0.17190988) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9446843) q[0];
sx q[0];
rz(-0.6854282) q[0];
sx q[0];
rz(-2.7759283) q[0];
rz(0.91224313) q[1];
sx q[1];
rz(-1.8625926) q[1];
sx q[1];
rz(-0.40547392) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.643754) q[0];
sx q[0];
rz(-0.23786834) q[0];
sx q[0];
rz(2.0980623) q[0];
rz(-2.3857152) q[2];
sx q[2];
rz(-1.4149932) q[2];
sx q[2];
rz(2.1550117) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7995389) q[1];
sx q[1];
rz(-0.56302445) q[1];
sx q[1];
rz(0.31837007) q[1];
rz(0.41678352) q[3];
sx q[3];
rz(-0.68282167) q[3];
sx q[3];
rz(-0.58931755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.071986467) q[2];
sx q[2];
rz(-1.7013841) q[2];
sx q[2];
rz(1.5427422) q[2];
rz(-0.37374464) q[3];
sx q[3];
rz(-1.6662686) q[3];
sx q[3];
rz(-1.4603978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5986346) q[0];
sx q[0];
rz(-2.3668508) q[0];
sx q[0];
rz(-2.4454818) q[0];
rz(-1.2593345) q[1];
sx q[1];
rz(-2.0399317) q[1];
sx q[1];
rz(-0.36516821) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8445963) q[0];
sx q[0];
rz(-2.0005324) q[0];
sx q[0];
rz(0.55083042) q[0];
rz(0.31814534) q[2];
sx q[2];
rz(-1.9650302) q[2];
sx q[2];
rz(2.1925558) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9112253) q[1];
sx q[1];
rz(-2.5960698) q[1];
sx q[1];
rz(-0.0026774252) q[1];
rz(2.5138084) q[3];
sx q[3];
rz(-1.9115698) q[3];
sx q[3];
rz(-0.22751156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.460707) q[2];
sx q[2];
rz(-1.1628217) q[2];
sx q[2];
rz(2.9648901) q[2];
rz(-1.7548615) q[3];
sx q[3];
rz(-1.0597798) q[3];
sx q[3];
rz(2.7669014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52727592) q[0];
sx q[0];
rz(-0.85955954) q[0];
sx q[0];
rz(1.5981307) q[0];
rz(1.3044926) q[1];
sx q[1];
rz(-1.0544798) q[1];
sx q[1];
rz(2.3354882) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5209608) q[0];
sx q[0];
rz(-1.3343108) q[0];
sx q[0];
rz(-2.3804401) q[0];
rz(-pi) q[1];
rz(0.97053846) q[2];
sx q[2];
rz(-1.1379062) q[2];
sx q[2];
rz(-0.25161703) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1404952) q[1];
sx q[1];
rz(-1.5211269) q[1];
sx q[1];
rz(-1.7752663) q[1];
rz(-0.73202912) q[3];
sx q[3];
rz(-2.424758) q[3];
sx q[3];
rz(-2.2562389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1690037) q[2];
sx q[2];
rz(-0.4946332) q[2];
sx q[2];
rz(-0.063610323) q[2];
rz(-2.5566067) q[3];
sx q[3];
rz(-1.7413185) q[3];
sx q[3];
rz(1.6271648) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99120283) q[0];
sx q[0];
rz(-1.7504033) q[0];
sx q[0];
rz(2.4124131) q[0];
rz(1.179262) q[1];
sx q[1];
rz(-2.3919892) q[1];
sx q[1];
rz(0.29456219) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8108494) q[0];
sx q[0];
rz(-1.0746351) q[0];
sx q[0];
rz(1.8627573) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0038347) q[2];
sx q[2];
rz(-1.4708418) q[2];
sx q[2];
rz(-1.2765027) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9443317) q[1];
sx q[1];
rz(-2.3880868) q[1];
sx q[1];
rz(1.7427518) q[1];
rz(-pi) q[2];
rz(-0.38334966) q[3];
sx q[3];
rz(-0.93689303) q[3];
sx q[3];
rz(1.0312361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4906759) q[2];
sx q[2];
rz(-1.4399521) q[2];
sx q[2];
rz(-0.73224625) q[2];
rz(-2.2722774) q[3];
sx q[3];
rz(-1.1704051) q[3];
sx q[3];
rz(1.4066345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9100087) q[0];
sx q[0];
rz(-1.5861479) q[0];
sx q[0];
rz(-2.3439132) q[0];
rz(2.8186467) q[1];
sx q[1];
rz(-0.99635092) q[1];
sx q[1];
rz(-1.7332044) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3638316) q[0];
sx q[0];
rz(-2.1182502) q[0];
sx q[0];
rz(1.8108032) q[0];
x q[1];
rz(2.9041821) q[2];
sx q[2];
rz(-2.0185197) q[2];
sx q[2];
rz(0.46502772) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3288142) q[1];
sx q[1];
rz(-1.3359038) q[1];
sx q[1];
rz(2.3702904) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8208169) q[3];
sx q[3];
rz(-1.7446221) q[3];
sx q[3];
rz(2.6717693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.8352167) q[2];
sx q[2];
rz(-2.0312467) q[2];
sx q[2];
rz(2.9442673) q[2];
rz(-0.80162445) q[3];
sx q[3];
rz(-2.1620731) q[3];
sx q[3];
rz(-0.42158034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1590969) q[0];
sx q[0];
rz(-1.6668586) q[0];
sx q[0];
rz(0.91484797) q[0];
rz(-2.9313056) q[1];
sx q[1];
rz(-2.5631914) q[1];
sx q[1];
rz(0.68967825) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1062968) q[0];
sx q[0];
rz(-1.6637422) q[0];
sx q[0];
rz(0.22648099) q[0];
x q[1];
rz(3.1269011) q[2];
sx q[2];
rz(-1.6882875) q[2];
sx q[2];
rz(1.9368287) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8685329) q[1];
sx q[1];
rz(-1.0976205) q[1];
sx q[1];
rz(2.563963) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1237064) q[3];
sx q[3];
rz(-0.59846557) q[3];
sx q[3];
rz(-2.726647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0939533) q[2];
sx q[2];
rz(-2.2932105) q[2];
sx q[2];
rz(-2.8116255) q[2];
rz(-1.0432358) q[3];
sx q[3];
rz(-1.4179351) q[3];
sx q[3];
rz(-0.51658336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.030815) q[0];
sx q[0];
rz(-1.3309706) q[0];
sx q[0];
rz(-1.0571085) q[0];
rz(0.2541751) q[1];
sx q[1];
rz(-1.1421685) q[1];
sx q[1];
rz(2.4370297) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.097226133) q[0];
sx q[0];
rz(-1.2597879) q[0];
sx q[0];
rz(1.0754271) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.754934) q[2];
sx q[2];
rz(-1.6312851) q[2];
sx q[2];
rz(-2.9764701) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5479991) q[1];
sx q[1];
rz(-1.3915466) q[1];
sx q[1];
rz(-2.0229193) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7514486) q[3];
sx q[3];
rz(-1.3657939) q[3];
sx q[3];
rz(-1.6222749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5083984) q[2];
sx q[2];
rz(-1.1722379) q[2];
sx q[2];
rz(-2.634826) q[2];
rz(0.1861598) q[3];
sx q[3];
rz(-0.23701826) q[3];
sx q[3];
rz(2.6059634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7645466) q[0];
sx q[0];
rz(-1.6896387) q[0];
sx q[0];
rz(1.1402546) q[0];
rz(2.2346732) q[1];
sx q[1];
rz(-0.74105558) q[1];
sx q[1];
rz(1.8190609) q[1];
rz(-2.604031) q[2];
sx q[2];
rz(-2.1765709) q[2];
sx q[2];
rz(-1.4449262) q[2];
rz(0.32313866) q[3];
sx q[3];
rz(-1.3165786) q[3];
sx q[3];
rz(1.0490882) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
