OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93249455) q[0];
sx q[0];
rz(-0.068109186) q[0];
sx q[0];
rz(0.77686247) q[0];
rz(-4.4486899) q[1];
sx q[1];
rz(-2.1721462) q[1];
sx q[1];
rz(8.1028508) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8213733) q[0];
sx q[0];
rz(-2.2801274) q[0];
sx q[0];
rz(-0.18289645) q[0];
rz(-0.58682449) q[2];
sx q[2];
rz(-2.5648562) q[2];
sx q[2];
rz(-1.5962034) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3275431) q[1];
sx q[1];
rz(-0.40990489) q[1];
sx q[1];
rz(1.4973212) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9453927) q[3];
sx q[3];
rz(-1.8561188) q[3];
sx q[3];
rz(-1.5232435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0778568) q[2];
sx q[2];
rz(-1.7813762) q[2];
sx q[2];
rz(0.35749164) q[2];
rz(-2.5743971) q[3];
sx q[3];
rz(-1.7287858) q[3];
sx q[3];
rz(-2.7020057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1883063) q[0];
sx q[0];
rz(-0.28306857) q[0];
sx q[0];
rz(-1.8081007) q[0];
rz(-2.7045344) q[1];
sx q[1];
rz(-2.5407365) q[1];
sx q[1];
rz(0.83998799) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14535689) q[0];
sx q[0];
rz(-2.0780434) q[0];
sx q[0];
rz(-2.9661687) q[0];
x q[1];
rz(-2.8905021) q[2];
sx q[2];
rz(-0.60949627) q[2];
sx q[2];
rz(-1.239038) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8312096) q[1];
sx q[1];
rz(-2.134517) q[1];
sx q[1];
rz(0.68117627) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0209818) q[3];
sx q[3];
rz(-1.5375759) q[3];
sx q[3];
rz(-2.0149505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4240894) q[2];
sx q[2];
rz(-1.6987897) q[2];
sx q[2];
rz(1.4060414) q[2];
rz(2.9794335) q[3];
sx q[3];
rz(-1.5806961) q[3];
sx q[3];
rz(2.6774008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22299396) q[0];
sx q[0];
rz(-3.0610237) q[0];
sx q[0];
rz(0.42174569) q[0];
rz(-0.84284198) q[1];
sx q[1];
rz(-1.3858567) q[1];
sx q[1];
rz(0.27580321) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11711794) q[0];
sx q[0];
rz(-1.3399276) q[0];
sx q[0];
rz(-0.33322115) q[0];
rz(-0.47068542) q[2];
sx q[2];
rz(-1.7690725) q[2];
sx q[2];
rz(-0.40775611) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0581857) q[1];
sx q[1];
rz(-2.8166933) q[1];
sx q[1];
rz(-2.9778381) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86902501) q[3];
sx q[3];
rz(-1.4956253) q[3];
sx q[3];
rz(-1.4787256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.80321035) q[2];
sx q[2];
rz(-0.79757491) q[2];
sx q[2];
rz(0.69022834) q[2];
rz(2.7817182) q[3];
sx q[3];
rz(-1.3709603) q[3];
sx q[3];
rz(-0.38351044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2279219) q[0];
sx q[0];
rz(-2.0991195) q[0];
sx q[0];
rz(-2.2699455) q[0];
rz(-2.7323885) q[1];
sx q[1];
rz(-2.6458461) q[1];
sx q[1];
rz(-0.31235487) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.271628) q[0];
sx q[0];
rz(-1.4532491) q[0];
sx q[0];
rz(-2.6948476) q[0];
rz(0.44584195) q[2];
sx q[2];
rz(-0.22805691) q[2];
sx q[2];
rz(-0.40846172) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.96906584) q[1];
sx q[1];
rz(-1.2102264) q[1];
sx q[1];
rz(1.3532624) q[1];
x q[2];
rz(0.64849017) q[3];
sx q[3];
rz(-1.4146418) q[3];
sx q[3];
rz(-2.3462669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.006134) q[2];
sx q[2];
rz(-1.6929071) q[2];
sx q[2];
rz(-0.82497605) q[2];
rz(-1.4247591) q[3];
sx q[3];
rz(-2.8202839) q[3];
sx q[3];
rz(-0.43535522) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7216126) q[0];
sx q[0];
rz(-1.5897607) q[0];
sx q[0];
rz(0.87442526) q[0];
rz(-2.8127316) q[1];
sx q[1];
rz(-1.7560274) q[1];
sx q[1];
rz(1.0985451) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80645567) q[0];
sx q[0];
rz(-2.7493434) q[0];
sx q[0];
rz(2.8492938) q[0];
rz(-pi) q[1];
rz(-2.3057865) q[2];
sx q[2];
rz(-1.5975738) q[2];
sx q[2];
rz(-1.4876897) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.85531536) q[1];
sx q[1];
rz(-2.554727) q[1];
sx q[1];
rz(-1.0081968) q[1];
rz(0.27733488) q[3];
sx q[3];
rz(-2.1430052) q[3];
sx q[3];
rz(-0.25857833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.99981368) q[2];
sx q[2];
rz(-1.9580611) q[2];
sx q[2];
rz(0.67406526) q[2];
rz(1.4874124) q[3];
sx q[3];
rz(-1.4451278) q[3];
sx q[3];
rz(-2.8580247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0312626) q[0];
sx q[0];
rz(-2.1188348) q[0];
sx q[0];
rz(-1.7559825) q[0];
rz(2.022187) q[1];
sx q[1];
rz(-1.4675354) q[1];
sx q[1];
rz(-1.3853692) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7840609) q[0];
sx q[0];
rz(-0.94793944) q[0];
sx q[0];
rz(-2.8297367) q[0];
rz(-0.15689705) q[2];
sx q[2];
rz(-2.6378535) q[2];
sx q[2];
rz(0.7157601) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1277571) q[1];
sx q[1];
rz(-2.7633939) q[1];
sx q[1];
rz(2.5950123) q[1];
rz(2.9117111) q[3];
sx q[3];
rz(-2.982989) q[3];
sx q[3];
rz(1.5029217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0783656) q[2];
sx q[2];
rz(-0.93457064) q[2];
sx q[2];
rz(2.7654977) q[2];
rz(2.0070576) q[3];
sx q[3];
rz(-2.7781656) q[3];
sx q[3];
rz(0.80662066) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019808708) q[0];
sx q[0];
rz(-0.60096318) q[0];
sx q[0];
rz(-0.10251775) q[0];
rz(0.58492297) q[1];
sx q[1];
rz(-2.1939317) q[1];
sx q[1];
rz(1.1835416) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4485156) q[0];
sx q[0];
rz(-1.7127556) q[0];
sx q[0];
rz(1.4407115) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66254424) q[2];
sx q[2];
rz(-1.3032315) q[2];
sx q[2];
rz(-2.2217158) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.592748) q[1];
sx q[1];
rz(-1.3164489) q[1];
sx q[1];
rz(-2.2266822) q[1];
rz(-0.46681301) q[3];
sx q[3];
rz(-1.8100421) q[3];
sx q[3];
rz(1.6166663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.79503235) q[2];
sx q[2];
rz(-0.9298032) q[2];
sx q[2];
rz(0.19980508) q[2];
rz(2.0810769) q[3];
sx q[3];
rz(-2.6001055) q[3];
sx q[3];
rz(-3.0919891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.878433) q[0];
sx q[0];
rz(-1.4819772) q[0];
sx q[0];
rz(2.8286381) q[0];
rz(-0.9494268) q[1];
sx q[1];
rz(-1.3525617) q[1];
sx q[1];
rz(1.688028) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0003371) q[0];
sx q[0];
rz(-2.9337058) q[0];
sx q[0];
rz(0.45335404) q[0];
rz(-0.61954478) q[2];
sx q[2];
rz(-1.6352374) q[2];
sx q[2];
rz(-3.0321414) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9447228) q[1];
sx q[1];
rz(-0.88366449) q[1];
sx q[1];
rz(-1.2792759) q[1];
rz(-pi) q[2];
x q[2];
rz(0.45519228) q[3];
sx q[3];
rz(-1.8774012) q[3];
sx q[3];
rz(0.85807499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7877385) q[2];
sx q[2];
rz(-0.9757897) q[2];
sx q[2];
rz(1.806951) q[2];
rz(-1.3245964) q[3];
sx q[3];
rz(-2.4748804) q[3];
sx q[3];
rz(-1.2142396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.036309328) q[0];
sx q[0];
rz(-0.61573017) q[0];
sx q[0];
rz(2.4435254) q[0];
rz(-2.2659194) q[1];
sx q[1];
rz(-2.0798637) q[1];
sx q[1];
rz(0.36044136) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74340133) q[0];
sx q[0];
rz(-2.0071021) q[0];
sx q[0];
rz(-2.4453336) q[0];
rz(-0.15881947) q[2];
sx q[2];
rz(-2.7544526) q[2];
sx q[2];
rz(-2.5617245) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4354168) q[1];
sx q[1];
rz(-0.98505322) q[1];
sx q[1];
rz(-2.6800568) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5927173) q[3];
sx q[3];
rz(-1.7143119) q[3];
sx q[3];
rz(1.1855857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.91161072) q[2];
sx q[2];
rz(-1.7155827) q[2];
sx q[2];
rz(-0.81725517) q[2];
rz(-0.83003712) q[3];
sx q[3];
rz(-0.54098141) q[3];
sx q[3];
rz(0.219492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63477409) q[0];
sx q[0];
rz(-2.2932678) q[0];
sx q[0];
rz(-0.032055227) q[0];
rz(1.3196779) q[1];
sx q[1];
rz(-1.3751043) q[1];
sx q[1];
rz(2.4874036) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7830156) q[0];
sx q[0];
rz(-2.9111283) q[0];
sx q[0];
rz(-2.5925267) q[0];
rz(2.6903321) q[2];
sx q[2];
rz(-2.3006744) q[2];
sx q[2];
rz(-0.40715363) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.75958222) q[1];
sx q[1];
rz(-2.0339662) q[1];
sx q[1];
rz(0.5447556) q[1];
rz(-pi) q[2];
rz(1.1919674) q[3];
sx q[3];
rz(-2.4740268) q[3];
sx q[3];
rz(0.51578705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0111982) q[2];
sx q[2];
rz(-1.0426714) q[2];
sx q[2];
rz(-2.1350258) q[2];
rz(2.448163) q[3];
sx q[3];
rz(-1.3207685) q[3];
sx q[3];
rz(1.8600195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7158159) q[0];
sx q[0];
rz(-2.4632813) q[0];
sx q[0];
rz(3.0671469) q[0];
rz(-2.2609932) q[1];
sx q[1];
rz(-2.0279299) q[1];
sx q[1];
rz(0.50150064) q[1];
rz(2.2702552) q[2];
sx q[2];
rz(-1.9819145) q[2];
sx q[2];
rz(-2.869538) q[2];
rz(1.952259) q[3];
sx q[3];
rz(-2.2305924) q[3];
sx q[3];
rz(2.3371405) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
