OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(2.2979484) q[0];
sx q[0];
rz(9.2568682) q[0];
rz(-1.9703938) q[1];
sx q[1];
rz(-0.29532239) q[1];
sx q[1];
rz(-0.056161031) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1355609) q[0];
sx q[0];
rz(-0.52838415) q[0];
sx q[0];
rz(2.0023268) q[0];
rz(-pi) q[1];
rz(0.21284717) q[2];
sx q[2];
rz(-2.2058862) q[2];
sx q[2];
rz(-2.0095306) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.15966283) q[1];
sx q[1];
rz(-0.61239457) q[1];
sx q[1];
rz(-2.0484522) q[1];
rz(-pi) q[2];
rz(1.5264741) q[3];
sx q[3];
rz(-1.830415) q[3];
sx q[3];
rz(1.1322024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.37796676) q[2];
sx q[2];
rz(-2.8597735) q[2];
sx q[2];
rz(-0.4326694) q[2];
rz(1.1928605) q[3];
sx q[3];
rz(-1.9038707) q[3];
sx q[3];
rz(2.7584934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6137961) q[0];
sx q[0];
rz(-2.6531117) q[0];
sx q[0];
rz(1.8288076) q[0];
rz(0.20547543) q[1];
sx q[1];
rz(-2.165129) q[1];
sx q[1];
rz(1.1516494) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8173556) q[0];
sx q[0];
rz(-1.865987) q[0];
sx q[0];
rz(0.8582219) q[0];
rz(-2.4685681) q[2];
sx q[2];
rz(-0.7012127) q[2];
sx q[2];
rz(-2.6075624) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.18560219) q[1];
sx q[1];
rz(-0.62528175) q[1];
sx q[1];
rz(-0.76693265) q[1];
rz(-pi) q[2];
rz(-0.50206708) q[3];
sx q[3];
rz(-1.9126529) q[3];
sx q[3];
rz(1.7934007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0097222086) q[2];
sx q[2];
rz(-1.6513609) q[2];
sx q[2];
rz(2.9197664) q[2];
rz(2.7644073) q[3];
sx q[3];
rz(-2.714034) q[3];
sx q[3];
rz(0.77243531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31056988) q[0];
sx q[0];
rz(-3.0492058) q[0];
sx q[0];
rz(-0.036852766) q[0];
rz(0.82551461) q[1];
sx q[1];
rz(-1.3157536) q[1];
sx q[1];
rz(3.085014) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73233561) q[0];
sx q[0];
rz(-1.1261254) q[0];
sx q[0];
rz(2.1952941) q[0];
rz(1.4372196) q[2];
sx q[2];
rz(-1.8250416) q[2];
sx q[2];
rz(2.996252) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11675662) q[1];
sx q[1];
rz(-0.36326888) q[1];
sx q[1];
rz(2.5519752) q[1];
rz(-pi) q[2];
rz(1.9583086) q[3];
sx q[3];
rz(-2.1790677) q[3];
sx q[3];
rz(-2.0351978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.27292192) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(-0.92612129) q[2];
rz(-0.55666322) q[3];
sx q[3];
rz(-2.848048) q[3];
sx q[3];
rz(-1.042897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91519231) q[0];
sx q[0];
rz(-2.4214348) q[0];
sx q[0];
rz(-0.74209374) q[0];
rz(2.0023951) q[1];
sx q[1];
rz(-0.4793872) q[1];
sx q[1];
rz(-0.46359584) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.050042) q[0];
sx q[0];
rz(-2.3766962) q[0];
sx q[0];
rz(-1.1372304) q[0];
rz(1.5092588) q[2];
sx q[2];
rz(-1.9029641) q[2];
sx q[2];
rz(-2.5663944) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78471781) q[1];
sx q[1];
rz(-0.80662913) q[1];
sx q[1];
rz(-0.23826092) q[1];
rz(-1.0158402) q[3];
sx q[3];
rz(-1.98588) q[3];
sx q[3];
rz(2.2500028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3670369) q[2];
sx q[2];
rz(-2.98428) q[2];
sx q[2];
rz(3.0920933) q[2];
rz(0.1285304) q[3];
sx q[3];
rz(-1.5934207) q[3];
sx q[3];
rz(-3.0226504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0054935) q[0];
sx q[0];
rz(-2.7042784) q[0];
sx q[0];
rz(-2.8438925) q[0];
rz(2.659335) q[1];
sx q[1];
rz(-0.75459701) q[1];
sx q[1];
rz(-0.94435) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34236318) q[0];
sx q[0];
rz(-0.21393299) q[0];
sx q[0];
rz(1.3571204) q[0];
x q[1];
rz(-0.40446754) q[2];
sx q[2];
rz(-0.81768113) q[2];
sx q[2];
rz(2.5596465) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3790834) q[1];
sx q[1];
rz(-3.0804539) q[1];
sx q[1];
rz(-1.5361384) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3145507) q[3];
sx q[3];
rz(-1.1503997) q[3];
sx q[3];
rz(-2.8375569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.258761) q[2];
sx q[2];
rz(-1.0927039) q[2];
sx q[2];
rz(3.0333701) q[2];
rz(0.0023068874) q[3];
sx q[3];
rz(-1.6121515) q[3];
sx q[3];
rz(-2.8172857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43679431) q[0];
sx q[0];
rz(-0.44610766) q[0];
sx q[0];
rz(0.58445245) q[0];
rz(-0.8862409) q[1];
sx q[1];
rz(-2.5247572) q[1];
sx q[1];
rz(-0.054919682) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14558218) q[0];
sx q[0];
rz(-1.3025563) q[0];
sx q[0];
rz(-3.0560188) q[0];
rz(-pi) q[1];
rz(1.8940077) q[2];
sx q[2];
rz(-1.5530469) q[2];
sx q[2];
rz(2.7780967) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.37327787) q[1];
sx q[1];
rz(-1.0265961) q[1];
sx q[1];
rz(2.0139704) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46338007) q[3];
sx q[3];
rz(-2.1043092) q[3];
sx q[3];
rz(2.2802071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.75446689) q[2];
sx q[2];
rz(-3.0209164) q[2];
sx q[2];
rz(-2.1248655) q[2];
rz(-2.5975442) q[3];
sx q[3];
rz(-2.7816911) q[3];
sx q[3];
rz(-1.8090766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6761557) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(0.28453919) q[0];
rz(-0.94447213) q[1];
sx q[1];
rz(-1.9453134) q[1];
sx q[1];
rz(2.231266) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.401424) q[0];
sx q[0];
rz(-1.6356042) q[0];
sx q[0];
rz(-0.054697371) q[0];
rz(1.1548642) q[2];
sx q[2];
rz(-1.8850733) q[2];
sx q[2];
rz(1.8679801) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.42156223) q[1];
sx q[1];
rz(-1.3898464) q[1];
sx q[1];
rz(2.6471495) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63776871) q[3];
sx q[3];
rz(-0.83003269) q[3];
sx q[3];
rz(0.72731599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3900782) q[2];
sx q[2];
rz(-0.063515924) q[2];
sx q[2];
rz(-2.2195623) q[2];
rz(2.5743124) q[3];
sx q[3];
rz(-1.4138979) q[3];
sx q[3];
rz(1.0197619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.89408016) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(0.12938736) q[0];
rz(-0.63240504) q[1];
sx q[1];
rz(-2.1148732) q[1];
sx q[1];
rz(2.8410889) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7104615) q[0];
sx q[0];
rz(-0.80695242) q[0];
sx q[0];
rz(1.0459082) q[0];
rz(-0.95327611) q[2];
sx q[2];
rz(-0.91070181) q[2];
sx q[2];
rz(2.2613139) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6943372) q[1];
sx q[1];
rz(-1.9112504) q[1];
sx q[1];
rz(-1.7561595) q[1];
rz(-pi) q[2];
rz(-1.9037876) q[3];
sx q[3];
rz(-2.3559642) q[3];
sx q[3];
rz(1.1803407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.58632103) q[2];
sx q[2];
rz(-0.99493146) q[2];
sx q[2];
rz(-2.3596181) q[2];
rz(-2.590495) q[3];
sx q[3];
rz(-1.3827773) q[3];
sx q[3];
rz(-2.9836695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5683811) q[0];
sx q[0];
rz(-1.9848354) q[0];
sx q[0];
rz(-0.12776275) q[0];
rz(-2.5993775) q[1];
sx q[1];
rz(-2.1844889) q[1];
sx q[1];
rz(-2.382747) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6700867) q[0];
sx q[0];
rz(-0.42450464) q[0];
sx q[0];
rz(1.612081) q[0];
rz(-pi) q[1];
rz(2.5326469) q[2];
sx q[2];
rz(-1.7296089) q[2];
sx q[2];
rz(-1.9132523) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.85941852) q[1];
sx q[1];
rz(-2.4255883) q[1];
sx q[1];
rz(0.23846682) q[1];
rz(-1.5893448) q[3];
sx q[3];
rz(-0.66569257) q[3];
sx q[3];
rz(-1.9513643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1252497) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(2.6515567) q[2];
rz(1.4222493) q[3];
sx q[3];
rz(-1.9336721) q[3];
sx q[3];
rz(1.0629883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35995099) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(-3.066257) q[0];
rz(-0.8967337) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(0.60992253) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41994914) q[0];
sx q[0];
rz(-0.83166612) q[0];
sx q[0];
rz(1.182611) q[0];
rz(-pi) q[1];
rz(2.5209849) q[2];
sx q[2];
rz(-2.6891516) q[2];
sx q[2];
rz(2.0394182) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8049106) q[1];
sx q[1];
rz(-1.8098127) q[1];
sx q[1];
rz(2.3906624) q[1];
rz(0.742357) q[3];
sx q[3];
rz(-1.9543813) q[3];
sx q[3];
rz(0.63391268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.23218368) q[2];
sx q[2];
rz(-0.82513088) q[2];
sx q[2];
rz(0.71371901) q[2];
rz(-0.37832007) q[3];
sx q[3];
rz(-2.648073) q[3];
sx q[3];
rz(-0.87987125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.338035) q[0];
sx q[0];
rz(-1.9914347) q[0];
sx q[0];
rz(1.5557355) q[0];
rz(2.4907885) q[1];
sx q[1];
rz(-1.4918068) q[1];
sx q[1];
rz(3.0202958) q[1];
rz(1.3748319) q[2];
sx q[2];
rz(-1.6008196) q[2];
sx q[2];
rz(0.65884789) q[2];
rz(1.0999023) q[3];
sx q[3];
rz(-1.3445911) q[3];
sx q[3];
rz(0.54683987) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
