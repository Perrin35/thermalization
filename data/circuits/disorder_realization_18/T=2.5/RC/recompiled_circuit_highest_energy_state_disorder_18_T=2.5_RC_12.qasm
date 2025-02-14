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
rz(-1.9189605) q[0];
sx q[0];
rz(-0.55813342) q[0];
sx q[0];
rz(0.65879917) q[0];
rz(0.88762033) q[1];
sx q[1];
rz(2.4467111) q[1];
sx q[1];
rz(12.477787) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50669955) q[0];
sx q[0];
rz(-1.1453712) q[0];
sx q[0];
rz(0.21197196) q[0];
rz(-2.7062662) q[2];
sx q[2];
rz(-1.0485888) q[2];
sx q[2];
rz(-2.2919523) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0483612) q[1];
sx q[1];
rz(-2.1856896) q[1];
sx q[1];
rz(-0.1860861) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1351571) q[3];
sx q[3];
rz(-2.2794323) q[3];
sx q[3];
rz(-2.9624727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.52458557) q[2];
sx q[2];
rz(-1.969939) q[2];
sx q[2];
rz(-0.48801547) q[2];
rz(-2.6739142) q[3];
sx q[3];
rz(-0.52507639) q[3];
sx q[3];
rz(0.16744965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2117598) q[0];
sx q[0];
rz(-0.81120315) q[0];
sx q[0];
rz(1.2540586) q[0];
rz(-2.5414741) q[1];
sx q[1];
rz(-1.3328726) q[1];
sx q[1];
rz(-0.41807237) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1644884) q[0];
sx q[0];
rz(-0.56796861) q[0];
sx q[0];
rz(-0.33645679) q[0];
x q[1];
rz(-0.16440161) q[2];
sx q[2];
rz(-0.65509701) q[2];
sx q[2];
rz(3.0384105) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.1723056) q[1];
sx q[1];
rz(-1.2551771) q[1];
sx q[1];
rz(3.0728042) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3206059) q[3];
sx q[3];
rz(-2.0625522) q[3];
sx q[3];
rz(3.0420835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.40689251) q[2];
sx q[2];
rz(-2.0105346) q[2];
sx q[2];
rz(3.0231754) q[2];
rz(0.40220574) q[3];
sx q[3];
rz(-1.4879358) q[3];
sx q[3];
rz(-1.3291298) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44620946) q[0];
sx q[0];
rz(-2.8948247) q[0];
sx q[0];
rz(-1.1545908) q[0];
rz(-2.0788976) q[1];
sx q[1];
rz(-2.1176391) q[1];
sx q[1];
rz(-1.1850932) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4615153) q[0];
sx q[0];
rz(-2.3993368) q[0];
sx q[0];
rz(-1.1552325) q[0];
x q[1];
rz(1.1923033) q[2];
sx q[2];
rz(-0.54735294) q[2];
sx q[2];
rz(-2.0277835) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.628988) q[1];
sx q[1];
rz(-1.1956685) q[1];
sx q[1];
rz(-0.30156044) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8293158) q[3];
sx q[3];
rz(-1.1791157) q[3];
sx q[3];
rz(1.2092115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2350754) q[2];
sx q[2];
rz(-2.2791028) q[2];
sx q[2];
rz(1.4441351) q[2];
rz(0.92153543) q[3];
sx q[3];
rz(-2.5369365) q[3];
sx q[3];
rz(3.0864033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7017355) q[0];
sx q[0];
rz(-3.0210962) q[0];
sx q[0];
rz(-3.0635656) q[0];
rz(-1.1174508) q[1];
sx q[1];
rz(-1.2470587) q[1];
sx q[1];
rz(1.4720346) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8142151) q[0];
sx q[0];
rz(-0.80670415) q[0];
sx q[0];
rz(-0.14847688) q[0];
rz(-pi) q[1];
rz(-2.1536768) q[2];
sx q[2];
rz(-1.1088409) q[2];
sx q[2];
rz(-2.5910026) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2591272) q[1];
sx q[1];
rz(-1.5992016) q[1];
sx q[1];
rz(1.9309429) q[1];
rz(0.72301282) q[3];
sx q[3];
rz(-1.2870871) q[3];
sx q[3];
rz(2.8151388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.83634079) q[2];
sx q[2];
rz(-2.4938816) q[2];
sx q[2];
rz(2.9913537) q[2];
rz(-0.59208208) q[3];
sx q[3];
rz(-1.838107) q[3];
sx q[3];
rz(-1.9534115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93064654) q[0];
sx q[0];
rz(-2.7879614) q[0];
sx q[0];
rz(0.76552248) q[0];
rz(-2.3402479) q[1];
sx q[1];
rz(-2.0739906) q[1];
sx q[1];
rz(1.9498922) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7758492) q[0];
sx q[0];
rz(-1.8738215) q[0];
sx q[0];
rz(0.55935212) q[0];
rz(2.043388) q[2];
sx q[2];
rz(-2.049438) q[2];
sx q[2];
rz(2.9413073) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0914407) q[1];
sx q[1];
rz(-1.5066083) q[1];
sx q[1];
rz(-2.2540171) q[1];
x q[2];
rz(2.3611154) q[3];
sx q[3];
rz(-0.37364081) q[3];
sx q[3];
rz(2.3821156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8754862) q[2];
sx q[2];
rz(-2.2849639) q[2];
sx q[2];
rz(-0.53492707) q[2];
rz(-0.62355012) q[3];
sx q[3];
rz(-2.0303191) q[3];
sx q[3];
rz(0.88304869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9542338) q[0];
sx q[0];
rz(-0.41340241) q[0];
sx q[0];
rz(2.2075388) q[0];
rz(-1.8961228) q[1];
sx q[1];
rz(-1.5555236) q[1];
sx q[1];
rz(-0.62560558) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0348822) q[0];
sx q[0];
rz(-0.83256522) q[0];
sx q[0];
rz(-2.9009079) q[0];
x q[1];
rz(1.4532188) q[2];
sx q[2];
rz(-2.46799) q[2];
sx q[2];
rz(-2.6308311) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6954591) q[1];
sx q[1];
rz(-1.2810105) q[1];
sx q[1];
rz(0.70420806) q[1];
rz(-pi) q[2];
rz(0.91208338) q[3];
sx q[3];
rz(-0.28277031) q[3];
sx q[3];
rz(1.8709489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.049456747) q[2];
sx q[2];
rz(-2.226604) q[2];
sx q[2];
rz(0.12506872) q[2];
rz(0.72758979) q[3];
sx q[3];
rz(-1.5672507) q[3];
sx q[3];
rz(-2.6653813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26559386) q[0];
sx q[0];
rz(-0.47530526) q[0];
sx q[0];
rz(-1.884961) q[0];
rz(-1.9427293) q[1];
sx q[1];
rz(-1.7190944) q[1];
sx q[1];
rz(1.8667603) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.066028) q[0];
sx q[0];
rz(-2.4504553) q[0];
sx q[0];
rz(2.4253571) q[0];
rz(0.9918757) q[2];
sx q[2];
rz(-2.5235368) q[2];
sx q[2];
rz(1.8628054) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.82060888) q[1];
sx q[1];
rz(-0.57772355) q[1];
sx q[1];
rz(0.89022762) q[1];
rz(-pi) q[2];
rz(1.9406464) q[3];
sx q[3];
rz(-0.33888926) q[3];
sx q[3];
rz(1.3047993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.690543) q[2];
sx q[2];
rz(-1.0131016) q[2];
sx q[2];
rz(2.1742382) q[2];
rz(1.5585772) q[3];
sx q[3];
rz(-1.5019006) q[3];
sx q[3];
rz(0.30174747) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3485182) q[0];
sx q[0];
rz(-2.5418042) q[0];
sx q[0];
rz(0.10979688) q[0];
rz(1.9684567) q[1];
sx q[1];
rz(-0.99183142) q[1];
sx q[1];
rz(1.0822302) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.968348) q[0];
sx q[0];
rz(-2.3254407) q[0];
sx q[0];
rz(-1.2307172) q[0];
rz(-pi) q[1];
rz(2.0556446) q[2];
sx q[2];
rz(-1.4746801) q[2];
sx q[2];
rz(-0.80730593) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6746862) q[1];
sx q[1];
rz(-1.162858) q[1];
sx q[1];
rz(0.80927421) q[1];
rz(-pi) q[2];
rz(1.1940597) q[3];
sx q[3];
rz(-2.1749745) q[3];
sx q[3];
rz(0.24482152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.75862306) q[2];
sx q[2];
rz(-2.0489645) q[2];
sx q[2];
rz(-0.36337241) q[2];
rz(-2.162497) q[3];
sx q[3];
rz(-0.64371124) q[3];
sx q[3];
rz(0.50160971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8363504) q[0];
sx q[0];
rz(-0.79638052) q[0];
sx q[0];
rz(-2.7889732) q[0];
rz(1.3615707) q[1];
sx q[1];
rz(-1.1905328) q[1];
sx q[1];
rz(0.1415267) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3196526) q[0];
sx q[0];
rz(-1.7403649) q[0];
sx q[0];
rz(-1.5362306) q[0];
rz(-pi) q[1];
rz(2.5472801) q[2];
sx q[2];
rz(-1.9611352) q[2];
sx q[2];
rz(0.28536404) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.80672164) q[1];
sx q[1];
rz(-1.1523243) q[1];
sx q[1];
rz(-1.970402) q[1];
x q[2];
rz(-0.64543313) q[3];
sx q[3];
rz(-1.107599) q[3];
sx q[3];
rz(-2.2768159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3129468) q[2];
sx q[2];
rz(-1.4742278) q[2];
sx q[2];
rz(-2.055577) q[2];
rz(-0.18276754) q[3];
sx q[3];
rz(-1.026231) q[3];
sx q[3];
rz(-0.97203794) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6163841) q[0];
sx q[0];
rz(-1.5807736) q[0];
sx q[0];
rz(-2.4973448) q[0];
rz(0.066990189) q[1];
sx q[1];
rz(-1.7552152) q[1];
sx q[1];
rz(-3.0279874) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8380801) q[0];
sx q[0];
rz(-1.3467731) q[0];
sx q[0];
rz(3.1372848) q[0];
rz(-pi) q[1];
rz(-3.058931) q[2];
sx q[2];
rz(-1.6726521) q[2];
sx q[2];
rz(1.5993303) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.56845957) q[1];
sx q[1];
rz(-1.4541469) q[1];
sx q[1];
rz(-0.17165143) q[1];
rz(-1.9285517) q[3];
sx q[3];
rz(-1.5309257) q[3];
sx q[3];
rz(-1.8893582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1429448) q[2];
sx q[2];
rz(-1.8871658) q[2];
sx q[2];
rz(-0.76510731) q[2];
rz(-0.1782002) q[3];
sx q[3];
rz(-1.5604115) q[3];
sx q[3];
rz(2.3044738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2160303) q[0];
sx q[0];
rz(-2.3860274) q[0];
sx q[0];
rz(-0.26269333) q[0];
rz(-2.8021011) q[1];
sx q[1];
rz(-1.9120293) q[1];
sx q[1];
rz(1.0614352) q[1];
rz(-0.65503623) q[2];
sx q[2];
rz(-1.8476386) q[2];
sx q[2];
rz(2.4551433) q[2];
rz(2.3332023) q[3];
sx q[3];
rz(-0.42338531) q[3];
sx q[3];
rz(0.12082621) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
