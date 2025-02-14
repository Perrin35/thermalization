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
rz(0.34613553) q[0];
sx q[0];
rz(4.7799568) q[0];
sx q[0];
rz(10.156375) q[0];
rz(-2.6919964) q[1];
sx q[1];
rz(-0.1875339) q[1];
sx q[1];
rz(2.7076758) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62747763) q[0];
sx q[0];
rz(-1.6058996) q[0];
sx q[0];
rz(2.3112464) q[0];
rz(0.97831867) q[2];
sx q[2];
rz(-2.0643532) q[2];
sx q[2];
rz(-1.6101642) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9286926) q[1];
sx q[1];
rz(-2.3688931) q[1];
sx q[1];
rz(2.1479285) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7308957) q[3];
sx q[3];
rz(-0.19354469) q[3];
sx q[3];
rz(-1.5610454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.71243858) q[2];
sx q[2];
rz(-1.6753847) q[2];
sx q[2];
rz(1.0038556) q[2];
rz(0.53576523) q[3];
sx q[3];
rz(-2.2771213) q[3];
sx q[3];
rz(-0.0011477688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4681604) q[0];
sx q[0];
rz(-0.68244451) q[0];
sx q[0];
rz(0.72714192) q[0];
rz(0.06761059) q[1];
sx q[1];
rz(-1.3550242) q[1];
sx q[1];
rz(-2.0179857) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7794117) q[0];
sx q[0];
rz(-1.0955278) q[0];
sx q[0];
rz(1.2527466) q[0];
x q[1];
rz(1.4482165) q[2];
sx q[2];
rz(-1.5059587) q[2];
sx q[2];
rz(-1.4139869) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20642463) q[1];
sx q[1];
rz(-1.4713396) q[1];
sx q[1];
rz(1.2897435) q[1];
rz(-pi) q[2];
rz(2.3109834) q[3];
sx q[3];
rz(-2.5334483) q[3];
sx q[3];
rz(0.45765314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2850538) q[2];
sx q[2];
rz(-1.5679789) q[2];
sx q[2];
rz(2.6598568) q[2];
rz(-0.5528062) q[3];
sx q[3];
rz(-1.0779287) q[3];
sx q[3];
rz(0.089574561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32473096) q[0];
sx q[0];
rz(-0.4751927) q[0];
sx q[0];
rz(-3.0480296) q[0];
rz(1.7612673) q[1];
sx q[1];
rz(-2.1880136) q[1];
sx q[1];
rz(-0.63724744) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1475315) q[0];
sx q[0];
rz(-2.0009181) q[0];
sx q[0];
rz(2.568666) q[0];
x q[1];
rz(0.60814823) q[2];
sx q[2];
rz(-1.7575193) q[2];
sx q[2];
rz(-0.0063467912) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.67575021) q[1];
sx q[1];
rz(-1.9681962) q[1];
sx q[1];
rz(-0.028156412) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8508857) q[3];
sx q[3];
rz(-1.5295431) q[3];
sx q[3];
rz(-2.8009529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20643413) q[2];
sx q[2];
rz(-2.5863402) q[2];
sx q[2];
rz(-0.68378249) q[2];
rz(0.23009662) q[3];
sx q[3];
rz(-1.4068539) q[3];
sx q[3];
rz(-0.42207119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049659599) q[0];
sx q[0];
rz(-0.50685087) q[0];
sx q[0];
rz(1.4403213) q[0];
rz(2.6662042) q[1];
sx q[1];
rz(-2.5071867) q[1];
sx q[1];
rz(-0.58116523) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10705778) q[0];
sx q[0];
rz(-0.85113445) q[0];
sx q[0];
rz(0.11373539) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0330795) q[2];
sx q[2];
rz(-1.3603185) q[2];
sx q[2];
rz(2.0049948) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6254723) q[1];
sx q[1];
rz(-2.5800018) q[1];
sx q[1];
rz(0.59188868) q[1];
x q[2];
rz(-2.851296) q[3];
sx q[3];
rz(-0.50945849) q[3];
sx q[3];
rz(-2.6313033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0541957) q[2];
sx q[2];
rz(-2.9554458) q[2];
sx q[2];
rz(-0.52097121) q[2];
rz(-1.6446796) q[3];
sx q[3];
rz(-1.729634) q[3];
sx q[3];
rz(1.514667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5533376) q[0];
sx q[0];
rz(-0.27565685) q[0];
sx q[0];
rz(2.0113373) q[0];
rz(-2.0626119) q[1];
sx q[1];
rz(-1.9422453) q[1];
sx q[1];
rz(-2.030453) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1905907) q[0];
sx q[0];
rz(-1.4050806) q[0];
sx q[0];
rz(1.1072876) q[0];
rz(-pi) q[1];
rz(-0.52972858) q[2];
sx q[2];
rz(-1.5324943) q[2];
sx q[2];
rz(-1.8977579) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.20566733) q[1];
sx q[1];
rz(-1.7280518) q[1];
sx q[1];
rz(-0.92604154) q[1];
rz(-pi) q[2];
rz(-0.91577282) q[3];
sx q[3];
rz(-1.1202328) q[3];
sx q[3];
rz(0.612606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.40372103) q[2];
sx q[2];
rz(-0.51509866) q[2];
sx q[2];
rz(-2.1007288) q[2];
rz(0.8199842) q[3];
sx q[3];
rz(-1.9085725) q[3];
sx q[3];
rz(-2.5177054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2417004) q[0];
sx q[0];
rz(-0.59881678) q[0];
sx q[0];
rz(-0.65648055) q[0];
rz(1.7680602) q[1];
sx q[1];
rz(-0.7936002) q[1];
sx q[1];
rz(-1.1318644) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0431932) q[0];
sx q[0];
rz(-0.88623673) q[0];
sx q[0];
rz(1.8266023) q[0];
x q[1];
rz(-2.4651338) q[2];
sx q[2];
rz(-2.0148811) q[2];
sx q[2];
rz(1.187834) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.51323503) q[1];
sx q[1];
rz(-1.5984922) q[1];
sx q[1];
rz(2.2818034) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5664212) q[3];
sx q[3];
rz(-1.7212015) q[3];
sx q[3];
rz(0.37585092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6953096) q[2];
sx q[2];
rz(-2.535847) q[2];
sx q[2];
rz(0.31965762) q[2];
rz(2.9122635) q[3];
sx q[3];
rz(-1.0234443) q[3];
sx q[3];
rz(0.91013175) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0730154) q[0];
sx q[0];
rz(-1.1443161) q[0];
sx q[0];
rz(1.9101494) q[0];
rz(3.0629509) q[1];
sx q[1];
rz(-1.6784724) q[1];
sx q[1];
rz(-0.15422779) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4336727) q[0];
sx q[0];
rz(-1.7485973) q[0];
sx q[0];
rz(0.34401293) q[0];
rz(-pi) q[1];
rz(-1.1604105) q[2];
sx q[2];
rz(-1.2912116) q[2];
sx q[2];
rz(-1.2321763) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7423305) q[1];
sx q[1];
rz(-0.62472099) q[1];
sx q[1];
rz(-2.4769251) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6811872) q[3];
sx q[3];
rz(-2.6889927) q[3];
sx q[3];
rz(0.79595882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.33588931) q[2];
sx q[2];
rz(-1.615639) q[2];
sx q[2];
rz(1.6501144) q[2];
rz(-1.7283745) q[3];
sx q[3];
rz(-1.7468942) q[3];
sx q[3];
rz(1.5707877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042628057) q[0];
sx q[0];
rz(-1.093981) q[0];
sx q[0];
rz(1.4549103) q[0];
rz(0.97575724) q[1];
sx q[1];
rz(-2.0819596) q[1];
sx q[1];
rz(1.3386493) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.094269116) q[0];
sx q[0];
rz(-1.1828831) q[0];
sx q[0];
rz(2.3045425) q[0];
rz(-pi) q[1];
rz(2.0744616) q[2];
sx q[2];
rz(-1.2758453) q[2];
sx q[2];
rz(2.7558212) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.4563096) q[1];
sx q[1];
rz(-1.7830308) q[1];
sx q[1];
rz(0.93385277) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1986794) q[3];
sx q[3];
rz(-2.2001079) q[3];
sx q[3];
rz(3.1222285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6300388) q[2];
sx q[2];
rz(-2.126882) q[2];
sx q[2];
rz(0.83941984) q[2];
rz(-1.2342341) q[3];
sx q[3];
rz(-1.0890361) q[3];
sx q[3];
rz(-0.031410005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79427528) q[0];
sx q[0];
rz(-1.7592156) q[0];
sx q[0];
rz(2.7224702) q[0];
rz(1.6436815) q[1];
sx q[1];
rz(-1.7828015) q[1];
sx q[1];
rz(-2.7083414) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5113034) q[0];
sx q[0];
rz(-2.0268906) q[0];
sx q[0];
rz(-3.0095741) q[0];
rz(-3.1243043) q[2];
sx q[2];
rz(-0.65708905) q[2];
sx q[2];
rz(1.8985871) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2280088) q[1];
sx q[1];
rz(-2.5435547) q[1];
sx q[1];
rz(2.0889455) q[1];
x q[2];
rz(-0.57859666) q[3];
sx q[3];
rz(-1.75226) q[3];
sx q[3];
rz(-2.6058448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.35932943) q[2];
sx q[2];
rz(-1.9908345) q[2];
sx q[2];
rz(-2.9456054) q[2];
rz(-2.0187142) q[3];
sx q[3];
rz(-2.4298318) q[3];
sx q[3];
rz(1.6897197) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6561683) q[0];
sx q[0];
rz(-2.4805785) q[0];
sx q[0];
rz(1.0957023) q[0];
rz(-0.51086673) q[1];
sx q[1];
rz(-2.8725862) q[1];
sx q[1];
rz(0.21024545) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73186648) q[0];
sx q[0];
rz(-1.7667701) q[0];
sx q[0];
rz(-2.1541697) q[0];
rz(2.2256081) q[2];
sx q[2];
rz(-1.4794) q[2];
sx q[2];
rz(2.2740728) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.89344674) q[1];
sx q[1];
rz(-2.4706744) q[1];
sx q[1];
rz(-2.4144961) q[1];
x q[2];
rz(2.5417762) q[3];
sx q[3];
rz(-1.9164403) q[3];
sx q[3];
rz(-1.4100072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.11253396) q[2];
sx q[2];
rz(-1.6363279) q[2];
sx q[2];
rz(-2.0564334) q[2];
rz(2.8339913) q[3];
sx q[3];
rz(-2.8234973) q[3];
sx q[3];
rz(-2.4881081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6899684) q[0];
sx q[0];
rz(-1.987048) q[0];
sx q[0];
rz(-1.2183627) q[0];
rz(-0.65144173) q[1];
sx q[1];
rz(-0.94964288) q[1];
sx q[1];
rz(-2.3655187) q[1];
rz(3.1262911) q[2];
sx q[2];
rz(-0.87991426) q[2];
sx q[2];
rz(2.8755398) q[2];
rz(0.92523706) q[3];
sx q[3];
rz(-2.400077) q[3];
sx q[3];
rz(-2.526068) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
