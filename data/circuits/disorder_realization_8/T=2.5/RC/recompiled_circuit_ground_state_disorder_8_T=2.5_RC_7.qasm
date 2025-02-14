OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72579223) q[0];
sx q[0];
rz(1.9771201) q[0];
sx q[0];
rz(10.676072) q[0];
rz(-0.30081055) q[1];
sx q[1];
rz(4.858074) q[1];
sx q[1];
rz(9.3088027) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1303367) q[0];
sx q[0];
rz(-2.7948871) q[0];
sx q[0];
rz(-1.1954855) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0865583) q[2];
sx q[2];
rz(-1.7770443) q[2];
sx q[2];
rz(-0.95824403) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8311316) q[1];
sx q[1];
rz(-0.77378213) q[1];
sx q[1];
rz(1.4813862) q[1];
rz(-pi) q[2];
rz(-0.2282397) q[3];
sx q[3];
rz(-2.5872218) q[3];
sx q[3];
rz(-1.7621463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7684324) q[2];
sx q[2];
rz(-0.065040437) q[2];
sx q[2];
rz(-1.4434641) q[2];
rz(-2.4453435) q[3];
sx q[3];
rz(-1.5158451) q[3];
sx q[3];
rz(-2.7878917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.26831216) q[0];
sx q[0];
rz(-0.074957632) q[0];
sx q[0];
rz(2.7798376) q[0];
rz(1.7573645) q[1];
sx q[1];
rz(-0.39159602) q[1];
sx q[1];
rz(2.1143544) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38798328) q[0];
sx q[0];
rz(-1.3902724) q[0];
sx q[0];
rz(-0.87708824) q[0];
rz(-pi) q[1];
rz(2.9774819) q[2];
sx q[2];
rz(-2.3375953) q[2];
sx q[2];
rz(-1.5259334) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9421922) q[1];
sx q[1];
rz(-2.2052599) q[1];
sx q[1];
rz(0.30997194) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6546685) q[3];
sx q[3];
rz(-1.7954651) q[3];
sx q[3];
rz(-0.39473907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0072173) q[2];
sx q[2];
rz(-0.5642887) q[2];
sx q[2];
rz(-2.1103653) q[2];
rz(2.4550896) q[3];
sx q[3];
rz(-1.406823) q[3];
sx q[3];
rz(0.321872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.4903851) q[0];
sx q[0];
rz(-1.729916) q[0];
sx q[0];
rz(1.7899293) q[0];
rz(-1.5216113) q[1];
sx q[1];
rz(-0.23623513) q[1];
sx q[1];
rz(1.1865541) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8299415) q[0];
sx q[0];
rz(-1.9970987) q[0];
sx q[0];
rz(1.9894129) q[0];
rz(2.6858052) q[2];
sx q[2];
rz(-1.6384587) q[2];
sx q[2];
rz(0.6099786) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6938942) q[1];
sx q[1];
rz(-1.5617626) q[1];
sx q[1];
rz(1.5771405) q[1];
rz(-pi) q[2];
x q[2];
rz(1.925674) q[3];
sx q[3];
rz(-2.4198341) q[3];
sx q[3];
rz(2.1510368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.97570193) q[2];
sx q[2];
rz(-1.8791135) q[2];
sx q[2];
rz(0.28124896) q[2];
rz(0.37505546) q[3];
sx q[3];
rz(-2.227759) q[3];
sx q[3];
rz(2.8121172) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3220045) q[0];
sx q[0];
rz(-0.81398886) q[0];
sx q[0];
rz(0.32459146) q[0];
rz(-1.548467) q[1];
sx q[1];
rz(-2.8209768) q[1];
sx q[1];
rz(2.2494907) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88282132) q[0];
sx q[0];
rz(-3.1207001) q[0];
sx q[0];
rz(-0.56751318) q[0];
rz(-pi) q[1];
rz(2.4302654) q[2];
sx q[2];
rz(-1.6454503) q[2];
sx q[2];
rz(1.1911281) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.39246757) q[1];
sx q[1];
rz(-1.3121079) q[1];
sx q[1];
rz(-1.2767747) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5209429) q[3];
sx q[3];
rz(-0.59403235) q[3];
sx q[3];
rz(-2.1332873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0804245) q[2];
sx q[2];
rz(-1.755038) q[2];
sx q[2];
rz(-2.5648153) q[2];
rz(0.12442496) q[3];
sx q[3];
rz(-2.3915229) q[3];
sx q[3];
rz(0.09593825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6782003) q[0];
sx q[0];
rz(-2.7233349) q[0];
sx q[0];
rz(-0.28045714) q[0];
rz(2.8434143) q[1];
sx q[1];
rz(-3.0715946) q[1];
sx q[1];
rz(-2.9895474) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9579118) q[0];
sx q[0];
rz(-0.059558161) q[0];
sx q[0];
rz(1.0899215) q[0];
rz(-1.399081) q[2];
sx q[2];
rz(-2.6125312) q[2];
sx q[2];
rz(0.79862874) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0714582) q[1];
sx q[1];
rz(-1.9862174) q[1];
sx q[1];
rz(-0.43674119) q[1];
rz(-1.5923813) q[3];
sx q[3];
rz(-2.7587771) q[3];
sx q[3];
rz(-0.60071731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8449479) q[2];
sx q[2];
rz(-1.351202) q[2];
sx q[2];
rz(0.039999261) q[2];
rz(-3.0039039) q[3];
sx q[3];
rz(-0.44860336) q[3];
sx q[3];
rz(-0.78534809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26507759) q[0];
sx q[0];
rz(-3.0410933) q[0];
sx q[0];
rz(-2.5608089) q[0];
rz(2.4526217) q[1];
sx q[1];
rz(-0.13801485) q[1];
sx q[1];
rz(1.8438011) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9286683) q[0];
sx q[0];
rz(-1.9868694) q[0];
sx q[0];
rz(-1.1088158) q[0];
x q[1];
rz(-1.6766729) q[2];
sx q[2];
rz(-2.5953802) q[2];
sx q[2];
rz(0.27413163) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7766469) q[1];
sx q[1];
rz(-0.45393911) q[1];
sx q[1];
rz(-0.52151068) q[1];
x q[2];
rz(2.5291555) q[3];
sx q[3];
rz(-2.4531657) q[3];
sx q[3];
rz(-1.5189309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.40655228) q[2];
sx q[2];
rz(-1.8031969) q[2];
sx q[2];
rz(-1.5945386) q[2];
rz(-2.7100587) q[3];
sx q[3];
rz(-1.1003234) q[3];
sx q[3];
rz(-2.4908454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4887017) q[0];
sx q[0];
rz(-0.015691375) q[0];
sx q[0];
rz(2.5456862) q[0];
rz(1.3504922) q[1];
sx q[1];
rz(-2.6563783) q[1];
sx q[1];
rz(2.5686666) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1112537) q[0];
sx q[0];
rz(-1.5736394) q[0];
sx q[0];
rz(-0.009441998) q[0];
rz(2.7155034) q[2];
sx q[2];
rz(-1.3182148) q[2];
sx q[2];
rz(1.9468228) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.10526722) q[1];
sx q[1];
rz(-0.67720264) q[1];
sx q[1];
rz(0.52097229) q[1];
x q[2];
rz(0.39471976) q[3];
sx q[3];
rz(-0.79666797) q[3];
sx q[3];
rz(1.7834237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4493745) q[2];
sx q[2];
rz(-2.2276679) q[2];
sx q[2];
rz(2.2268028) q[2];
rz(1.5335013) q[3];
sx q[3];
rz(-1.1398818) q[3];
sx q[3];
rz(2.1515414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51946259) q[0];
sx q[0];
rz(-1.7161481) q[0];
sx q[0];
rz(-3.07716) q[0];
rz(2.875115) q[1];
sx q[1];
rz(-0.12424145) q[1];
sx q[1];
rz(1.8775833) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85487932) q[0];
sx q[0];
rz(-0.42511308) q[0];
sx q[0];
rz(0.69010587) q[0];
rz(-pi) q[1];
rz(-1.1830215) q[2];
sx q[2];
rz(-0.90803781) q[2];
sx q[2];
rz(-2.6135664) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.19169584) q[1];
sx q[1];
rz(-2.765871) q[1];
sx q[1];
rz(-0.57447489) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76412075) q[3];
sx q[3];
rz(-0.74148889) q[3];
sx q[3];
rz(-1.6958069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5233351) q[2];
sx q[2];
rz(-1.9674415) q[2];
sx q[2];
rz(2.0972283) q[2];
rz(2.4339645) q[3];
sx q[3];
rz(-0.39053598) q[3];
sx q[3];
rz(1.0467168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3656965) q[0];
sx q[0];
rz(-1.5713659) q[0];
sx q[0];
rz(2.7112992) q[0];
rz(0.73768342) q[1];
sx q[1];
rz(-0.084641181) q[1];
sx q[1];
rz(2.784909) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9707613) q[0];
sx q[0];
rz(-0.30730844) q[0];
sx q[0];
rz(-0.31913646) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8356933) q[2];
sx q[2];
rz(-1.6000433) q[2];
sx q[2];
rz(-2.4550329) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7995924) q[1];
sx q[1];
rz(-2.1260602) q[1];
sx q[1];
rz(-2.0702122) q[1];
rz(-1.6453708) q[3];
sx q[3];
rz(-1.4514203) q[3];
sx q[3];
rz(2.125691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.77013993) q[2];
sx q[2];
rz(-2.3944103) q[2];
sx q[2];
rz(2.746197) q[2];
rz(1.6793518) q[3];
sx q[3];
rz(-1.7489) q[3];
sx q[3];
rz(-0.028954884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8545561) q[0];
sx q[0];
rz(-2.365878) q[0];
sx q[0];
rz(2.4391644) q[0];
rz(2.9632945) q[1];
sx q[1];
rz(-0.65419227) q[1];
sx q[1];
rz(-1.6184695) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0337153) q[0];
sx q[0];
rz(-1.3340257) q[0];
sx q[0];
rz(1.9997755) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4236529) q[2];
sx q[2];
rz(-0.098966397) q[2];
sx q[2];
rz(2.5731125) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0023481) q[1];
sx q[1];
rz(-1.7114471) q[1];
sx q[1];
rz(2.572525) q[1];
rz(-pi) q[2];
rz(0.63796343) q[3];
sx q[3];
rz(-0.59989446) q[3];
sx q[3];
rz(-1.1092345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.496333) q[2];
sx q[2];
rz(-2.1767904) q[2];
sx q[2];
rz(-2.1920835) q[2];
rz(-1.1358787) q[3];
sx q[3];
rz(-0.94616008) q[3];
sx q[3];
rz(-2.5778263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5180494) q[0];
sx q[0];
rz(-1.2508871) q[0];
sx q[0];
rz(2.293806) q[0];
rz(1.4354979) q[1];
sx q[1];
rz(-1.2789627) q[1];
sx q[1];
rz(0.36437558) q[1];
rz(-0.87068229) q[2];
sx q[2];
rz(-1.489594) q[2];
sx q[2];
rz(1.7832047) q[2];
rz(-0.061312231) q[3];
sx q[3];
rz(-1.8777962) q[3];
sx q[3];
rz(2.5971436) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
