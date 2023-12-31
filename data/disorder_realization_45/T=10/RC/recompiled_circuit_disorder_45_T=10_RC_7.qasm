OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.835445) q[0];
sx q[0];
rz(-0.68343502) q[0];
sx q[0];
rz(0.47877065) q[0];
rz(0.03102826) q[1];
sx q[1];
rz(5.1217084) q[1];
sx q[1];
rz(6.9245467) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3135932) q[0];
sx q[0];
rz(-1.8090973) q[0];
sx q[0];
rz(0.39512623) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2864561) q[2];
sx q[2];
rz(-0.54358608) q[2];
sx q[2];
rz(-2.2216575) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3798843) q[1];
sx q[1];
rz(-1.1176846) q[1];
sx q[1];
rz(-2.3945827) q[1];
rz(-pi) q[2];
rz(2.6061329) q[3];
sx q[3];
rz(-0.70264953) q[3];
sx q[3];
rz(-3.1010166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.550094) q[2];
sx q[2];
rz(-1.9248362) q[2];
sx q[2];
rz(-2.6634898) q[2];
rz(-1.452662) q[3];
sx q[3];
rz(-2.0958459) q[3];
sx q[3];
rz(-2.5959192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96145445) q[0];
sx q[0];
rz(-2.366876) q[0];
sx q[0];
rz(1.0189198) q[0];
rz(1.6628751) q[1];
sx q[1];
rz(-2.5264085) q[1];
sx q[1];
rz(0.63308024) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4396673) q[0];
sx q[0];
rz(-0.81322008) q[0];
sx q[0];
rz(-1.4600091) q[0];
x q[1];
rz(-0.74866809) q[2];
sx q[2];
rz(-2.0249172) q[2];
sx q[2];
rz(1.5926966) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1770597) q[1];
sx q[1];
rz(-0.53828439) q[1];
sx q[1];
rz(-0.86175491) q[1];
rz(-0.31140621) q[3];
sx q[3];
rz(-0.91324556) q[3];
sx q[3];
rz(-0.94535512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7818266) q[2];
sx q[2];
rz(-0.40010139) q[2];
sx q[2];
rz(3.0241372) q[2];
rz(0.30101267) q[3];
sx q[3];
rz(-1.3344701) q[3];
sx q[3];
rz(1.7787748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85787073) q[0];
sx q[0];
rz(-0.71325934) q[0];
sx q[0];
rz(0.088407956) q[0];
rz(-1.1075426) q[1];
sx q[1];
rz(-2.2231367) q[1];
sx q[1];
rz(-2.450768) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.847825) q[0];
sx q[0];
rz(-0.17058897) q[0];
sx q[0];
rz(2.6831496) q[0];
rz(-pi) q[1];
rz(-2.8969953) q[2];
sx q[2];
rz(-0.78002143) q[2];
sx q[2];
rz(-0.56842677) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1289127) q[1];
sx q[1];
rz(-1.7528273) q[1];
sx q[1];
rz(-0.026489594) q[1];
rz(2.7615943) q[3];
sx q[3];
rz(-1.2424386) q[3];
sx q[3];
rz(-0.85193714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.92418015) q[2];
sx q[2];
rz(-1.2574544) q[2];
sx q[2];
rz(-3.0351191) q[2];
rz(1.7051833) q[3];
sx q[3];
rz(-0.36542106) q[3];
sx q[3];
rz(1.6586554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0687662) q[0];
sx q[0];
rz(-0.23385736) q[0];
sx q[0];
rz(2.6191214) q[0];
rz(2.8126295) q[1];
sx q[1];
rz(-1.6429699) q[1];
sx q[1];
rz(2.5879588) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75008167) q[0];
sx q[0];
rz(-0.9010074) q[0];
sx q[0];
rz(-2.1704587) q[0];
rz(-pi) q[1];
x q[1];
rz(2.373898) q[2];
sx q[2];
rz(-0.97937095) q[2];
sx q[2];
rz(0.64085811) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5428634) q[1];
sx q[1];
rz(-0.79303629) q[1];
sx q[1];
rz(1.383177) q[1];
rz(-pi) q[2];
rz(-1.5997361) q[3];
sx q[3];
rz(-0.77416285) q[3];
sx q[3];
rz(2.418747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.557495) q[2];
sx q[2];
rz(-0.9157052) q[2];
sx q[2];
rz(0.49475691) q[2];
rz(2.2385712) q[3];
sx q[3];
rz(-1.5938063) q[3];
sx q[3];
rz(1.8516699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6506127) q[0];
sx q[0];
rz(-2.3968093) q[0];
sx q[0];
rz(-1.0850798) q[0];
rz(-0.96013534) q[1];
sx q[1];
rz(-1.9624058) q[1];
sx q[1];
rz(2.9575612) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3071614) q[0];
sx q[0];
rz(-0.30858985) q[0];
sx q[0];
rz(-0.35085268) q[0];
rz(-2.7765772) q[2];
sx q[2];
rz(-1.3114309) q[2];
sx q[2];
rz(2.5600381) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9102238) q[1];
sx q[1];
rz(-2.0407045) q[1];
sx q[1];
rz(-2.4889357) q[1];
x q[2];
rz(1.7720513) q[3];
sx q[3];
rz(-0.85490899) q[3];
sx q[3];
rz(1.9191238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.90157834) q[2];
sx q[2];
rz(-2.7928536) q[2];
sx q[2];
rz(-1.1408172) q[2];
rz(-2.7187637) q[3];
sx q[3];
rz(-1.6641649) q[3];
sx q[3];
rz(1.1269349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35247701) q[0];
sx q[0];
rz(-1.1081835) q[0];
sx q[0];
rz(-0.88678962) q[0];
rz(2.9011762) q[1];
sx q[1];
rz(-2.1227032) q[1];
sx q[1];
rz(-2.9930847) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1214971) q[0];
sx q[0];
rz(-1.4098865) q[0];
sx q[0];
rz(0.73385977) q[0];
x q[1];
rz(-1.6806904) q[2];
sx q[2];
rz(-1.3434778) q[2];
sx q[2];
rz(2.0109039) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.099523274) q[1];
sx q[1];
rz(-1.6798786) q[1];
sx q[1];
rz(-1.8314349) q[1];
rz(-2.4466483) q[3];
sx q[3];
rz(-2.5616025) q[3];
sx q[3];
rz(0.65141962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.85577661) q[2];
sx q[2];
rz(-2.6591876) q[2];
sx q[2];
rz(2.6532069) q[2];
rz(0.48940247) q[3];
sx q[3];
rz(-0.92178744) q[3];
sx q[3];
rz(1.2667806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69328904) q[0];
sx q[0];
rz(-1.332809) q[0];
sx q[0];
rz(0.81800246) q[0];
rz(1.8891634) q[1];
sx q[1];
rz(-2.7493582) q[1];
sx q[1];
rz(-2.0163527) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.274652) q[0];
sx q[0];
rz(-0.38892239) q[0];
sx q[0];
rz(-1.2961943) q[0];
rz(-1.2721887) q[2];
sx q[2];
rz(-1.4223137) q[2];
sx q[2];
rz(-0.4543002) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6553073) q[1];
sx q[1];
rz(-0.40973046) q[1];
sx q[1];
rz(-2.3182931) q[1];
x q[2];
rz(2.1358228) q[3];
sx q[3];
rz(-0.98494782) q[3];
sx q[3];
rz(2.1573531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0685048) q[2];
sx q[2];
rz(-2.1606074) q[2];
sx q[2];
rz(-2.4970064) q[2];
rz(1.4792431) q[3];
sx q[3];
rz(-2.1846266) q[3];
sx q[3];
rz(-1.7361599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.1709764) q[0];
sx q[0];
rz(-3.0727486) q[0];
sx q[0];
rz(-1.53565) q[0];
rz(1.9203141) q[1];
sx q[1];
rz(-1.6284643) q[1];
sx q[1];
rz(2.1309526) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0897652) q[0];
sx q[0];
rz(-2.3065789) q[0];
sx q[0];
rz(-2.9249973) q[0];
rz(-pi) q[1];
rz(2.8433617) q[2];
sx q[2];
rz(-2.0474307) q[2];
sx q[2];
rz(-2.6963866) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6992221) q[1];
sx q[1];
rz(-1.3500299) q[1];
sx q[1];
rz(-1.3352331) q[1];
rz(1.6628077) q[3];
sx q[3];
rz(-1.1717456) q[3];
sx q[3];
rz(-1.240318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5236139) q[2];
sx q[2];
rz(-0.31948677) q[2];
sx q[2];
rz(-2.4475205) q[2];
rz(0.56898919) q[3];
sx q[3];
rz(-1.4469955) q[3];
sx q[3];
rz(0.93769658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(3.0555608) q[0];
sx q[0];
rz(-1.9984351) q[0];
sx q[0];
rz(2.6468497) q[0];
rz(-2.5231979) q[1];
sx q[1];
rz(-1.6463552) q[1];
sx q[1];
rz(-0.075597413) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8537124) q[0];
sx q[0];
rz(-1.259385) q[0];
sx q[0];
rz(2.1561554) q[0];
rz(1.8972626) q[2];
sx q[2];
rz(-0.18507659) q[2];
sx q[2];
rz(1.2125804) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1765715) q[1];
sx q[1];
rz(-2.2551943) q[1];
sx q[1];
rz(-1.4375163) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62065403) q[3];
sx q[3];
rz(-1.7194347) q[3];
sx q[3];
rz(1.312048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.33346924) q[2];
sx q[2];
rz(-1.418891) q[2];
sx q[2];
rz(1.3405651) q[2];
rz(-0.30424413) q[3];
sx q[3];
rz(-2.2777568) q[3];
sx q[3];
rz(-2.5568967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8463523) q[0];
sx q[0];
rz(-1.3051935) q[0];
sx q[0];
rz(0.17679581) q[0];
rz(1.8999752) q[1];
sx q[1];
rz(-0.63443628) q[1];
sx q[1];
rz(0.44100824) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4917883) q[0];
sx q[0];
rz(-1.7735964) q[0];
sx q[0];
rz(-2.9359398) q[0];
rz(1.6454234) q[2];
sx q[2];
rz(-1.7241038) q[2];
sx q[2];
rz(-2.9849433) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2446574) q[1];
sx q[1];
rz(-1.999563) q[1];
sx q[1];
rz(2.7662617) q[1];
x q[2];
rz(1.8701843) q[3];
sx q[3];
rz(-2.1736439) q[3];
sx q[3];
rz(0.30764461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.56197721) q[2];
sx q[2];
rz(-0.56869555) q[2];
sx q[2];
rz(0.30187541) q[2];
rz(-0.89312303) q[3];
sx q[3];
rz(-1.8538657) q[3];
sx q[3];
rz(2.9368403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72398913) q[0];
sx q[0];
rz(-1.2932734) q[0];
sx q[0];
rz(-1.4751157) q[0];
rz(3.1148615) q[1];
sx q[1];
rz(-1.6865128) q[1];
sx q[1];
rz(-1.7105688) q[1];
rz(-2.0886291) q[2];
sx q[2];
rz(-2.3452407) q[2];
sx q[2];
rz(-1.8635545) q[2];
rz(2.5857153) q[3];
sx q[3];
rz(-1.1365969) q[3];
sx q[3];
rz(2.668276) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
