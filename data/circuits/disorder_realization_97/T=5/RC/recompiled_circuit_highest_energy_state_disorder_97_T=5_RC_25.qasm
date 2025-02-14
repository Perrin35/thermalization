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
rz(1.3069557) q[0];
sx q[0];
rz(4.0210273) q[0];
sx q[0];
rz(9.9183912) q[0];
rz(-1.2797132) q[1];
sx q[1];
rz(3.7682025) q[1];
sx q[1];
rz(11.087853) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0042838) q[0];
sx q[0];
rz(-2.2861436) q[0];
sx q[0];
rz(-0.8887995) q[0];
rz(-pi) q[1];
rz(-1.2463208) q[2];
sx q[2];
rz(-0.68973422) q[2];
sx q[2];
rz(-1.5943499) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.263674) q[1];
sx q[1];
rz(-1.3527737) q[1];
sx q[1];
rz(2.430357) q[1];
x q[2];
rz(-0.36840393) q[3];
sx q[3];
rz(-0.96782902) q[3];
sx q[3];
rz(-1.6507698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0066068) q[2];
sx q[2];
rz(-1.1400433) q[2];
sx q[2];
rz(2.0331649) q[2];
rz(-1.8290352) q[3];
sx q[3];
rz(-2.8969942) q[3];
sx q[3];
rz(-0.31497064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7600064) q[0];
sx q[0];
rz(-0.8257603) q[0];
sx q[0];
rz(-0.015856892) q[0];
rz(0.99041692) q[1];
sx q[1];
rz(-0.76737338) q[1];
sx q[1];
rz(2.2784065) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6619157) q[0];
sx q[0];
rz(-2.4004395) q[0];
sx q[0];
rz(-2.9984498) q[0];
x q[1];
rz(0.77570373) q[2];
sx q[2];
rz(-0.45368567) q[2];
sx q[2];
rz(2.016748) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9941123) q[1];
sx q[1];
rz(-1.0801472) q[1];
sx q[1];
rz(1.7685686) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7402418) q[3];
sx q[3];
rz(-0.72899216) q[3];
sx q[3];
rz(2.9251826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.70011675) q[2];
sx q[2];
rz(-0.2773383) q[2];
sx q[2];
rz(0.46879834) q[2];
rz(2.4539454) q[3];
sx q[3];
rz(-1.5117437) q[3];
sx q[3];
rz(0.30011737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0013292) q[0];
sx q[0];
rz(-2.2183473) q[0];
sx q[0];
rz(-0.73626751) q[0];
rz(-2.7983792) q[1];
sx q[1];
rz(-0.88756573) q[1];
sx q[1];
rz(-0.18377486) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7385173) q[0];
sx q[0];
rz(-2.3949497) q[0];
sx q[0];
rz(-3.0439506) q[0];
rz(-pi) q[1];
rz(-2.5415475) q[2];
sx q[2];
rz(-0.79216865) q[2];
sx q[2];
rz(-0.60952696) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0597093) q[1];
sx q[1];
rz(-0.78468152) q[1];
sx q[1];
rz(-2.7126524) q[1];
rz(-2.7815357) q[3];
sx q[3];
rz(-1.0601794) q[3];
sx q[3];
rz(-1.0086446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.78202128) q[2];
sx q[2];
rz(-0.51681334) q[2];
sx q[2];
rz(-0.51592958) q[2];
rz(-2.0333717) q[3];
sx q[3];
rz(-2.5907232) q[3];
sx q[3];
rz(-1.9903323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16294031) q[0];
sx q[0];
rz(-1.2115703) q[0];
sx q[0];
rz(-0.51280713) q[0];
rz(0.45979744) q[1];
sx q[1];
rz(-1.5922981) q[1];
sx q[1];
rz(-0.76382452) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5483185) q[0];
sx q[0];
rz(-1.700907) q[0];
sx q[0];
rz(-1.6862555) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8945208) q[2];
sx q[2];
rz(-1.9678332) q[2];
sx q[2];
rz(2.4415602) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0415548) q[1];
sx q[1];
rz(-0.44719346) q[1];
sx q[1];
rz(2.2639422) q[1];
rz(-pi) q[2];
rz(1.4277589) q[3];
sx q[3];
rz(-1.9120354) q[3];
sx q[3];
rz(-0.18254292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0189455) q[2];
sx q[2];
rz(-0.15572369) q[2];
sx q[2];
rz(2.8141008) q[2];
rz(2.859595) q[3];
sx q[3];
rz(-1.3982541) q[3];
sx q[3];
rz(-1.5871083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.17664385) q[0];
sx q[0];
rz(-2.7556941) q[0];
sx q[0];
rz(-0.97736812) q[0];
rz(-2.7727959) q[1];
sx q[1];
rz(-2.2803523) q[1];
sx q[1];
rz(1.4126973) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4810886) q[0];
sx q[0];
rz(-1.7605283) q[0];
sx q[0];
rz(-1.0718179) q[0];
x q[1];
rz(-0.013864442) q[2];
sx q[2];
rz(-1.0540451) q[2];
sx q[2];
rz(2.2114829) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.32294867) q[1];
sx q[1];
rz(-1.539577) q[1];
sx q[1];
rz(1.853456) q[1];
rz(-pi) q[2];
rz(-2.5630412) q[3];
sx q[3];
rz(-1.2935841) q[3];
sx q[3];
rz(-1.8892242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29350027) q[2];
sx q[2];
rz(-2.3257181) q[2];
sx q[2];
rz(0.14981848) q[2];
rz(-0.39854974) q[3];
sx q[3];
rz(-1.5236676) q[3];
sx q[3];
rz(2.985305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.3715816) q[0];
sx q[0];
rz(-3.0122029) q[0];
sx q[0];
rz(-0.55321252) q[0];
rz(-0.54168701) q[1];
sx q[1];
rz(-1.0463511) q[1];
sx q[1];
rz(-2.6936626) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89066896) q[0];
sx q[0];
rz(-0.4145997) q[0];
sx q[0];
rz(2.7756734) q[0];
rz(2.7121353) q[2];
sx q[2];
rz(-1.1633368) q[2];
sx q[2];
rz(0.40312672) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.74400157) q[1];
sx q[1];
rz(-2.7251456) q[1];
sx q[1];
rz(1.2214425) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6277091) q[3];
sx q[3];
rz(-1.2384733) q[3];
sx q[3];
rz(3.027738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5003659) q[2];
sx q[2];
rz(-1.9772823) q[2];
sx q[2];
rz(0.79916239) q[2];
rz(-0.3847807) q[3];
sx q[3];
rz(-2.81541) q[3];
sx q[3];
rz(-2.1414115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.36050972) q[0];
sx q[0];
rz(-2.0979083) q[0];
sx q[0];
rz(-2.5497896) q[0];
rz(-2.7599755) q[1];
sx q[1];
rz(-0.38002574) q[1];
sx q[1];
rz(-2.9891678) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29983175) q[0];
sx q[0];
rz(-2.526847) q[0];
sx q[0];
rz(2.1257945) q[0];
x q[1];
rz(0.32976361) q[2];
sx q[2];
rz(-2.641039) q[2];
sx q[2];
rz(-2.9422174) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.72327282) q[1];
sx q[1];
rz(-1.5525463) q[1];
sx q[1];
rz(-0.34679522) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2493208) q[3];
sx q[3];
rz(-1.0425869) q[3];
sx q[3];
rz(1.5367374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3064208) q[2];
sx q[2];
rz(-2.1495337) q[2];
sx q[2];
rz(1.6888118) q[2];
rz(-2.6460904) q[3];
sx q[3];
rz(-2.317954) q[3];
sx q[3];
rz(0.56408322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38967663) q[0];
sx q[0];
rz(-3.1041807) q[0];
sx q[0];
rz(-0.16146846) q[0];
rz(0.031919315) q[1];
sx q[1];
rz(-2.4999764) q[1];
sx q[1];
rz(-1.8643103) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.204435) q[0];
sx q[0];
rz(-1.3124518) q[0];
sx q[0];
rz(0.72465557) q[0];
x q[1];
rz(-0.23060282) q[2];
sx q[2];
rz(-1.1586894) q[2];
sx q[2];
rz(1.7526363) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1254221) q[1];
sx q[1];
rz(-1.4628264) q[1];
sx q[1];
rz(-0.30534621) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2767918) q[3];
sx q[3];
rz(-1.6075896) q[3];
sx q[3];
rz(-2.5599673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6892467) q[2];
sx q[2];
rz(-2.2301058) q[2];
sx q[2];
rz(-2.7455043) q[2];
rz(-2.7416157) q[3];
sx q[3];
rz(-0.59498274) q[3];
sx q[3];
rz(-2.2649435) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9006627) q[0];
sx q[0];
rz(-3.1234968) q[0];
sx q[0];
rz(2.990429) q[0];
rz(2.1647272) q[1];
sx q[1];
rz(-0.38115373) q[1];
sx q[1];
rz(0.74473286) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.789294) q[0];
sx q[0];
rz(-1.2937102) q[0];
sx q[0];
rz(0.21488551) q[0];
rz(-pi) q[1];
rz(-2.9989002) q[2];
sx q[2];
rz(-1.9798793) q[2];
sx q[2];
rz(-0.2689223) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7691466) q[1];
sx q[1];
rz(-0.92248864) q[1];
sx q[1];
rz(-1.5689956) q[1];
rz(1.8119393) q[3];
sx q[3];
rz(-1.9134812) q[3];
sx q[3];
rz(1.2872651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.62290827) q[2];
sx q[2];
rz(-1.238995) q[2];
sx q[2];
rz(-0.94804478) q[2];
rz(-1.0575804) q[3];
sx q[3];
rz(-1.6653929) q[3];
sx q[3];
rz(2.0675366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093216151) q[0];
sx q[0];
rz(-2.8792448) q[0];
sx q[0];
rz(2.9076305) q[0];
rz(-1.1853064) q[1];
sx q[1];
rz(-2.2133841) q[1];
sx q[1];
rz(1.8918461) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7980745) q[0];
sx q[0];
rz(-1.8600704) q[0];
sx q[0];
rz(-0.57855655) q[0];
rz(-0.62057497) q[2];
sx q[2];
rz(-1.4016376) q[2];
sx q[2];
rz(0.38693869) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0873336) q[1];
sx q[1];
rz(-2.1076084) q[1];
sx q[1];
rz(2.3786503) q[1];
rz(-pi) q[2];
rz(-0.21548157) q[3];
sx q[3];
rz(-2.4487447) q[3];
sx q[3];
rz(-0.97164916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.027111) q[2];
sx q[2];
rz(-1.1610843) q[2];
sx q[2];
rz(1.7005881) q[2];
rz(0.13127413) q[3];
sx q[3];
rz(-0.461853) q[3];
sx q[3];
rz(-0.24391267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.092939) q[0];
sx q[0];
rz(-0.96554148) q[0];
sx q[0];
rz(1.1337793) q[0];
rz(-1.0406021) q[1];
sx q[1];
rz(-2.2941209) q[1];
sx q[1];
rz(2.6170731) q[1];
rz(-1.16083) q[2];
sx q[2];
rz(-0.66902918) q[2];
sx q[2];
rz(1.4794028) q[2];
rz(-0.20411251) q[3];
sx q[3];
rz(-0.64328803) q[3];
sx q[3];
rz(-2.358123) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
